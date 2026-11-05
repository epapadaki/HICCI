"""
Paper-style Discussion of the curated findings, grounded on real PubMed.

Pipeline (retrieval-augmented generation):
  1. Gather the app's own computed findings — the curated ligand-receptor
     network, the busiest hubs, and the functional-enrichment highlights.
  2. Retrieve real PubMed literature for the top curated axes using the app's
     existing entity-aware, cached search (so nothing is invented and results
     work offline after the first fetch).
  3. Ask a configured LLM to write a Discussion that may cite ONLY the
     retrieved references, by their numbers, and must not invent PMIDs.

Needs a Gemini API key, read (as in GeneCurator) from the environment variable
GEMINI_API_KEY, normally set in a .env file in the project root, which git
ignores. The key is never written in the code or entered in the app. Without it
the endpoint reports that it is not configured and writes nothing — the model is
never called. The model is constrained to the supplied abstracts so it cannot
fabricate citations.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import time
from typing import Any
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

from . import assessments, curated, functional_enrichment, project_context, settings

from utils.pubmed_literature import search_pubmed_pair_entity_aware

DEFAULT_MODEL = {
    "gemini": "gemini-2.5-flash",
    "openai": "gpt-4o-mini",
    "anthropic": "claude-3-5-sonnet-20241022",
}

# Preference order when auto-picking an available Gemini model.
_GEMINI_PREF = ["gemini-2.5-flash", "gemini-2.0-flash", "gemini-flash-latest",
                "gemini-2.5-pro", "gemini-2.0-pro", "gemini-1.5-flash", "gemini-pro-latest"]

# ─────────────────────────────── API KEY (.env) ───────────────────────────────
# As in GeneCurator, the key is never written in the code. Whoever runs the
# platform puts it in a .env file in the project root (ignored by git):
#     GEMINI_API_KEY=your-key
# Optionally GEMINI_MODEL=gemini-2.5-flash (blank = auto-pick an available model).
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _load_dotenv() -> None:
    """Load PROJECT_ROOT/.env into the environment (python-dotenv if installed,
    otherwise a minimal KEY=VALUE reader). Existing variables are not overridden."""
    env_file = PROJECT_ROOT / ".env"
    try:
        from dotenv import load_dotenv
        load_dotenv(env_file)
        return
    except ImportError:
        pass
    try:
        for line in env_file.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, _, v = line.partition("=")
            k, v = k.strip(), v.strip().strip('"').strip("'")
            if k and k not in os.environ:
                os.environ[k] = v
    except OSError:
        pass


_load_dotenv()

_CACHE: dict[str, dict[str, Any]] = {}


# ------------------------------------------------------------------- config ---
def _llm_config() -> dict[str, str]:
    api_key = os.environ.get("GEMINI_API_KEY", "").strip()
    model = os.environ.get("GEMINI_MODEL", "").strip()
    return {"provider": "gemini", "api_key": api_key, "model": model}


def status() -> dict[str, Any]:
    cfg = _llm_config()
    return {"available": bool(cfg["api_key"]), "provider": cfg["provider"], "model": cfg["model"]}


# ---------------------------------------------------------------- findings ---
def _experiment_context() -> dict[str, str]:
    try:
        ctx = project_context.load()
    except Exception:
        return {}
    groups = ctx.get("experimental_groups") or []
    if isinstance(groups, list):
        groups = ", ".join(str(g) for g in groups if g)
    return {
        "organism": str(ctx.get("organism") or "").strip(),
        "tissue": str(ctx.get("tissue_or_organ") or "").strip(),
        "disease": str(ctx.get("disease_or_condition") or "").strip(),
        "groups": str(groups or "").strip(),
        "notes": str(ctx.get("notes") or "").strip(),
    }


def _context_terms(ctx: dict[str, str]) -> list[str]:
    """Keywords from the experiment's disease/tissue, to bias reference ranking."""
    terms: list[str] = []
    for key in ("disease", "tissue"):
        for w in re.split(r"[^A-Za-z0-9]+", ctx.get(key, "").lower()):
            if len(w) >= 4 and w not in terms:
                terms.append(w)
    return terms


def _findings(cutoff: float) -> dict[str, Any]:
    # Enrichment of the curated genes against the candidate ligand-receptor
    # background (what the curated subset is specifically enriched for).
    enr = functional_enrichment.functional_enrichment(cutoff)
    rew = curated.lr_rewiring(cutoff)
    hubs = curated.molecular_hubs(cutoff)

    enrich = []
    for src in enr.get("sources", []):
        terms = [t["name"] for t in src.get("after", [])[:6]]
        if terms:
            enrich.append({"library": src["label"], "terms": terms})

    return {
        "context": _experiment_context(),
        "before_n": enr.get("before_n", 0),
        "after_n": enr.get("after_n", 0),
        "n_interactions": rew.get("n_interactions", 0),
        "n_ligands": rew.get("n_ligands", 0),
        "n_receptors": rew.get("n_receptors", 0),
        "pairs": rew.get("pairs", [])[:10],
        "top_ligands": [h["name"] for h in hubs.get("ligands", [])[:6]],
        "top_receptors": [h["name"] for h in hubs.get("receptors", [])[:6]],
        "enrichment": enrich,
    }


# -------------------------------------------------------------- references ---
def _authors_short(authors) -> str:
    if not authors:
        return ""
    first = str(authors[0])
    return first + (" et al." if len(authors) > 1 else "")


def _references(pairs: list[dict[str, Any]], ctx_terms: list[str] | None = None,
                max_pairs: int = 8, per_pair: int = 4, cap: int = 18) -> list[dict[str, Any]]:
    ctx_terms = ctx_terms or []
    cand: list[tuple[int, int, dict[str, Any]]] = []  # (relevance, order, ref)
    seen: set[str] = set()
    order = 0
    for p in pairs[:max_pairs]:
        try:
            res = search_pubmed_pair_entity_aware(p["ligand"], p["receptor"], max_results=per_pair)
        except Exception:
            continue
        for art in (res.get("articles") or [])[:per_pair]:
            pmid = str(art.get("pmid") or "")
            if not pmid or pmid in seen:
                continue
            seen.add(pmid)
            abstract = (art.get("abstract") or "").strip().replace("\n", " ")
            hay = (str(art.get("title") or "") + " " + abstract + " "
                   + " ".join(art.get("mesh_terms") or [])).lower()
            score = sum(1 for t in ctx_terms if t in hay)
            if len(abstract) > 900:
                abstract = abstract[:899] + "…"
            cand.append((score, order, {
                "pmid": pmid,
                "title": art.get("title") or "",
                "authors": _authors_short(art.get("authors")),
                "year": art.get("year"),
                "journal": art.get("journal") or "",
                "url": art.get("pubmed_url") or f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "pair": f'{p["ligand"]}–{p["receptor"]}',
                "abstract": abstract,
            }))
            order += 1
    # Context-relevant articles first (disease/tissue match), then original order.
    # When nothing matches, everything scores 0 and the original order is kept, so
    # biasing never removes references — it only reorders.
    cand.sort(key=lambda x: (-x[0], x[1]))
    refs = []
    for i, (_score, _o, ref) in enumerate(cand[:cap]):
        ref = dict(ref)
        ref["n"] = i + 1
        refs.append(ref)
    return refs


# ----------------------------------------------------------------- prompt ---
_SYSTEM = (
    "You are a scientific writer drafting the Discussion section of a bioinformatics "
    "application paper. The study runs cell-cell communication (ligand-receptor) "
    "inference on single-cell RNA-seq and adds expert human-in-the-loop curation to "
    "remove noise and keep biologically credible interactions.\n\n"
    "Frame the whole Discussion in the specific EXPERIMENT CONTEXT given below "
    "(organism, tissue, disease/condition). Interpret the findings for that system — "
    "name the disease/tissue where relevant — rather than writing generically. Do not "
    "introduce unrelated diseases unless a provided reference is directly relevant.\n\n"
    "This Discussion must cover BOTH sets of results in balance: (a) the curated "
    "ligand-receptor network — the main ligand-receptor axes, the busiest hubs, and how "
    "many interactions were kept — and (b) the functional enrichment of the curated genes, "
    "i.e. which biological processes and pathways the curated subset is over-represented for "
    "relative to the background of all candidate ligand-receptor genes (so generic signalling "
    "terms are controlled for). Weave the two together; do not treat either as an afterthought.\n\n"
    "Write a focused Discussion of 4 paragraphs, about 350-450 words, markdown. Interpret the "
    "curated network and the enriched processes/pathways, say what the curated subset is "
    "specifically enriched for, connect to the literature, and give one honest limitations paragraph "
    "(over-representation reflects annotation, not activity; predictions come from RNA "
    "expression, not proven binding; curation is expert judgement).\n\n"
    "STYLE — write like a person, not a language model:\n"
    "- Use PASSIVE VOICE throughout. Do NOT use the first person: no 'we', 'our', 'our study', "
    "'in this study', 'in this work', 'here we', 'the present study'. State what was done and "
    "found impersonally (e.g. 'the curated set was enriched for…', not 'we found…').\n"
    "- Do NOT use these words/phrases: 'robust framework', 'crucial role', 'plays a key/"
    "crucial role', 'underscores', 'highlights', 'sheds light', 'it is important to note/"
    "acknowledge', 'notably', 'furthermore', 'moreover', 'importantly', 'delve', 'intricate', "
    "'nuanced', 'compelling', 'valuable insights', 'critical insights', 'warrant further "
    "investigation', 'pave the way', 'in conclusion', 'overall'.\n"
    "- No formulaic opener (never start with 'The functional enrichment analysis provided "
    "insights into…'). Open on a concrete result. No wrap-up sentence about future work.\n"
    "- Vary sentence length; avoid making every paragraph a topic sentence followed by three "
    "parallel supporting sentences. Use at most one or two connective adverbs in the whole text.\n"
    "- Be concrete and specific (name the disease/tissue, the genes, the numbers) rather than "
    "general. Cut hedging and filler; state things plainly and only hedge where the data "
    "genuinely require it.\n\n"
    "CITATIONS: you may cite ONLY the numbered references provided, using their numbers in "
    "square brackets, e.g. [3] or [3,7]. Never invent a citation, a PMID, an author or a "
    "finding that is not supported by a provided abstract. If the references do not support "
    "a claim, make the claim without a citation or omit it. Do not add a reference list — it "
    "is rendered separately."
)


def _user_prompt(f: dict[str, Any], refs: list[dict[str, Any]]) -> str:
    lines = []
    c = f.get("context") or {}
    if any(c.get(k) for k in ("organism", "tissue", "disease", "groups", "notes")):
        lines.append("## Experiment context\n")
        if c.get("organism"):
            lines.append(f"Organism: {c['organism']}.")
        if c.get("tissue"):
            lines.append(f"Tissue/organ: {c['tissue']}.")
        if c.get("disease"):
            lines.append(f"Disease/condition: {c['disease']}.")
        if c.get("groups"):
            lines.append(f"Experimental groups: {c['groups']}.")
        if c.get("notes"):
            lines.append(f"Notes: {c['notes']}.")
        lines.append("")
    lines.append("## Functional enrichment to discuss\n")
    lines.append(f"Curated gene set: {f['after_n']} genes kept after review, out of "
                 f"{f['before_n']} candidate ligand-receptor genes.")
    if f["enrichment"]:
        lines.append("\nTerms over-represented in the curated genes, tested against the background of all "
                     "candidate ligand-receptor genes (so terms common to ligands/receptors in general are "
                     "controlled for) — top per library:")
        for e in f["enrichment"]:
            if e["terms"]:
                lines.append(f"  - {e['library']}: {', '.join(e['terms'])}")
    lines.append("\n## Curated interaction network\n")
    lines.append(f"Curated network: {f['n_interactions']} kept ligand-receptor interactions, "
                 f"{f['n_ligands']} ligands to {f['n_receptors']} receptors.")
    if f["top_ligands"] or f["top_receptors"]:
        lines.append(f"Busiest ligands: {', '.join(f['top_ligands'])}. "
                     f"Busiest receptors: {', '.join(f['top_receptors'])}.")
    if f["pairs"]:
        lines.append("Top curated axes (ligand -> receptor, kept count): "
                     + "; ".join(f"{p['ligand']}->{p['receptor']} (n={p.get('count',0)})" for p in f["pairs"]))
    lines.append("\n## Numbered references (cite only these, by number)\n")
    for r in refs:
        head = f"[{r['n']}] {r['authors']} ({r['year']}). {r['title']} {r['journal']}. PMID {r['pmid']}. Pair: {r['pair']}."
        lines.append(head)
        if r["abstract"]:
            lines.append(f"      Abstract: {r['abstract']}")
    lines.append("\nWrite the Discussion now.")
    return "\n".join(lines)


# -------------------------------------------------------------- llm clients ---
def _http_json(url: str, headers: dict[str, str], payload: dict[str, Any], timeout: int = 150) -> dict[str, Any]:
    body = json.dumps(payload).encode("utf-8")
    req = Request(url, data=body, headers={**headers, "Content-Type": "application/json"}, method="POST")
    with urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _http_get_json(url: str, timeout: int = 30) -> dict[str, Any]:
    with urlopen(Request(url, method="GET"), timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _gemini_pick_model(api_key: str, want: str = "") -> str:
    """Return a model the key can actually use for generateContent (model names
    churn, so resolve against the live ListModels instead of trusting a default)."""
    out = _http_get_json(f"https://generativelanguage.googleapis.com/v1beta/models?key={api_key}")
    avail = [m.get("name", "").split("/")[-1] for m in out.get("models", [])
             if "generateContent" in (m.get("supportedGenerationMethods") or [])]
    if want and want in avail:
        return want
    for o in _GEMINI_PREF:
        if o in avail:
            return o
    for n in avail:
        if "flash" in n:
            return n
    return avail[0] if avail else (want or DEFAULT_MODEL["gemini"])


def _call_openai(cfg: dict[str, str], system: str, user: str) -> str:
    out = _http_json(
        "https://api.openai.com/v1/chat/completions",
        {"Authorization": "Bearer " + cfg["api_key"]},
        {
            "model": cfg["model"],
            "messages": [{"role": "system", "content": system}, {"role": "user", "content": user}],
            "temperature": 0.6,
        },
    )
    return out["choices"][0]["message"]["content"].strip()


def _call_anthropic(cfg: dict[str, str], system: str, user: str) -> str:
    out = _http_json(
        "https://api.anthropic.com/v1/messages",
        {"x-api-key": cfg["api_key"], "anthropic-version": "2023-06-01"},
        {
            "model": cfg["model"],
            "max_tokens": 4000,
            "temperature": 0.6,
            "system": system,
            "messages": [{"role": "user", "content": user}],
        },
    )
    parts = out.get("content", [])
    return "".join(b.get("text", "") for b in parts if b.get("type") == "text").strip()


def _gemini_generate(cfg: dict[str, str], model: str, system: str, user: str) -> str:
    url = (f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent"
           f"?key={cfg['api_key']}")
    out = _http_json(url, {}, {
        "systemInstruction": {"parts": [{"text": system}]},
        "contents": [{"role": "user", "parts": [{"text": user}]}],
        "generationConfig": {"temperature": 0.6, "maxOutputTokens": 8192},
    })
    cands = out.get("candidates", [])
    if not cands:
        raise RuntimeError("Gemini returned no text (possibly blocked or quota). " + json.dumps(out)[:200])
    parts = cands[0].get("content", {}).get("parts", [])
    return "".join(p.get("text", "") for p in parts).strip()


def _call_gemini(cfg: dict[str, str], system: str, user: str) -> str:
    # Try the resolved model, then a couple of fallbacks, retrying on transient
    # overload (503) or rate limit (429) with a short backoff. Model names that
    # the key cannot use (404) are skipped to the next candidate.
    primary = cfg["model"] or DEFAULT_MODEL["gemini"]
    candidates = [primary] + [m for m in _GEMINI_PREF if m != primary]
    last_err: Exception | None = None
    for model in candidates[:4]:
        for attempt in range(3):
            try:
                text = _gemini_generate(cfg, model, system, user)
                cfg["model"] = model  # report the model that actually answered
                return text
            except HTTPError as e:
                last_err = e
                if e.code in (503, 429):
                    time.sleep(1.5 * (attempt + 1))
                    continue
                if e.code == 404:
                    break  # this model is not available on the key; try the next
                raise
            except URLError:
                raise
    if last_err is not None:
        raise last_err
    raise RuntimeError("Gemini is temporarily unavailable; please try again.")


def _call_llm(cfg: dict[str, str], system: str, user: str) -> str:
    if cfg["provider"] == "gemini":
        return _call_gemini(cfg, system, user)
    if cfg["provider"] == "anthropic":
        return _call_anthropic(cfg, system, user)
    return _call_openai(cfg, system, user)


# ----------------------------------------------------------------- public ---
def discussion(cutoff: float = 95.0, force_refresh: bool = False) -> dict[str, Any]:
    cfg = _llm_config()
    if not cfg["api_key"]:
        return {
            "available": False,
            "note": ("GEMINI_API_KEY not set. Add it to a .env file in the project root "
                     "(GEMINI_API_KEY=your-key) and restart the platform. "
                     "The text is generated only from your enrichment results and the retrieved PubMed abstracts."),
        }

    # Gemini model names change over releases; resolve to one the key can use.
    if cfg["provider"] == "gemini":
        try:
            cfg["model"] = _gemini_pick_model(cfg["api_key"], cfg["model"])
        except Exception:
            pass

    # The cache must also key on the expert decisions: the curated set (and so
    # the whole Discussion) changes when the expert accepts/rejects interactions,
    # even at the same cutoff. Without this, a stale text is returned after curation.
    try:
        decided = assessments.assessed_keys()
        fp = hashlib.md5(
            "|".join(f"{k}:{decided[k]}" for k in sorted(decided)).encode("utf-8")
        ).hexdigest()[:12]
    except Exception:
        fp = "na"
    key = f"{cfg['provider']}|{cfg['model']}|{float(cutoff)}|{fp}"
    if not force_refresh and key in _CACHE:
        return _CACHE[key]

    f = _findings(cutoff)
    if not f["pairs"]:
        return {"available": True, "note": "No curated interactions at this cutoff to discuss."}

    refs = _references(f["pairs"], _context_terms(f.get("context") or {}))
    try:
        text = _call_llm(cfg, _SYSTEM, _user_prompt(f, refs))
    except HTTPError as e:
        detail = e.read().decode("utf-8", "replace")[:300]
        return {"available": True, "error": f"The LLM request failed (HTTP {e.code}). {detail}"}
    except URLError as e:
        return {"available": True, "error": f"Could not reach the LLM API ({e}). Check this machine's internet connection and the API key."}
    except Exception as e:  # noqa: BLE001
        return {"available": True, "error": f"Discussion generation failed: {e}"}

    # Keep only the references the text actually cites, then renumber both the
    # in-text citations and the reference list to a clean contiguous 1..k so the
    # numbers always line up. Fall back to all retrieved if the model cited none.
    cited = {int(n) for grp in re.findall(r"\[(\d+(?:\s*,\s*\d+)*)\]", text)
             for n in re.findall(r"\d+", grp)}
    kept = [r for r in refs if r["n"] in cited]
    if kept:
        remap = {r["n"]: i + 1 for i, r in enumerate(kept)}  # old number -> new (1..k)

        def _renum(m):
            nums = [str(remap[int(x)]) for x in re.findall(r"\d+", m.group(1)) if int(x) in remap]
            return "[" + ", ".join(nums) + "]" if nums else ""

        text = re.sub(r"\[(\d+(?:\s*,\s*\d+)*)\]", _renum, text)
        used = [{**r, "n": remap[r["n"]]} for r in kept]
    else:
        used = refs

    out = {
        "available": True,
        "cutoff": float(cutoff),
        "provider": cfg["provider"],
        "model": cfg["model"],
        "text": text,
        "references": [{k: r[k] for k in ("n", "pmid", "title", "authors", "year", "journal", "url", "pair")} for r in used],
        "n_refs": len(used),
        "generated_at": time.strftime("%Y-%m-%d %H:%M"),
        "note": "",
    }
    _CACHE[key] = out
    return out
