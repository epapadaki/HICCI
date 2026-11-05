"""
Functional enrichment, before vs after curation — fully local and offline.

For every candidate gene (ligand or receptor of any candidate interaction) and
for the curated gene set (kept after expert review), run an over-representation
analysis against local gene-set libraries (GO Biological Process, KEGG,
Reactome, WikiPathways) and return the top terms per source. The result is a
before/after pair per source, mirroring the gene-curator enrichment view.

The statistics (hypergeometric test + Benjamini-Hochberg FDR) run locally with
no network access. The gene-set libraries are plain GMT files under
``data/genesets/``; download them once with ``scripts/fetch_genesets.py`` (or
drop the four Enrichr .gmt files there by hand). Nothing is fetched at query
time, so the view works on a machine with no internet.
"""

from __future__ import annotations

from math import exp, lgamma, log10
from pathlib import Path
from typing import Any

from . import curated, review_queue

# webccc/data/genesets  (functional_enrichment.py is web/backend/...)
DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "genesets"

# (code, label, gmt filename written by scripts/fetch_genesets.py)
# The same 15 libraries the GeneCurator tool reports, grouped GO / pathways /
# disease / drug / cell type.
SOURCES = [
    ("GO:BP", "GO Biological Process", "GO_Biological_Process.gmt"),
    ("GO:CC", "GO Cellular Component", "GO_Cellular_Component.gmt"),
    ("GO:MF", "GO Molecular Function", "GO_Molecular_Function.gmt"),
    ("KEGG", "KEGG", "KEGG.gmt"),
    ("REAC", "Reactome", "Reactome.gmt"),
    ("WP", "WikiPathways", "WikiPathways.gmt"),
    ("BIOPLANET", "BioPlanet", "BioPlanet.gmt"),
    ("HALLMARK", "MSigDB Hallmark", "MSigDB_Hallmark.gmt"),
    ("DISGENET", "DisGeNET (disease)", "DisGeNET.gmt"),
    ("DISEASE_UP", "Disease Perturbations (up)", "Disease_Perturbations_up.gmt"),
    ("DISEASE_DN", "Disease Perturbations (down)", "Disease_Perturbations_down.gmt"),
    ("DSIGDB", "DSigDB (drugs)", "DSigDB.gmt"),
    ("DRUG_UP", "Drug Perturbations (up)", "Drug_Perturbations_up.gmt"),
    ("DRUG_DN", "Drug Perturbations (down)", "Drug_Perturbations_down.gmt"),
]

# term size limits for a sensible over-representation analysis
_MIN_TERM = 3
_MAX_TERM = 2000
_MIN_OVERLAP = 2
_FDR = 0.05
_TOP = 10

_LIB_CACHE: dict[str, dict[str, frozenset] | None] = {}


# ---------------------------------------------------------------- gene sets ---
def _split(gene: str) -> list[str]:
    out = []
    for part in str(gene).replace(" complex", "").split("+"):
        g = part.strip()
        if g:
            out.append(g)
    return out


def _genes_from(df) -> list[str]:
    if df is None or getattr(df, "empty", True):
        return []
    s: set[str] = set()
    for col in ("ligand", "receptor"):
        if col not in df.columns:
            continue
        for v in df[col].astype(str):
            s.update(_split(v))
    return sorted({g.upper() for g in s if g})


def _load_gmt(path: Path) -> dict[str, frozenset]:
    sets: dict[str, frozenset] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            genes = {p.split(",")[0].strip().upper() for p in parts[2:]}
            genes = {g for g in genes if g}
            if name and _MIN_TERM <= len(genes) <= _MAX_TERM:
                sets[name] = frozenset(genes)
    return sets


def _library(fname: str) -> dict[str, frozenset] | None:
    if fname not in _LIB_CACHE:
        path = DATA_DIR / fname
        _LIB_CACHE[fname] = _load_gmt(path) if path.exists() else None
    return _LIB_CACHE[fname]


# --------------------------------------------------------------- statistics ---
def _logC(a: int, b: int) -> float:
    if b < 0 or b > a:
        return float("-inf")
    return lgamma(a + 1) - lgamma(b + 1) - lgamma(a - b + 1)


def _hyper_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(N population, K successes, n draws)."""
    lo = max(0, n - (N - K))
    hi = min(K, n)
    if k <= lo:
        return 1.0
    if k > hi:
        return 0.0
    denom = _logC(N, n)
    terms = [_logC(K, i) + _logC(N - K, n - i) - denom for i in range(k, hi + 1)]
    mx = max(terms)
    s = sum(exp(t - mx) for t in terms)
    return min(1.0, exp(mx) * s)


def _finalize(rows: list[tuple[str, float, int]]) -> list[dict[str, Any]]:
    """Benjamini-Hochberg across the tested terms, keep FDR < threshold, top N."""
    if not rows:
        return []
    rows.sort(key=lambda r: r[1])
    m = len(rows)
    fdr = [0.0] * m
    prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, rows[i][1] * m / (i + 1))
        fdr[i] = min(1.0, prev)
    out: list[dict[str, Any]] = []
    for (name, _p, k), q_val in zip(rows, fdr):
        if q_val > _FDR:
            continue
        out.append({"name": name, "neglog10p": round(-log10(max(q_val, 1e-300)), 2), "n": int(k)})
        if len(out) >= _TOP:
            break
    return out


def _enrich_source(lib: dict[str, frozenset], query: list[str],
                   background: list[str]) -> list[dict[str, Any]]:
    """Inferential ORA: `query` within the candidate-L-R `background` (both
    restricted to genes annotated in this library). This is the statistically
    appropriate test — it controls for ligands/receptors being intrinsically
    signalling-heavy, so it reflects what the curated subset is *specifically*
    enriched for."""
    universe: set[str] = set()
    for gs in lib.values():
        universe |= gs
    bg = set(background) & universe
    q = set(query) & bg
    N, n = len(bg), len(q)
    if n < _MIN_TERM or N < _MIN_TERM:
        return []
    rows: list[tuple[str, float, int]] = []
    for name, gs in lib.items():
        K = len(gs & bg)
        if K < 1:
            continue
        k = len(gs & q)
        if k < _MIN_OVERLAP:
            continue
        rows.append((name, _hyper_sf(k, N, K, n), k))
    return _finalize(rows)


def _enrich_vs_genome(lib: dict[str, frozenset], query: list[str]) -> list[dict[str, Any]]:
    """Descriptive ORA against the whole library universe (genome-wide background).
    For an L-R gene set this is biased toward generic signalling terms, so it is
    used ONLY for a clearly-labelled 'compositional shift' view (before vs after),
    never as the inferential result."""
    universe: set[str] = set()
    for gs in lib.values():
        universe |= gs
    q = set(query) & universe
    N, n = len(universe), len(q)
    if n < _MIN_TERM or N == 0:
        return []
    rows: list[tuple[str, float, int]] = []
    for name, gs in lib.items():
        k = len(gs & q)
        if k < _MIN_OVERLAP:
            continue
        rows.append((name, _hyper_sf(k, N, len(gs), n), k))
    return _finalize(rows)


# ------------------------------------------------------------------- public ---
def functional_enrichment(cutoff: float = 95.0, descriptive: bool = False) -> dict[str, Any]:
    # Inferential result: background = all candidate ligand-receptor genes, query =
    # curated genes. Optional descriptive view (genome background) adds the
    # before/after compositional shift when `descriptive` is requested.
    background = _genes_from(review_queue._events())
    after = _genes_from(curated._kept(cutoff))

    sources: list[dict[str, Any]] = []
    missing: list[str] = []
    any_lib = False
    for code, label, fname in SOURCES:
        lib = _library(fname)
        if lib is None:
            missing.append(label)
            sources.append({"source": code, "label": label, "installed": False, "after": []})
            continue
        any_lib = True
        src = {
            "source": code,
            "label": label,
            "installed": True,
            "after": _enrich_source(lib, after, background),
        }
        if descriptive:
            src["before_genome"] = _enrich_vs_genome(lib, background)
            src["after_genome"] = _enrich_vs_genome(lib, after)
        sources.append(src)

    if not any_lib:
        return {
            "cutoff": float(cutoff),
            "before_n": len(background),
            "after_n": len(after),
            "reached": False,
            "sources": sources,
            "note": (
                "Gene-set libraries are not installed yet. Run  "
                "python scripts/fetch_genesets.py  once (it downloads the 15 libraries "
                "— GO, pathways, disease, drug and cell-type — into data/genesets/), then "
                "reload this page. No internet on this machine? Download the libraries "
                "from Enrichr in a browser and drop the .gmt files into data/genesets/ — "
                "nothing else is fetched at run time."
            ),
        }

    note = (
        "Over-representation of the curated, expert-reviewed genes against the background "
        "of all candidate ligand–receptor genes (not the whole genome), so terms that are "
        "common to ligands/receptors in general are controlled for and the result reflects "
        "what the curated subset is specifically enriched for. Hypergeometric test, "
        "Benjamini–Hochberg FDR < 0.05 per library; runs fully offline."
    )
    if missing:
        note += " Not installed: " + ", ".join(missing) + " (run scripts/fetch_genesets.py to add them)."

    return {
        "cutoff": float(cutoff),
        "before_n": len(background),
        "after_n": len(after),
        "descriptive": bool(descriptive),
        "reached": True,
        "sources": sources,
        "note": note,
    }
