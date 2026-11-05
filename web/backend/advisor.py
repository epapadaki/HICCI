"""
Which CCI tools suit a dataset?

Transparent rules based on the dataset's own properties (checked when it was
prepared) and on documented properties of each tool. There is no single best
tool: benchmarks show that tools and their ligand-receptor databases disagree
substantially, so the advice is about suitability and trade-offs, and it
suggests combining tools with different databases.

AI assistance (optional): when a Gemini key is configured (GEMINI_API_KEY in
.env), a language model reviews the same dataset facts, the documented tool
profiles and the rule-based assessment, and proposes one method with a short
rationale. Its choice is validated: it must be one of the five tools and not a
tool the rules mark as unsuitable for this data; otherwise the rule-based
recommendation is kept. The rule-based recommendation is always returned too,
and the user makes the final choice.
"""

from __future__ import annotations

import json
from typing import Any

from . import datasets
from utils import gemini_client

_AI_CACHE: dict[str, dict[str, Any]] = {}

_AI_SYSTEM = (
    "You are a bioinformatics expert advising which cell-cell communication inference method "
    "suits a single-cell RNA-seq dataset. You are given the dataset facts, the documented profile "
    "of each available method, a rule-based suitability assessment and the rule-based recommendation.\n"
    "RULES:\n"
    "- Recommend exactly ONE method, chosen only from the methods listed.\n"
    "- Never recommend a method whose suitability level is 'poor'.\n"
    "- Base the choice on the dataset facts and the documented method properties only; do not claim "
    "that any method is more accurate in general, since benchmarks show substantial disagreement between methods.\n"
    "- Give 2-4 short sentences of rationale that refer to concrete facts of this dataset "
    "(e.g. population sizes, immune share, dataset size, need for complexes or pathway annotation).\n"
    "- Optionally name one complementary method with a different database and one sentence why, or null.\n"
    "- Output valid JSON only: {\"tool\": str, \"why\": [str, ...], "
    "\"complementary\": {\"tool\": str, \"why\": str} or null}"
)


def _ai_recommend(ds, facts: dict[str, Any], per_type: dict[str, Any], items: list[dict[str, Any]],
                  rule_rec: str, rule_why: list[str]) -> dict[str, Any] | None:
    """AI-assisted recommendation, validated against the rules. None if unavailable or invalid."""
    if not gemini_client.available():
        return None
    levels = {i["tool"]: i["level"] for i in items}
    payload = {
        "dataset": {"label": ds.label, "organism": ds.meta.get("organism"), "tissue": ds.meta.get("tissue"),
                    **facts, "cells_per_population": per_type},
        "methods": [{"tool": i["tool"], "suitability": i["level"], "database": i["database"],
                     "strength": i["strength"], "caveats": i["reasons"]} for i in items],
        "rule_based_recommendation": {"tool": rule_rec, "why": rule_why},
    }
    key = json.dumps(payload, sort_keys=True, default=str) + "|" + gemini_client.model_name()
    if key in _AI_CACHE:
        return _AI_CACHE[key]
    try:
        data, model = gemini_client.generate_json(_AI_SYSTEM, json.dumps(payload, default=str), temperature=0.2)
    except Exception:
        return None
    tool = str(data.get("tool") or "").strip()
    match = next((t for t in TOOLS if t.lower() == tool.lower()), None)
    if not match or levels.get(match) == "poor":
        return None                                  # invalid choice: keep the rule-based one
    why = [str(w)[:400] for w in (data.get("why") or []) if str(w).strip()][:4]
    if not why:
        return None
    comp = data.get("complementary") or None
    if isinstance(comp, dict):
        ct = next((t for t in TOOLS if t.lower() == str(comp.get("tool") or "").strip().lower()), None)
        comp = ({"tool": ct, "why": str(comp.get("why") or "")[:300]}
                if ct and ct != match and levels.get(ct) != "poor" else None)
    else:
        comp = None
    out = {"tool": match, "why": why, "complementary": comp, "source": "ai", "model": model}
    _AI_CACHE[key] = out
    return out

TOOLS = ["CellPhoneDB", "NATMI", "CellChat", "SingleCellSignalR", "ICELLNET"]

PROFILE = {
    "CellPhoneDB": {"complexes": True, "database": "CellPhoneDB v5 (curated, with complexes)",
                    "strength": "Curated database with multi-subunit complexes; fast."},
    "NATMI": {"complexes": False, "database": "connectomeDB2020 (lrc2p)",
              "strength": "Expression-weight and specificity scores; fast; gene-level pairs."},
    "CellChat": {"complexes": True, "database": "CellChatDB v2 (with co-factors and pathway annotation)",
                 "strength": "Permutation p-values and signalling-pathway annotation; slower."},
    "SingleCellSignalR": {"complexes": False, "database": "LRdb (SingleCellSignalR v1)",
                          "strength": "Simple, regularised LRscore; gene-level pairs."},
    "ICELLNET": {"complexes": True, "database": "ICELLNET database (human, experimentally demonstrated, with complexes)",
                 "strength": "Strictly curated interactions grouped into families (cytokines, chemokines, checkpoints…)."},
}


def advise(ds_id: str, installed: dict[str, bool] | None = None) -> dict[str, Any] | None:
    ds = datasets.get(ds_id)
    if not ds:
        return None
    v = ds.meta.get("validation") or {}
    n_cells = int(v.get("n_cells") or 0)
    n_genes = int(v.get("n_genes") or 0)
    per_type = v.get("cells_per_type") or {}
    n_types = len(per_type)
    smallest = min(per_type.values()) if per_type else None
    small = sorted([k for k, n in per_type.items() if n < 10])
    raw = (v.get("counts") or {}).get("looks_like_raw_counts")
    organism = str(ds.meta.get("organism") or "").lower()
    dense_gb = n_cells * n_genes * 8 / 1e9 if n_cells and n_genes else 0

    # How immune is this dataset? (used to recommend ICELLNET, which is curated for
    # immune communication). Based on the cell-type names the user annotated.
    IMMUNE_TERMS = ("t cell", "t-cell", "cd4", "cd8", "treg", "regulatory t",
                    "b cell", "b-cell", "plasma", "nk", "natural killer",
                    "macrophage", "monocyte", "dendritic", "myeloid", "mast",
                    "neutrophil", "granulocyte", "lymphocyte", "kupffer",
                    "microglia", "langerhans", "leukocyte", "immune")
    def _is_immune(name: str) -> bool:
        n = str(name).lower()
        return any(term in n for term in IMMUNE_TERMS)
    immune_types = sorted([k for k in per_type if _is_immune(k)])
    immune_frac = (len(immune_types) / n_types) if n_types else 0.0

    notes: list[str] = []
    if not v:
        notes.append("This dataset has not been checked yet, so the advice relies on general properties only.")
    if organism and organism not in ("human", "homo sapiens", "hsa"):
        notes.append(f"The tools are set up for human data; this dataset is “{ds.meta.get('organism')}”.")
    if raw:
        notes.append("The matrix looks like raw counts. The tools compare expression levels, so a normalised "
                     "matrix (or layer) is preferable; choose it when preparing the dataset.")
    if small:
        notes.append(f"{len(small)} population(s) have fewer than 10 cells ({', '.join(small[:5])}{'…' if len(small) > 5 else ''}); "
                     "their predictions are unreliable with any tool.")

    items = []
    for t in TOOLS:
        reasons, level = [], "good"
        if t == "CellChat":
            if small:
                reasons.append("CellChat drops populations with fewer than 10 cells when filtering communication.")
            if n_types > 30 or n_cells > 50000:
                reasons.append("Many populations or cells: the permutation test will take long.")
                level = "caution"
        if t == "SingleCellSignalR" and dense_gb > 4:
            reasons.append(f"Version 1 holds the whole matrix in memory as a dense table (about {dense_gb:.0f} GB for this dataset).")
            level = "caution" if dense_gb <= 16 else "poor"
        if t == "NATMI" and n_cells > 100000:
            reasons.append("NATMI reads the matrix as a text table; very large datasets are slow.")
            level = "caution"
        if t == "ICELLNET":
            reasons.append("Covers fewer interactions (strict curation, focused on families such as cytokines and checkpoints): "
                           "strong for immune questions, less complete elsewhere.")
        if not PROFILE[t]["complexes"]:
            reasons.append("Works on single genes: multi-subunit receptors (for example integrins, cytokine receptors) are split.")
        if organism and organism not in ("human", "homo sapiens", "hsa"):
            level = "poor"
        items.append({"tool": t, "level": level, "database": PROFILE[t]["database"],
                      "strength": PROFILE[t]["strength"], "reasons": reasons})

    # One recommended tool for this dataset, chosen from the dataset's own
    # properties with transparent rules. There is no single best tool; each rule
    # reflects a documented strength of the tool it points to.
    levels = {i["tool"]: i["level"] for i in items}
    heavy = n_types > 30 or n_cells > 50000
    huge = n_cells > 120000
    tiny = (0 < n_cells <= 2500) and dense_gb <= 1.0
    why: list[str] = []
    if not per_type:
        rec = "CellChat"
        why.append("This dataset has not been checked yet; CellChat is the default choice: it gives a p-value for each "
                   "interaction and annotates it with a signalling pathway.")
    elif immune_frac >= 0.6:
        rec = "ICELLNET"
        why.append(f"{len(immune_types)} of {n_types} populations are immune ({immune_frac*100:.0f}%): "
                   f"{', '.join(immune_types[:6])}{'…' if len(immune_types) > 6 else ''}.")
        why.append("ICELLNET is curated specifically for immune communication (cytokines, chemokines, checkpoints) "
                   "and its interactions are experimentally validated and grouped into families.")
        why.append("It keeps multi-subunit complexes whole.")
    elif small:
        rec = "CellPhoneDB"
        why.append(f"{len(small)} population(s) have fewer than 10 cells; CellPhoneDB keeps them "
                   "(CellChat would drop those).")
        why.append("Its curated database keeps multi-subunit complexes whole, and it is fast.")
    elif huge:
        rec = "NATMI"
        why.append(f"Very large dataset ({n_cells:,} cells); NATMI computes specificity-weighted edge scores "
                   "quickly and scales where CellChat's permutation test would be slow.")
        why.append("Note: it works on single genes, so multi-subunit receptors are split.")
    elif heavy:
        rec = "CellPhoneDB"
        why.append(f"Large dataset ({n_cells:,} cells / {n_types} populations); CellPhoneDB is fast while "
                   "CellChat's permutation test would take long, and it keeps complexes.")
    elif tiny:
        rec = "SingleCellSignalR"
        why.append(f"Small dataset ({n_cells:,} cells); a simple, regularised LRscore is quick and sufficient here, "
                   "and the matrix fits easily in memory.")
        why.append("Note: it works on single genes, so multi-subunit receptors are split.")
    else:
        rec = "CellChat"
        why += ["Every population has enough cells for CellChat's permutation test, which gives a p-value for each interaction.",
                "Its database models complexes and co-factors and annotates each interaction with a signalling pathway.",
                "The dataset size suits it."]
    # Never recommend a tool that is unsuitable for this data (e.g. SingleCellSignalR
    # out of memory, or any tool on non-human data): fall back to the best available.
    if levels.get(rec) == "poor":
        rec = next((t for t in ["CellChat", "CellPhoneDB", "ICELLNET", "NATMI", "SingleCellSignalR"] if levels.get(t) != "poor"), rec)
    WHEN = {
        "CellPhoneDB": "Fast, curated complexes; reports every expressed ligand–receptor pair with a mean-expression score.",
        "CellChat": "Statistical test and signalling-pathway annotation; slower on large datasets.",
        "ICELLNET": "Strictly curated, human-validated interactions; best for immune questions (cytokines, chemokines, checkpoints).",
        "NATMI": "Specificity-weighted scores; fast; works on single genes.",
        "SingleCellSignalR": "Simple regularised score; works on single genes.",
    }
    facts = {"cells": n_cells, "genes": n_genes, "populations": n_types, "smallest_population": smallest,
             "raw_counts": raw, "immune_populations": immune_types}
    rule_based = {"tool": rec, "why": why, "source": "rules"}
    recommended = _ai_recommend(ds, facts, per_type, items, rec, why) or rule_based
    rec = recommended["tool"]
    alternatives = [{"tool": t, "when": WHEN[t], "level": levels.get(t)} for t in TOOLS if t != rec]
    suggested = [rec]
    return {
        "dataset": ds.label,
        "facts": {"cells": n_cells, "genes": n_genes, "populations": n_types, "smallest_population": smallest,
                  "raw_counts": raw, "organism": ds.meta.get("organism")},
        "notes": notes,
        "tools": items,
        "recommended": recommended,
        "rule_based": rule_based,
        "ai_available": gemini_client.available(),
        "alternatives": alternatives,
        "suggested": suggested,
        "caution": ("There is no single best tool: comparisons of methods and resources show substantial "
                    "disagreement between them (Dimitrov et al., Nature Communications 2022). This advice is "
                    "about suitability for your data, not accuracy."),
    }
