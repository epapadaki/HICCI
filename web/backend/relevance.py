"""
How close is each publication to the user's experiment?

Deterministic and transparent: it compares the organism, tissue, cell types
and condition extracted from PubMed metadata/abstracts with the project
description and with the sender/receiver populations of the selected event.

Each aspect gets one of:
    same / different / partly  - both sides are known
    not_reported               - the publication metadata does not state it
    not_set                    - the user has not described the experiment
Related-but-differently-named cell types (e.g. macrophages vs "Myeloid") are
NOT inferred; the interface says so instead of guessing.
"""

from __future__ import annotations

import re
from typing import Any

_ORGANISM_SYNONYMS = {
    "Human": ("human", "humans", "homo sapiens", "patient", "patients", "hsa"),
    "Mouse": ("mouse", "mice", "murine", "mus musculus", "mmu"),
    "Rat": ("rat", "rats", "rattus norvegicus"),
    "Zebrafish": ("zebrafish", "danio rerio"),
    "Pig": ("pig", "pigs", "swine", "porcine", "sus scrofa"),
    "Primate": ("primate", "primates", "macaque", "macaca"),
}

_GENERIC_CELL_WORDS = {"cell", "cells", "type", "types", "like", "positive", "negative", "high", "other"}

_MODEL_NOTES = {
    "cell line": "Includes cell-line experiments, which may not reflect behaviour in tissue.",
    "xenograft": "Includes xenograft models.",
    "organoid": "Includes organoid models.",
    "in vitro": "Includes in vitro experiments.",
}


def normalize_organism(value: str | None) -> str | None:
    text = (value or "").strip().lower()
    if not text:
        return None
    for canonical, words in _ORGANISM_SYNONYMS.items():
        if text in words or text == canonical.lower():
            return canonical
    return value.strip().capitalize()


def _tokens(text: str) -> set[str]:
    out = set()
    for part in re.split(r"[^A-Za-z0-9]+", (text or "").lower()):
        if len(part) < 4 or part in _GENERIC_CELL_WORDS:
            continue
        out.add(part[:-1] if part.endswith("s") and len(part) > 4 else part)
    return out


def _tissue_terms(tissue: str) -> list[str]:
    terms = []
    full = tissue.strip().lower()
    if full:
        terms.append(full)
    for tok in re.split(r"[^a-z0-9]+", full):
        if len(tok) >= 4:
            terms.append(tok)
            if len(tok) >= 8:
                terms.append(tok[:7])   # endometrium -> endomet (endometrial, …)
    return list(dict.fromkeys(terms))


def _article_text(article: dict[str, Any]) -> str:
    parts = [str(article.get("title") or ""), str(article.get("abstract") or "")]
    parts += [str(x) for x in article.get("tissue_context") or []]
    return " ".join(parts).lower()


def compare(
    article: dict[str, Any],
    ctx: dict[str, Any],
    sender: str,
    receiver: str,
) -> dict[str, Any]:
    items: list[dict[str, Any]] = []
    text = _article_text(article)

    # --- Organism -------------------------------------------------------
    yours = normalize_organism(ctx.get("organism"))
    study_orgs = [normalize_organism(o) for o in article.get("organisms") or []]
    study_orgs = [o for o in dict.fromkeys(study_orgs) if o]
    if not yours:
        status = "not_set"
    elif not study_orgs:
        status = "not_reported"
    elif study_orgs == [yours]:
        status = "same"
    elif yours in study_orgs:
        status = "partly"
    else:
        status = "different"
    items.append({
        "aspect": "Organism", "status": status,
        "study": ", ".join(study_orgs) or None, "yours": yours,
        "note": "The study also uses other organisms." if status == "partly" else None,
    })

    # --- Tissue ---------------------------------------------------------
    tissue = str(ctx.get("tissue_or_organ") or "").strip()
    study_tissue = [str(t) for t in article.get("tissue_context") or []]
    if not tissue:
        status = "not_set"
    elif any(term in text for term in _tissue_terms(tissue)):
        status = "same"
    elif study_tissue:
        status = "different"
    else:
        status = "not_reported"
    items.append({
        "aspect": "Tissue", "status": status,
        "study": ", ".join(study_tissue) or None, "yours": tissue or None, "note": None,
    })

    # --- Cell types (event populations vs studied cell types) -------------
    study_cells = [str(c) for c in article.get("cell_types") or []]
    matches = []
    for population in (sender, receiver):
        pop_tokens = _tokens(population)
        for cell in study_cells:
            if pop_tokens & _tokens(cell):
                matches.append(f"{population} ↔ {cell}")
    matches = list(dict.fromkeys(matches))
    if not study_cells:
        status = "not_reported"
        note = None
    elif matches:
        status = "same"
        note = None
    else:
        status = "no_name_match"
        note = None  # explained once for all publications in the interface
    items.append({
        "aspect": "Cell types", "status": status,
        "study": ", ".join(study_cells) or None,
        "yours": f"{sender} → {receiver}",
        "matches": matches, "note": note,
    })

    # --- Condition ------------------------------------------------------
    condition_source = " ".join(
        [str(ctx.get("disease_or_condition") or "")]
        + [str(g) for g in ctx.get("experimental_groups") or []]
    ).strip()
    if not condition_source:
        status, hits = "not_set", []
    else:
        hits = sorted(t for t in _tokens(condition_source) if t in text)
        status = "mentioned" if hits else "not_mentioned"
    items.append({
        "aspect": "Condition", "status": status,
        "study": None, "yours": condition_source or None,
        "note": ("Terms found in the abstract: " + ", ".join(hits)) if hits else None,
    })

    # --- Study system (informational) -------------------------------------
    systems = [str(s) for s in article.get("experimental_context") or []]
    lowered = " ".join(systems).lower()
    notes = [msg for key, msg in _MODEL_NOTES.items() if key in lowered]
    items.append({
        "aspect": "Study system", "status": "info" if systems else "not_reported",
        "study": ", ".join(systems) or None, "yours": None,
        "note": " ".join(notes) or None,
    })

    return {"items": items, "summary": _summary(items)}


def _summary(items: list[dict[str, Any]]) -> str:
    by = {i["aspect"]: i for i in items}
    parts = []

    org = by["Organism"]
    if org["status"] == "same":
        parts.append(f"the same organism ({org['study']})")
    elif org["status"] == "partly":
        parts.append(f"the same organism among others ({org['study']})")
    elif org["status"] == "different":
        parts.append(f"a different organism ({org['study']})")

    tis = by["Tissue"]
    if tis["status"] == "same":
        parts.append("the same tissue")
    elif tis["status"] == "different":
        parts.append(f"a different tissue or condition ({tis['study']})")

    cells = by["Cell types"]
    if cells["status"] == "same":
        parts.append("cell types that share names with your populations")
    elif cells["status"] == "no_name_match":
        parts.append("cell types with different names from your populations")

    if not parts:
        return "Not enough information in the publication metadata to compare it with your experiment."

    sentence = "The study concerns " + ", ".join(parts[:-1])
    sentence = (sentence + " and " + parts[-1]) if len(parts) > 1 else "The study concerns " + parts[0]
    missing = [i["aspect"].lower() for i in items if i["status"] == "not_reported" and i["aspect"] != "Study system"]
    if missing:
        sentence += ". Not stated in the metadata: " + ", ".join(missing)
    return sentence + "."
