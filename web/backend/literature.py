"""
Where has this ligand-receptor connection been studied?

Uses the existing entity-aware PubMed search (complexes are searched by their
subunits, never collapsed to one gene) and the existing local, rule-based
evidence classifier (no LLM, no API key). Results are cached by the
underlying modules, so repeated views are instant and work offline.

Publications where the two molecules only co-occur are kept but flagged, and
are never presented as evidence that the interaction exists.
"""

from __future__ import annotations

import re
from typing import Any

from . import config
from . import index_store, project_context, relevance

from utils.pubmed_literature import search_pubmed_pair_entity_aware
from utils.publication_evidence import classify_publication_evidence
from utils import gemini_client

EVIDENCE_LABELS = {
    "direct_binding": "Direct binding",
    "perturbation_evidence": "Perturbation (block, knockout, knockdown)",
    "functional_response": "Functional response",
    "expression_or_coexpression": "Expression / co-expression",
    "computational_prediction": "Computational prediction",
    "pair_co_mention_only": "Co-mention only",
    "insufficient_information": "Not enough information",
}

STUDY_MODE_LABELS = {
    "experimental": "Experimental",
    "computational": "Computational",
    "mixed": "Experimental and computational",
    "unclear": "Unclear",
}

DIRECTION_LABELS = {
    "activation": "activation",
    "inhibition": "inhibition",
    "modulation": "modulation",
    "association": "association only",
    "unclear": "unclear",
}


def _whole_word(term: str, text: str) -> bool:
    """Guard against keyword hits inside other words (e.g. 'STAT' in 'state')."""
    if not term:
        return False
    flags = 0 if term.isupper() else re.IGNORECASE
    return re.search(r"(?<![A-Za-z0-9])" + re.escape(term) + r"(?![A-Za-z0-9])", text, flags) is not None


def literature(result_id: int, *, force_refresh: bool = False) -> dict[str, Any] | None:
    rec = index_store.get_result(result_id)
    if not rec:
        return None

    lig = index_store.entity_display(rec["ligand_entity_key"], rec.get("m_ligand"))
    rcp = index_store.entity_display(rec["receptor_entity_key"], rec.get("m_receptor"))

    if lig["kind"] in ("proxy", "other", "unknown") or rcp["kind"] in ("proxy", "other", "unknown"):
        return {
            "status": "unsupported",
            "message": (
                "Literature search needs both partners resolved to genes or verified "
                "complex subunits. "
                + ("The ligand is a small molecule estimated from enzyme expression. " if lig["kind"] == "proxy" else "")
                + ("The receptor is a CellPhoneDB entity without verified subunits. " if rcp["kind"] == "other" else "")
                + "Searching by its label could return misleading results, so no search was run."
            ),
            "publications": [],
            "pair_pathways": [],
        }

    try:
        search = search_pubmed_pair_entity_aware(
            lig["label"],
            rcp["label"],
            ligand_entity_key=lig["key"],
            receptor_entity_key=rcp["key"],
            ligand_required_subunits=lig["subunits"],
            receptor_required_subunits=rcp["subunits"],
            max_results=config.PUBMED_MAX_RESULTS,
            force_refresh=force_refresh,
        )
    except Exception as exc:
        return {
            "status": "error",
            "message": (
                "PubMed could not be searched. Check the internet connection and try "
                f"again. Details: {exc}"
            ),
            "publications": [],
            "pair_pathways": [],
        }

    articles = search.get("articles") or []
    classes: dict[str, dict[str, Any]] = {}
    classified_by = None
    if articles:
        result = None
        if gemini_client.available():
            # AI-assisted classification (Gemini); any failure falls back to the rules.
            try:
                result = classify_publication_evidence(lig["label"], rcp["label"], articles, use_llm=True)
            except Exception:
                result = None
        if result is None:
            try:
                result = classify_publication_evidence(lig["label"], rcp["label"], articles)
            except Exception:
                result = None
        if result:
            classes = {str(c.get("pmid")): c for c in result.get("classifications") or []}
            classified_by = {"ai": not result.get("local_rule_based", False),
                             "provider": result.get("provider") or "Local rule-based",
                             "model": result.get("model")}

    ctx = project_context.load()
    publications, pair_pathways = [], []

    for art in articles:
        pmid = str(art.get("pmid") or "")
        cls = classes.get(pmid, {})
        types = cls.get("evidence_types") or []
        primary = cls.get("primary_evidence_type")
        pair_status = cls.get("pair_specific_evidence") or "unclear"
        co_mention = primary == "pair_co_mention_only" or cls.get("evidence_scope") == "co_mention_only"

        publications.append({
            "pmid": pmid or None,
            "doi": art.get("doi") or None,
            "title": art.get("title"),
            "year": art.get("year"),
            "journal": art.get("journal"),
            "url": art.get("pubmed_url") or (f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None),
            "doi_url": f"https://doi.org/{art['doi']}" if art.get("doi") else None,
            "organisms": art.get("organisms") or [],
            "tissue": art.get("tissue_context") or [],
            "cell_types": art.get("cell_types") or [],
            "study_system": art.get("experimental_context") or [],
            "what_studied": cls.get("what_was_studied") or art.get("pair_evidence_snippet"),
            "evidence": {
                "primary": EVIDENCE_LABELS.get(primary, "Not classified"),
                "types": [EVIDENCE_LABELS.get(t, t) for t in types],
                "pair_specific": pair_status,
                "direct_binding": cls.get("direct_binding_evidence"),
                "study_mode": STUDY_MODE_LABELS.get(cls.get("study_mode"), None),
                "method": cls.get("classification_reason"),
            },
            "pair_specific": pair_status == "supported" and not co_mention,
            "co_mention_only": co_mention,
            "abstract_only": True,
            "relevance": relevance.compare(art, ctx, rec["sender"], rec["receiver"]),
        })

        source_text = " ".join([str(art.get("title") or ""), str(art.get("abstract") or "")])
        for pw in cls.get("pair_specific_pathways") or []:
            name = str(pw.get("pathway_name") or "").strip()
            if not name or not _whole_word(name, source_text):
                continue
            pair_pathways.append({
                "name": name,
                "direction": DIRECTION_LABELS.get(pw.get("effect_direction"), "unclear"),
                "context": pw.get("biological_context"),
                "pmid": pmid,
                "title": art.get("title"),
                "year": art.get("year"),
                "sentence": pw.get("source_sentence"),
                "scope": "pair",
            })

    # Pair-specific studies first, co-mentions last, newest first within groups.
    def _rank(p: dict[str, Any]) -> tuple:
        group = 0 if p["pair_specific"] else (2 if p["co_mention_only"] else 1)
        return (group, -int(p["year"]) if str(p.get("year") or "").isdigit() else 0)

    publications.sort(key=_rank)

    qmeta = search.get("query_metadata") or {}
    return {
        "status": "ok",
        "query": search.get("query"),
        "query_scope": qmeta.get("scope") if isinstance(qmeta, dict) else None,
        "query_mode": qmeta.get("query_mode") if isinstance(qmeta, dict) else None,
        "cached": bool(search.get("cached")),
        "fetched_at": search.get("fetched_at"),
        "publications": publications,
        "pair_pathways": pair_pathways,
        "project_context_filled": project_context.is_filled(ctx),
        "classified_by": classified_by,
        "notes": [
            "Based on PubMed titles, abstracts and metadata only; full texts were not read.",
            ("Evidence types are assigned by an AI model (" + str(classified_by.get("model") or "Gemini")
             + ") from the title and abstract, with output checked against a fixed schema. "
             "Check the publication before relying on them.")
            if (classified_by or {}).get("ai") else
            "Evidence types are assigned by local keyword rules. Check the publication before relying on them.",
            "A publication that mentions both molecules is not, by itself, evidence that they interact.",
        ],
    }
