"""
Biological documentation for one candidate interaction.

This module only assembles information. The scientific lookups are the
existing, tested functions in ``utils/`` (local UniProt, Reactome and SIGNOR
snapshots). Every item carries a ``scope`` so the interface can show whether
it concerns the specific ligand-receptor pair or receptor signalling in
general.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

from . import config  # noqa: F401  (puts the project root on sys.path)
from . import index_store, score_definitions

from utils.biological_context import get_biological_context
from utils.reactome_local import get_receptor_pathways
from utils.signor_downstream_local import get_downstream_causal_neighborhood

SCOPE_PAIR = "pair"                  # evidence about this ligand-receptor pair
SCOPE_RECEPTOR = "receptor_general"  # evidence about the receptor in general

MAX_SIGNOR_DIRECT = 40


# ---------------------------------------------------------------------------
# Header / computational result
# ---------------------------------------------------------------------------
def overview(result_id: int) -> dict[str, Any] | None:
    rec = index_store.get_result(result_id)
    if not rec:
        return None

    method = rec["method"]
    signaling = None
    if rec.get("signaling_pathway") or rec.get("signaling_annotation"):
        signaling = {
            "pathway": rec.get("signaling_pathway"),
            "annotation": rec.get("signaling_annotation"),
            "evidence": rec.get("signaling_evidence"),
            "note": (
                f"Category assigned by {method} from its ligand-receptor database. "
                "It groups related ligand-receptor pairs. It is not a list of "
                "intracellular pathways activated in the receiver, and it does not "
                "show that this pathway is active in your data."
            ),
        }

    return {
        "id": rec["id"],
        "event_key": rec["event_key"],
        "sender": rec["sender"],
        "receiver": rec["receiver"],
        "ligand": index_store.entity_display(rec["ligand_entity_key"], rec.get("m_ligand")),
        "receptor": index_store.entity_display(rec["receptor_entity_key"], rec.get("m_receptor")),
        "method": method,
        "method_labels": {"ligand": rec.get("m_ligand"), "receptor": rec.get("m_receptor")},
        "score": rec.get("score"),
        "score_info": score_definitions.describe_score(method, rec.get("score_name")),
        "global_percentile": rec.get("global_percentile"),
        "cell_pair_percentile": rec.get("cell_pair_percentile"),
        "percentile_definition": score_definitions.PERCENTILE_DEFINITION,
        "method_count": rec.get("method_count"),
        "other_methods": [
            {**o, "score_info": score_definitions.describe_score(o["method"], o.get("score_name"))}
            for o in rec.get("other_methods", [])
        ],
        "signaling": signaling,
        "computational_note": (
            "This is a computational candidate: the method found expression of the "
            "ligand in the sender and of the receptor in the receiver consistent with "
            "communication through this pair. A high score does not show that the "
            "proteins bind, that the receptor is activated, or that any downstream "
            "pathway is active. scRNA-seq measures RNA; binding happens between the "
            "encoded proteins."
        ),
    }


# ---------------------------------------------------------------------------
# Section A: biological meaning and pathways
# ---------------------------------------------------------------------------
def _molecule(bc: dict[str, Any], side: str) -> dict[str, Any]:
    subunits = []
    for ann in bc.get(f"{side}_subunit_annotations") or []:
        if not ann or str(ann.get("found")) in ("False", "false"):
            continue
        subunits.append({
            "gene": ann.get("gene_symbol"),
            "uniprot": ann.get("entry"),
            "protein_name": ann.get("protein_name"),
            "function": ann.get("function"),
            "location": ann.get("subcellular_location") if ann.get("subcellular_location") not in (None, "None") else None,
        })
    return {
        "label": bc.get(side),
        "entity_key": bc.get(f"{side}_entity_key"),
        "entity_type": bc.get(f"{side}_entity_type"),
        "required_subunits": bc.get(f"{side}_required_subunits") or [],
        "summary": bc.get(f"{side}_description"),
        "subunits": subunits,
    }


def _direction(effect: str) -> str:
    e = (effect or "").lower()
    if e.startswith("up-regulates"):
        return "activation"
    if e.startswith("down-regulates"):
        return "inhibition"
    if e:
        return "other"
    return "unknown"


@lru_cache(maxsize=512)
def _reactome_for(gene: str) -> dict[str, Any]:
    try:
        res = get_receptor_pathways(gene)
    except Exception as exc:  # unresolved symbol, missing snapshot, …
        return {"gene": gene, "pathways": [], "error": str(exc)}
    pathways, seen = [], set()
    for p in res.get("pathways") or []:
        rid = p.get("reactome_id")
        if rid in seen:          # the snapshot can list a pathway more than once
            continue
        seen.add(rid)
        pathways.append({
            "id": rid,
            "name": p.get("pathway_name"),
            "description": p.get("description"),
            "url": f"https://reactome.org/content/detail/{rid}" if rid else None,
            "link_to_receptor": f"{gene} is annotated as a participant in this pathway.",
            "direction": "not_specified",
            "context": "Curated human pathway (Reactome)",
            "scope": SCOPE_RECEPTOR,
        })
    return {"gene": gene, "pathways": pathways, "error": None}


@lru_cache(maxsize=512)
def _signor_for(gene: str) -> dict[str, Any]:
    try:
        res = get_downstream_causal_neighborhood(gene, max_depth=2)
    except Exception as exc:
        return {"gene": gene, "direct": [], "indirect_count": 0, "error": str(exc)}

    direct, seen, indirect = [], set(), 0
    for rel in res.get("relations") or []:
        if rel.get("depth") != 1:
            indirect += 1
            continue
        key = (rel.get("target_entity"), rel.get("effect"))
        if key in seen:
            continue
        seen.add(key)
        pmid = str(rel.get("pmid") or "").strip()
        direct.append({
            "target": rel.get("target_entity"),
            "target_type": rel.get("target_type"),
            "effect": rel.get("effect"),
            "direction": _direction(rel.get("effect")),
            "mechanism": rel.get("mechanism") or None,
            "cell_context": rel.get("cell_context") or None,
            "tissue_context": rel.get("tissue_context") or None,
            "sentence": rel.get("evidence_sentence") or None,
            "signor_id": rel.get("signor_id"),
            "pmid": pmid or None,
            "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid.isdigit() else None,
            "scope": SCOPE_RECEPTOR,
        })
    # Activation/inhibition first, then others; keep the list readable.
    order = {"activation": 0, "inhibition": 1, "other": 2, "unknown": 3}
    direct.sort(key=lambda r: (order[r["direction"]], str(r["target"])))
    return {
        "gene": gene,
        "direct": direct[:MAX_SIGNOR_DIRECT],
        "direct_total": len(direct),
        "indirect_count": indirect,
        "error": None,
    }


def biology(result_id: int) -> dict[str, Any] | None:
    rec = index_store.get_result(result_id)
    if not rec:
        return None

    bc = get_biological_context(
        ligand=rec["ligand"],
        receptor=rec["receptor"],
        sender_cell_type=rec["sender"],
        receiver_cell_type=rec["receiver"],
        ligand_entity_key=rec["ligand_entity_key"],
        receptor_entity_key=rec["receptor_entity_key"],
    )

    ligand = _molecule(bc, "ligand")
    receptor = _molecule(bc, "receptor")
    receptor_genes = receptor["required_subunits"]

    return {
        "ligand": ligand,
        "receptor": receptor,
        "pair": {
            "text": bc.get("possible_role"),
            "scope": SCOPE_PAIR,
            "note": (
                "Built from the curated annotations of each molecule. Evidence that "
                "concerns this exact pair comes from the literature (next tab)."
            ),
        },
        "receptor_pathways": {
            "genes": receptor_genes,
            "reactome": [_reactome_for(g) for g in receptor_genes],
            "signor": [_signor_for(g) for g in receptor_genes],
            "scope": SCOPE_RECEPTOR,
            "note": (
                "These pathways are linked to the receptor"
                + (" subunits" if len(receptor_genes) > 1 else "")
                + " in general, whatever ligand binds it. They are not attributed to "
                f"{ligand['label']} specifically, and they are not shown to be active "
                "in your data."
            ),
            "unresolved_note": (
                None if receptor_genes else
                "The receptor could not be resolved to verified genes, so no "
                "receptor pathways are shown rather than guessing a single gene."
            ),
        },
        "sources": {
            "uniprot": bc.get("source"),
            "reactome": "Reactome (local snapshot)",
            "signor": "SIGNOR (local snapshot)",
        },
        "caution": bc.get("caution"),
    }


def warm_up() -> None:
    """Load the local UniProt snapshot once at startup so the first click is fast."""
    try:
        get_biological_context(ligand="APP", receptor="CD74",
                               ligand_entity_key="GENE:APP", receptor_entity_key="GENE:CD74")
        _reactome_for("CD74")
        _signor_for("CD74")
    except Exception:
        pass
