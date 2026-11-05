"""
A biological picture of one ligand-receptor interaction.

Everything drawn comes from a named source, and "known mechanism" is kept
apart from "predicted in your sample":

  mechanism type  CellChatDB annotation of the exact pair (Secreted Signaling,
                  Cell-Cell Contact, ECM-Receptor, Non-protein Signaling); else
                  the annotation of other CellChatDB pairs with the same
                  ligand; else the ligand's UniProt subcellular location
  receptor        all subunits of the receptor entity (complexes kept whole),
                  plus CellChatDB co-stimulatory / co-inhibitory receptors
  signalling      direct downstream relations of the receptor subunits in
                  SIGNOR: receptor signalling in general, not shown to be
                  active in the user's data
  your data       detection of the partners in the sender / receiver
                  populations (only if the expression matrix is prepared)
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import pandas as pd

from . import config, documentation, expression, index_store, reactome_reactions

from utils.biological_context import get_biological_context

CELLCHAT_DB = config.PROJECT_ROOT / "data" / "biological_context" / "cellchat" / "cellchatdb_human_interaction_local.tsv"
TYPES = {
    "Secreted Signaling": "secreted",
    "Cell-Cell Contact": "contact",
    "ECM-Receptor": "ecm",
    "Non-protein Signaling": "nonprotein",
}
LABELS = {
    "secreted": "Secreted ligand",
    "contact": "Cell–cell contact",
    "ecm": "Extracellular matrix – receptor",
    "nonprotein": "Non-protein (small-molecule) signalling",
    "unknown": "Mechanism not annotated",
}


@lru_cache(maxsize=1)
def _cellchatdb() -> pd.DataFrame:
    if not CELLCHAT_DB.exists():
        return pd.DataFrame()
    db = pd.read_csv(CELLCHAT_DB, sep="\t", dtype=str, keep_default_na=False)
    db["_rec_set"] = db["receptor.symbol"].map(lambda s: frozenset(x.strip().upper() for x in s.split(",") if x.strip()))
    db["_lig"] = db["ligand.symbol"].str.upper()
    return db


def _split(s: str) -> list[str]:
    return [x.strip() for x in str(s or "").replace(";", ",").split(",") if x.strip()]


def _classify(lig: dict, rcp: dict) -> dict[str, Any]:
    db = _cellchatdb()
    lig_genes = [g.upper() for g in lig["subunits"]] or [lig["label"].upper()]
    rec_set = frozenset(g.upper() for g in rcp["subunits"])
    if not db.empty:
        exact = db[db["_lig"].isin(lig_genes) & (db["_rec_set"] == rec_set)] if rec_set else db.iloc[0:0]
        if len(exact):
            r = exact.iloc[0]
            return {"type": TYPES.get(r["annotation"], "unknown"), "source": "CellChatDB (this ligand–receptor pair)",
                    "pathway": r["pathway_name"] or None, "evidence": r["evidence"] or None,
                    "co_stimulatory": _split(r["co_A_receptor"]), "co_inhibitory": _split(r["co_I_receptor"]),
                    "agonist": r["agonist"] or None, "antagonist": r["antagonist"] or None}
        same_lig = db[db["_lig"].isin(lig_genes)]
        if len(same_lig):
            ann = same_lig["annotation"].value_counts().index[0]
            return {"type": TYPES.get(ann, "unknown"), "source": "CellChatDB (other pairs with this ligand)",
                    "pathway": same_lig["pathway_name"].value_counts().index[0] or None, "evidence": None,
                    "co_stimulatory": [], "co_inhibitory": [], "agonist": None, "antagonist": None}
    return {"type": "unknown", "source": None, "pathway": None, "evidence": None,
            "co_stimulatory": [], "co_inhibitory": [], "agonist": None, "antagonist": None}


def mechanism(result_id: int) -> dict[str, Any] | None:
    rec = index_store.get_result(result_id)
    if not rec:
        return None
    lig = index_store.entity_display(rec["ligand_entity_key"], rec.get("m_ligand"))
    rcp = index_store.entity_display(rec["receptor_entity_key"], rec.get("m_receptor"))
    cls = _classify(lig, rcp)

    if lig["kind"] == "proxy":
        cls.update(type="nonprotein", source="CellPhoneDB (small molecule estimated from the enzyme that makes it)")
    if cls["type"] == "unknown":
        # Fall back on where UniProt locates the ligand protein.
        bc = get_biological_context(ligand=rec["ligand"], receptor=rec["receptor"],
                                    ligand_entity_key=rec["ligand_entity_key"], receptor_entity_key=rec["receptor_entity_key"])
        locs = " ".join(str(a.get("subcellular_location") or "") for a in bc.get("ligand_subunit_annotations") or [])
        low = locs.lower()
        if "extracellular matrix" in low:
            cls.update(type="ecm", source="UniProt subcellular location of the ligand")
        elif "secreted" in low:
            cls.update(type="secreted", source="UniProt subcellular location of the ligand")
        elif "cell membrane" in low or "membrane" in low:
            cls.update(type="contact", source="UniProt subcellular location of the ligand (membrane protein)")

    # Downstream: direct SIGNOR relations of each receptor subunit (receptor-general).
    steps = []
    for g in rcp["subunits"]:
        sg = documentation._signor_for(g)
        for rel in sg.get("direct", [])[:6]:
            if rel["direction"] not in ("activation", "inhibition"):
                continue
            steps.append({"from": g, "target": rel["target"], "direction": rel["direction"],
                          "mechanism": rel.get("mechanism"), "phenotype": (rel.get("target_type") or "").lower() == "phenotype",
                          "pmid": rel.get("pmid"), "signor_id": rel.get("signor_id")})
    seen, uniq = set(), []
    for s_ in steps:
        k = (s_["from"], s_["target"])
        if k not in seen:
            seen.add(k); uniq.append(s_)
    uniq = uniq[:5]

    data = None
    ev = expression.dataset_evidence(result_id)
    if ev and ev.get("available"):
        def best(side):
            sd = ev[side]
            if sd.get("joint"):
                return sd["joint"]["detected"]
            found = [x["detected"] for x in sd["subunits"] if x.get("found")]
            return min(found) if found else None
        data = {"ligand_detected": best("ligand"), "receptor_detected": best("receptor"),
                "dataset": ev.get("dataset")}

    try:
        reactome = reactome_reactions.reactions(rec, lig, rcp)
    except Exception as exc:      # a missing or damaged snapshot must not break the picture
        reactome = {"available": False, "message": f"Reactome reactions could not be read: {exc}"}

    return {
        "reactome": reactome,
        "sender": rec["sender"], "receiver": rec["receiver"],
        "ligand": lig, "receptor": rcp,
        "type": cls["type"], "type_label": LABELS[cls["type"]], "type_source": cls["source"],
        "pathway": cls["pathway"], "evidence": cls["evidence"],
        "co_stimulatory": cls["co_stimulatory"], "co_inhibitory": cls["co_inhibitory"],
        "signalling": uniq,
        "signalling_note": ("Direct downstream relations of the receptor in SIGNOR. They describe the receptor in "
                            "general, whatever ligand binds it, and are not shown to be active in your data."
                            if uniq else
                            "SIGNOR has no curated direct downstream relations for this receptor. Reactome pathways "
                            "that include it are listed below."),
        "your_data": data,
        "caution": ("The picture shows the known mechanism from databases. That these two populations "
                    "communicate this way in your sample is a computational prediction from RNA."),
    }
