"""
Curated Reactome reactions for a ligand-receptor pair.

A reaction qualifies when it involves the ligand AND at least one receptor
subunit (UniProt accessions, human). Reactions that include every receptor
subunit come first, then reactions whose name describes binding. The diagram
itself is served by Reactome (the browser loads it from reactome.org).
"""

from __future__ import annotations

import sqlite3
from contextlib import closing
from typing import Any

from . import config

from utils.biological_context import get_biological_context

DB = config.PROJECT_ROOT / "data" / "biological_context" / "reactome" / "reactome_reactions.sqlite"
EXPORT = "https://reactome.org/ContentService/exporter/reaction/{id}.{ext}"
DETAIL = "https://reactome.org/content/detail/{id}"


def _entries(annotations) -> dict[str, str]:
    """gene symbol -> UniProt accession"""
    out = {}
    for a in annotations or []:
        if a and a.get("entry") and str(a.get("found")) not in ("False", "false"):
            out[str(a.get("gene_symbol"))] = str(a["entry"])
    return out


def reactions(rec: dict[str, Any], lig: dict[str, Any], rcp: dict[str, Any], limit: int = 5) -> dict[str, Any]:
    if not DB.exists():
        return {"available": False, "message": ("Reaction-level Reactome data is not installed. Run "
                                                "“python web/tools/update_reactome_reactions.py” once (internet needed).")}
    if lig["kind"] == "proxy":
        return {"available": True, "items": [], "message": ("The ligand is a small molecule estimated from the enzyme "
                "that makes it; reactions of the enzyme would describe its synthesis, not receptor binding.")}
    bc = get_biological_context(ligand=rec["ligand"], receptor=rec["receptor"],
                                ligand_entity_key=rec["ligand_entity_key"], receptor_entity_key=rec["receptor_entity_key"])
    L, R = _entries(bc.get("ligand_subunit_annotations")), _entries(bc.get("receptor_subunit_annotations"))
    if not L or not R:
        return {"available": True, "items": [], "message": "The ligand or receptor could not be matched to UniProt."}
    with closing(sqlite3.connect(f"file:{DB.as_posix()}?mode=ro", uri=True)) as conn:
        def rx_for(acc):
            return {r[0]: r[1] for r in conn.execute(
                "SELECT reaction_id, reaction_name FROM uniprot_reactions WHERE uniprot_id = ?", (acc,))}
        lig_rx = {}
        for acc in L.values():
            lig_rx.update(rx_for(acc))
        rec_rx = {g: rx_for(acc) for g, acc in R.items()}
        created = conn.execute("SELECT value FROM metadata WHERE key='created'").fetchone()
    items = []
    for rid, name in lig_rx.items():
        covers = [g for g, rx in rec_rx.items() if rid in rx]
        if not covers:
            continue
        items.append({
            "id": rid, "name": name, "covers": covers,
            "missing": [g for g in R if g not in covers],
            "binding": any(w in name.lower() for w in ("bind", "binding", "complex", "interact")),
            "svg": EXPORT.format(id=rid, ext="svg"), "png": EXPORT.format(id=rid, ext="png"),
            "url": DETAIL.format(id=rid),
        })
    items.sort(key=lambda x: (len(x["missing"]), not x["binding"], len(x["name"])))
    return {
        "available": True, "items": items[:limit], "total": len(items),
        "snapshot": created[0] if created else None,
        "message": None if items else "No Reactome reaction involves both the ligand and the receptor.",
        "note": ("Curated human reaction from Reactome. It shows the molecular event; it does not show that it "
                 "happens in your sample. The diagram is loaded from reactome.org (internet needed)."),
    }
