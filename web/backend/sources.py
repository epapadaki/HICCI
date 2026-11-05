"""
Versions and dates of the local knowledge bases and tool databases, for the
Help page (and for anyone citing the platform).
"""

from __future__ import annotations

import sqlite3
import time
from contextlib import closing
from pathlib import Path
from typing import Any

from . import config

BC = config.PROJECT_ROOT / "data" / "biological_context"


def _date(p: Path) -> str | None:
    try:
        return time.strftime("%Y-%m-%d", time.localtime(p.stat().st_mtime))
    except OSError:
        return None


def _meta(db: Path) -> dict[str, str]:
    if not db.exists():
        return {}
    try:
        with closing(sqlite3.connect(f"file:{db.as_posix()}?mode=ro", uri=True)) as c:
            return {str(k): str(v) for k, v in c.execute("SELECT key, value FROM metadata")}
    except sqlite3.Error:
        return {}


def sources() -> list[dict[str, Any]]:
    out = []
    uni = BC / "uniprot_human_reviewed.tsv"
    out.append({"name": "UniProtKB/Swiss-Prot (human, reviewed)", "use": "Function and location of ligands and receptors",
                "url": "https://www.uniprot.org", "local": uni.exists(), "date": _date(uni), "version": None})
    rx = BC / "reactome" / "reactome_local.sqlite"; m = _meta(rx)
    out.append({"name": "Reactome (pathways)", "use": "Pathways that include each receptor",
                "url": "https://reactome.org", "local": rx.exists(), "date": (m.get("imported_at_utc") or "")[:10] or _date(rx),
                "version": None})
    rr = BC / "reactome" / "reactome_reactions.sqlite"; m = _meta(rr)
    out.append({"name": "Reactome (reactions)", "use": "Curated reaction diagrams of ligand–receptor pairs",
                "url": "https://reactome.org", "local": rr.exists(), "date": (m.get("created") or "")[:10] or _date(rr), "version": None})
    sg = BC / "signor" / "signor_local.sqlite"; m = _meta(sg)
    out.append({"name": "SIGNOR", "use": "Direct downstream relations of receptors",
                "url": "https://signor.uniroma2.it", "local": sg.exists(), "date": (m.get("imported_at_utc") or "")[:10] or _date(sg),
                "version": (m.get("source_file") or "").replace("_release.txt", " release").replace(".txt", "") or None})
    cc = BC / "cellchat" / "cellchatdb_human_interaction_local.tsv"
    out.append({"name": "CellChatDB (human)", "use": "Interaction type (secreted, contact, ECM) and co-receptors; database of CellChat",
                "url": "https://github.com/jinworks/CellChat", "local": cc.exists(), "date": _date(cc), "version": None})
    out.append({"name": "CellPhoneDB database v5.0.0", "use": "Database of CellPhoneDB (downloaded on the first run)",
                "url": "https://www.cellphonedb.org", "local": (config.PROJECT_ROOT / "data" / "cellphonedb" / "database" / "v5.0.0").exists(),
                "date": None, "version": "v5.0.0"})
    out.append({"name": "NATMI lrc2p", "use": "Ligand–receptor list used by NATMI",
                "url": "https://github.com/forrest-lab/NATMI", "local": (config.PROJECT_ROOT / "external" / "NATMI" / "lrdbs" / "lrc2p.csv").exists(),
                "date": None, "version": None})
    out.append({"name": "SingleCellSignalR LRdb", "use": "Ligand–receptor list used by SingleCellSignalR (bundled with the R package, v1)",
                "url": "https://bioconductor.org/packages/SingleCellSignalR", "local": None, "date": None, "version": "v1"})
    ic = BC / "icellnet" / "ICELLNETdb.tsv"
    out.append({"name": "ICELLNET database", "use": "Curated ligand–receptor list used by ICELLNET (human, experimentally demonstrated)",
                "url": "https://github.com/soumelis-lab/ICELLNET", "local": ic.exists(), "date": _date(ic), "version": None})
    rs = BC / "response"
    out.append({"name": "PROGENy (via OmniPath)", "use": "Pathway-responsive genes for the receiver response (14 pathways)",
                "url": "https://saezlab.github.io/progeny/", "local": (rs / "progeny_human_top500.tsv").exists(),
                "date": _date(rs / "progeny_human_top500.tsv"), "version": None})
    out.append({"name": "CollecTRI (via OmniPath)", "use": "Transcription factor targets for the receiver response",
                "url": "https://github.com/saezlab/CollecTRI", "local": (rs / "collectri_human.tsv").exists(),
                "date": _date(rs / "collectri_human.tsv"), "version": None})
    out.append({"name": "PubMed (NCBI E-utilities)", "use": "Publications on each ligand–receptor pair (searched online, then cached)",
                "url": "https://pubmed.ncbi.nlm.nih.gov", "local": True, "date": None, "version": None})
    return out
