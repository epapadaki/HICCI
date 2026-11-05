"""
Analysis records (provenance) written by tools/run_methods.py, and a draft
"Methods" paragraph generated from them for manuscripts.

The paragraph only states what the record contains; anything not recorded is
left out rather than guessed.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import config, datasets


def record(ds_id: str) -> dict[str, Any] | None:
    ds = datasets.get(ds_id)
    if not ds:
        return None
    rel = ds.meta.get("provenance")
    path = (config.PROJECT_ROOT / rel) if rel else None
    if not path or not path.exists():
        return {"available": False,
                "message": ("No analysis record. The CCI results of this dataset were produced "
                            "outside the platform, or before records were introduced. Running the "
                            "methods from the platform creates one.")}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {"available": False, "message": f"The analysis record could not be read: {exc}"}
    return {"available": True, "record": data, "methods_text": methods_text(data),
            "file": str(path.relative_to(config.PROJECT_ROOT)).replace("\\", "/")}


def _join(parts: list[str]) -> str:
    parts = [p for p in parts if p]
    return ", ".join(parts)


def methods_text(rec: dict[str, Any]) -> str:
    sw = rec.get("software", {})
    py = sw.get("python", {})
    r = sw.get("r", {})
    par = rec.get("parameters", {})
    db = rec.get("databases", {})
    cpdb_par = par.get("cellphonedb", {})
    cc = par.get("cellchat", {})

    cpdb = "CellPhoneDB" + (f" v{py['cellphonedb']}" if py.get("cellphonedb") else "")
    cpdb_details = _join([
        f"database {db.get('cellphonedb', {}).get('version')}" if db.get("cellphonedb", {}).get("version") else "",
        {"non statistical method": "non-statistical analysis",
         "statistical method": "statistical analysis"}.get(cpdb_par.get("analysis", ""), cpdb_par.get("analysis", "")),
        f"expression threshold {cpdb_par['expression_threshold']:g}" if cpdb_par.get("expression_threshold") is not None else "",
    ])
    cellchat = "CellChat" + (f" v{r['CellChat']}" if r.get("CellChat") else "")
    cc_details = _join([
        f"R {r['R']}" if r.get("R") else "",
        cc.get("database") or "",
        f"{cc['expression_average']} expression averaging" if cc.get("expression_average") else "",
        f"{cc['nboot']} permutations" if cc.get("nboot") else "",
        f"interactions retained at p < {cc['pvalue_threshold']}" if cc.get("pvalue_threshold") else "",
    ])
    natmi_par = par.get("natmi", {})
    natmi_details = _join([
        f"ligand–receptor resource {natmi_par['interaction_db']}" if natmi_par.get("interaction_db") else "",
    ])

    methods = [
        cpdb + (f" ({cpdb_details})" if cpdb_details else ""),
        "NATMI" + (f" ({natmi_details})" if natmi_details else ""),
        cellchat + (f" ({cc_details})" if cc_details else ""),
    ]
    date = (rec.get("finished") or rec.get("started") or "")[:10]
    text = (
        "Cell–cell communication was inferred with " + ", ".join(methods[:-1]) + " and " + methods[-1] + ". "
        "The methods were run through the CCC Platform"
        + (f" on {date}" if date else "")
        + (f" (Python {py['python']})" if py.get("python") else "")
        + ", and their predictions were harmonised into a complex-aware consensus in which each candidate "
          "interaction is defined by its sender population, receiver population, ligand entity and receptor "
          "entity, with multi-subunit complexes kept intact."
    )
    kb = rec.get("knowledge_bases", {})
    kb_names = [n for n, label in (("uniprot", "UniProtKB/Swiss-Prot"), ("reactome", "Reactome"), ("signor", "SIGNOR"))
                if kb.get(n, {}).get("exists", True) and kb.get(n)]
    if kb_names:
        labels = {"uniprot": "UniProtKB/Swiss-Prot", "reactome": "Reactome", "signor": "SIGNOR"}
        dates = [f"{labels[n]} (local file dated {kb[n]['modified'][:10]})" for n in kb_names if kb[n].get("modified")]
        if dates:
            text += (" Biological documentation of each candidate interaction used local copies of "
                     + ", ".join(dates[:-1]) + (" and " if len(dates) > 1 else "") + dates[-1]
                     + ", together with PubMed searches.")
    return text
