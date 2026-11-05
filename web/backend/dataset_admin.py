"""
Adding datasets from the web app.

Inspection and preparation run tools/prepare_dataset.py in the analysis
environment (settings.analysis_python), so the web app never needs anndata.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path
from typing import Any

from . import config, datasets, jobs, settings

TOOL = config.WEB_ROOT / "tools" / "prepare_dataset.py"
RUN_TOOL = config.WEB_ROOT / "tools" / "run_methods.py"
ID_RE = re.compile(r"^[a-z0-9][a-z0-9_-]{1,40}$")


def _resolve(path: str) -> Path:
    p = Path(settings.clean_path(path))
    return p if p.is_absolute() else (config.PROJECT_ROOT / p)


def _python() -> str:
    exe = settings.analysis_python()
    if not exe or not Path(exe).exists():
        raise ValueError("Set the analysis environment first (the Python that has anndata installed).")
    return exe


def suggest_id(label: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "_", (label or "").lower()).strip("_")
    return (s or "dataset")[:40]


def inspect(h5ad: str, metadata: str | None = None) -> dict[str, Any]:
    h5 = _resolve(h5ad)
    if not h5.exists():
        raise ValueError(f"File not found: {h5}")
    if h5.suffix.lower() != ".h5ad":
        raise ValueError("The file must be an .h5ad file.")
    cmd = [_python(), str(TOOL), "--inspect", "--h5ad", str(h5)]
    if metadata:
        meta = _resolve(metadata)
        if not meta.exists():
            raise ValueError(f"File not found: {meta}")
        cmd += ["--metadata", str(meta)]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=900, cwd=str(config.PROJECT_ROOT))
    except subprocess.TimeoutExpired:
        raise ValueError("Reading the file took too long.")
    for line in res.stdout.splitlines():
        if line.startswith(jobs.RESULT_MARK):
            out = json.loads(line[len(jobs.RESULT_MARK):])
            out["h5ad"] = str(h5)
            return out
    detail = (res.stderr or res.stdout).strip().splitlines()[-1:] or ["no output"]
    raise ValueError("The file could not be read: " + detail[0])


def prepare(form: dict[str, Any]) -> dict[str, Any]:
    label = str(form.get("label") or "").strip()
    ds_id = str(form.get("id") or suggest_id(label)).strip().lower()
    if not label:
        raise ValueError("Give the dataset a name.")
    if not ID_RE.match(ds_id):
        raise ValueError("The identifier may contain lowercase letters, digits, - and _ (2–41 characters).")
    existing = datasets.get(ds_id)
    if existing and not form.get("replace"):
        raise ValueError(f"A dataset with the identifier “{ds_id}” already exists.")
    if not form.get("cell_type_column"):
        raise ValueError("Choose the column with the cell types.")
    if jobs.running("prepare"):
        raise ValueError("Another dataset is being prepared. Wait for it to finish.")

    h5 = _resolve(form.get("h5ad") or "")
    cmd = [_python(), str(TOOL), "--json", "--id", ds_id, "--h5ad", str(h5),
           "--cell-type-column", str(form["cell_type_column"]), "--label", label]
    for key, flag in (("metadata", "--metadata"), ("organism", "--organism"), ("tissue", "--tissue"),
                      ("sample_column", "--sample-column"), ("group_column", "--group-column"),
                      ("layer", "--layer"), ("gene_symbol_column", "--gene-symbol-column")):
        val = str(form.get(key) or "").strip()
        if val:
            cmd += [flag, str(_resolve(val)) if key == "metadata" else val]
    if form.get("use_raw"):
        cmd.append("--use-raw")
    return jobs.start("prepare", f"Preparing “{label}”", cmd, cwd=str(config.PROJECT_ROOT))


def details(ds_id: str) -> dict[str, Any] | None:
    ds = datasets.get(ds_id)
    if not ds:
        return None
    report = ds.meta.get("validation") or {}
    return {**ds.summary(), "validation": report, "source": ds.meta.get("source") or {},
            "cell_type_column": ds.meta.get("cell_type_column"),
            "sample_column": ds.meta.get("sample_column"),
            "group_column": ds.meta.get("group_column")}


def run_methods(ds_id: str, tools: list[str] | None = None) -> dict[str, Any]:
    ds = datasets.get(ds_id)
    if not ds:
        raise ValueError("Dataset not found.")
    if ds.expression_dir is None:
        raise ValueError("Prepare the dataset’s expression matrix first.")
    if jobs.running("methods"):
        raise ValueError("CCI tools are already running for a dataset. Wait for them to finish.")
    known = ["CellPhoneDB", "NATMI", "CellChat", "SingleCellSignalR", "ICELLNET"]
    tools = [t for t in known if t in (tools or known[:3])]
    if not tools:
        raise ValueError("Choose at least one tool.")
    check = settings.check_methods()
    missing = [t for t in tools if not check["tools"][t]["ok"]]
    if missing:
        raise ValueError("Not ready to run: " + "; ".join(f"{t}: {check['tools'][t]['message']}" for t in missing))
    cmd = [_python(), str(RUN_TOOL), "--id", ds.id, "--methods", ",".join(tools)]
    if check.get("rscript"):
        cmd += ["--rscript", check["rscript"]]
    if "CellPhoneDB" in tools:
        tp = settings.tools_python()
        if tp:
            cmd += ["--tools-python", tp]
    return jobs.start("methods", f"Running {', '.join(tools)} for “{ds.label}”", cmd,
                      cwd=str(config.PROJECT_ROOT), dataset=ds.id)


def running_jobs() -> list[dict[str, Any]]:
    return jobs.running()


def discover_files(limit: int = 40) -> dict[str, Any]:
    """Single-cell files (and possible annotation files) inside the project's data folder.

    Files that are already the source of an existing dataset (including the
    shipped example) are skipped, so the upload form only offers new files.
    """
    from . import datasets as _datasets  # late import to avoid a cycle
    used: set[str] = set()
    for d in _datasets.all_datasets():
        src = (d.meta or {}).get("source") or {}
        for key in ("counts", "metadata"):
            val = src.get(key)
            if val:
                used.add(str(val).replace("\\", "/").lstrip("./"))

    root = config.PROJECT_ROOT / "data"
    h5ad, annot = [], []
    if root.exists():
        for f in root.rglob("*"):
            if not f.is_file():
                continue
            rel = str(f.relative_to(config.PROJECT_ROOT)).replace("\\", "/")
            if rel.lstrip("./") in used:
                continue  # already imported as a dataset
            if f.suffix.lower() == ".h5ad":
                st = f.stat()
                h5ad.append({"path": rel, "name": f.name, "folder": str(f.parent.relative_to(config.PROJECT_ROOT)).replace("\\", "/"),
                             "size_mb": round(st.st_size / 1e6, 1)})
            elif f.suffix.lower() in (".tsv", ".csv", ".txt") and any(k in f.name.lower() for k in ("meta", "annot", "obs", "cell_type", "celltype")):
                annot.append({"path": rel, "name": f.name, "folder": str(f.parent.relative_to(config.PROJECT_ROOT)).replace("\\", "/")})
    h5ad.sort(key=lambda x: x["path"]); annot.sort(key=lambda x: x["path"])
    return {"h5ad": h5ad[:limit], "annotation": annot[:limit], "folder": "data"}
