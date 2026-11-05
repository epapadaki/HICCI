"""
Run the CCI methods (CellPhoneDB, NATMI, CellChat) for one dataset and build
the complex-aware consensus that the web app reads.

The platform's own pipeline scripts (scripts/*.py, scripts/run_cellchat.R) are
used unchanged. They expect a fixed project layout, so each dataset gets its
own working copy of that layout:

    datasets/<id>/run/
        scripts/                       copies of the pipeline scripts
        data/cellphonedb/              the dataset's .h5ad and cell annotation
                                       (under the file names the scripts expect)
        data/cellphonedb/database/     CellPhoneDB database (shared, copied in)
        data/biological_context/cellchat/   CellChatDB snapshot (copied in)
        results/                       all method outputs

Run in the analysis environment:

    python web/tools/run_methods.py --id <dataset> --rscript "C:\\Program Files\\R\\R-4.6.1\\bin\\Rscript.exe"
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# Windows pipes default to the cp1252 code page, which cannot print many
# characters that appear in the methods' output. Always use UTF-8.
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rlauncher  # noqa: E402
import pylauncher  # noqa: E402
SCRIPTS = PROJECT_ROOT / "scripts"
NATMI = PROJECT_ROOT / "external" / "NATMI" / "ExtractEdges.py"
SHARED_CPDB_DB = PROJECT_ROOT / "data" / "cellphonedb" / "database" / "v5.0.0"
CELLCHAT_SNAPSHOT = PROJECT_ROOT / "data" / "biological_context" / "cellchat"
STEP_MARK = "@@STEP@@ "
RESULT_MARK = "@@RESULT@@ "

# File names the pipeline scripts expect.
H5AD_NAME = "endometrium_example_counts.h5ad"
META_NAME = "endometrium_example_meta.tsv"

PIPELINE_SCRIPTS = [
    "run_cellphonedb.py", "standardize_cellphonedb.py", "enrich_cellphonedb_standardized_complexes.py",
    "prepare_natmi_input.py", "standardize_natmi.py",
    "prepare_cellchat_input.py", "run_cellchat.R", "standardize_cellchat.py", "enrich_cellchat_raw_complexes.py",
    "build_complex_aware_candidate_consensus.py",
    "run_singlecellsignalr.R", "standardize_singlecellsignalr.py",
    "run_icellnet.R", "standardize_icellnet.py",
]
TOOLS = ["CellPhoneDB", "NATMI", "CellChat", "SingleCellSignalR", "ICELLNET"]
R_TOOLS = {"CellChat", "SingleCellSignalR", "ICELLNET"}
ICELLNET_DB = PROJECT_ROOT / "data" / "biological_context" / "icellnet" / "ICELLNETdb.tsv"

def consensus_rel(tools: list[str]) -> str:
    slug = "__".join(t.lower() for t in tools)
    return f"results/consensus_candidate_complex_aware/consensus_{slug}.csv"


# ---------------------------------------------------------------------------
# Provenance: what was run, with which versions, databases and parameters.
# Written to datasets/<id>/run/results/provenance.json, on success and on failure.
# ---------------------------------------------------------------------------
PROV: dict = {"schema": 1, "status": "running", "steps": [], "software": {}, "databases": {},
              "parameters": {}, "input": {}, "scripts": {}, "knowledge_bases": {}}
_PROV_PATH: list = []          # set once the run folder is known
_OUTPUT: dict = {}             # step output, for extracting parameters
_RSCRIPT: list = []


def sha256(path: Path, limit: int = 2 * 1024 ** 3) -> str | None:
    """SHA-256 of a file (skipped above 2 GB, where only size and date are recorded)."""
    import hashlib
    try:
        if path.stat().st_size > limit:
            return None
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def file_record(path: Path) -> dict:
    try:
        st = path.stat()
    except OSError:
        return {"path": str(path), "exists": False}
    return {"path": str(path), "size_bytes": st.st_size,
            "modified": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime)),
            "sha256": sha256(path)}


def python_versions() -> dict:
    import importlib.metadata as md
    import platform
    out = {"python": platform.python_version(), "platform": platform.platform()}
    for pkg in ("cellphonedb", "anndata", "pandas", "numpy", "scipy", "h5py"):
        try:
            out[pkg] = md.version(pkg)
        except md.PackageNotFoundError:
            pass
    return out


def r_versions(rscript: str) -> dict:
    code = ("v <- function(p) tryCatch(as.character(utils::packageVersion(p)), error=function(e) NA); "
            "cat(paste0('R=', R.version[['major']], '.', R.version[['minor']]), "
            "paste0('CellChat=', v('CellChat')), paste0('NMF=', v('NMF')), "
            "paste0('presto=', v('presto')), paste0('future=', v('future')), "
            "paste0('SingleCellSignalR=', v('SingleCellSignalR')), paste0('icellnet=', v('icellnet')), sep='\\n')")
    try:
        res = subprocess.run(rlauncher.command(rscript, ["-e", code]), capture_output=True,
                             text=True, encoding="utf-8", errors="replace", timeout=300)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"error": str(exc)}
    out = {}
    for line in res.stdout.splitlines():
        if "=" in line:
            k, _, v = line.strip().partition("=")
            if v and v != "NA":
                out[k] = v
    out["where"] = "WSL" if rlauncher.is_wsl(rscript) else "native"
    return out


def write_provenance() -> None:
    if not _PROV_PATH:
        return
    PROV["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    try:
        _PROV_PATH[0].parent.mkdir(parents=True, exist_ok=True)
        _PROV_PATH[0].write_text(json.dumps(PROV, indent=2, ensure_ascii=False), encoding="utf-8")
    except OSError as exc:
        print(f"  (could not write the provenance record: {exc})", flush=True)


def step(n: int, total: int, title: str) -> None:
    print(f"{STEP_MARK}{n}/{total} {title}", flush=True)
    print(f"\n=== Step {n} of {total}: {title} ===", flush=True)
    PROV["steps"].append({"step": n, "title": title, "started": time.strftime("%Y-%m-%d %H:%M:%S")})


def fail(message: str) -> None:
    print(f"\nERROR: {message}", flush=True)
    PROV["status"] = "failed"
    PROV["error"] = message
    if PROV["steps"]:
        PROV["failed_step"] = PROV["steps"][-1]["title"]
    if _RSCRIPT and "r" not in PROV["software"]:
        PROV["software"]["r"] = r_versions(_RSCRIPT[0])
    write_provenance()
    sys.exit(1)


def run(cmd: list[str], cwd: Path, what: str) -> None:
    env = dict(os.environ)
    # Scripts import the platform's utils package from the real project folder.
    env["PYTHONPATH"] = str(PROJECT_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    env["PYTHONUNBUFFERED"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"
    started = time.time()
    proc = subprocess.Popen(cmd, cwd=str(cwd), env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, encoding="utf-8", errors="replace", bufsize=1)
    captured = []
    for line in proc.stdout:
        captured.append(line)
        print("  " + line.rstrip(), flush=True)
    code = proc.wait()
    _OUTPUT[what] = "".join(captured[-400:])
    if PROV["steps"]:
        PROV["steps"][-1]["seconds"] = round(time.time() - started, 1)
    if code != 0:
        fail(f"{what} stopped with exit code {code}. The lines above show the reason.")
    print(f"  ({what}: {time.time() - started:.0f} s)", flush=True)


def link_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        dst.unlink()
    try:
        os.link(src, dst)          # no extra disk space; same drive only
    except OSError:
        shutil.copy2(src, dst)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--id", required=True)
    ap.add_argument("--rscript", help="Path to Rscript (needed for CellChat and SingleCellSignalR)")
    ap.add_argument("--tools-python", dest="tools_python",
                    help="Python for CellPhoneDB (optional). A 'wsl:<linux path>' value runs it "
                         "inside WSL — use this when a locked-down Windows blocks SciPy's DLLs. "
                         "If omitted, CellPhoneDB uses the same Python as the rest of the pipeline.")
    ap.add_argument("--methods", default="CellPhoneDB,NATMI,CellChat",
                    help="Comma-separated tools to run: " + ", ".join(TOOLS))
    ap.add_argument("--cores", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    a = ap.parse_args(argv)

    import pandas as pd

    folder = PROJECT_ROOT / "datasets" / a.id
    ds_file = folder / "dataset.json"
    if not ds_file.exists():
        fail(f"Dataset '{a.id}' not found ({ds_file}).")
    meta = json.loads(ds_file.read_text(encoding="utf-8"))
    organism = str(meta.get("organism") or "human").lower()
    if organism not in ("human", "homo sapiens", "hsa"):
        fail(f"The methods are set up for human data (CellPhoneDB v5, CellChatDB human); this dataset is '{organism}'.")
    src_rel = (meta.get("source") or {}).get("counts")
    if not src_rel:
        fail("The dataset has no source .h5ad recorded. Prepare it on the Datasets page first.")
    h5ad = Path(src_rel) if Path(src_rel).is_absolute() else PROJECT_ROOT / src_rel
    cells_tsv = folder / "expression" / "cells.tsv"
    if not h5ad.exists():
        fail(f"The .h5ad file was not found: {h5ad}")
    if not cells_tsv.exists():
        fail("The dataset's expression matrix has not been prepared yet.")
    tools = [t for t in TOOLS if t in [x.strip() for x in a.methods.split(",")]]
    if not tools:
        fail("Choose at least one tool.")
    if R_TOOLS & set(tools) and (not a.rscript or not rlauncher.exists(a.rscript)):
        fail("Rscript was not found. Set it on the Datasets page (needed for " + " and ".join(sorted(R_TOOLS & set(tools))) + ").")
    for s in PIPELINE_SCRIPTS:
        if not (SCRIPTS / s).exists():
            fail(f"Pipeline script missing: scripts/{s}")
    if "NATMI" in tools and not NATMI.exists():
        fail(f"NATMI was not found at {NATMI}")
    if "CellPhoneDB" in tools and a.tools_python and not pylauncher.exists(a.tools_python):
        fail(f"The CellPhoneDB Python was not found: {a.tools_python}")

    # Steps: setup, the chosen tools, consensus, registration.
    # R tools share one prepared expression input (counted once).
    n_steps = {"CellPhoneDB": 3, "NATMI": 3, "CellChat": 3, "SingleCellSignalR": 2, "ICELLNET": 2}
    total = 1 + sum(n_steps[t] for t in tools) + (1 if R_TOOLS & set(tools) else 0) + 2
    slug = "__".join(t.lower() for t in tools)
    run_dir = folder / "runs" / slug
    py = sys.executable
    _RSCRIPT.append(a.rscript)
    PROV.update(dataset=a.id, dataset_label=meta.get("label"), tools=tools, started=time.strftime("%Y-%m-%d %H:%M:%S"),
                platform_tool="CCC Platform, web/tools/run_methods.py")
    PROV["software"]["python"] = python_versions()
    PROV["input"] = {"h5ad": file_record(h5ad), "cells_tsv": file_record(cells_tsv),
                     "cell_type_column": meta.get("cell_type_column")}
    PROV["scripts"] = {s: sha256(SCRIPTS / s) for s in PIPELINE_SCRIPTS}
    PROV["scripts"]["NATMI/ExtractEdges.py"] = sha256(NATMI)

    step(1, total, "Setting up the working folder")
    if (run_dir / "results").exists():
        stamp = time.strftime("%Y%m%d-%H%M%S")
        keep = folder / "runs" / f"{slug}_previous_{stamp}"
        keep.mkdir(parents=True, exist_ok=True)
        shutil.move(str(run_dir / "results"), str(keep / "results"))
        print(f"  Earlier results of the same tools kept in runs/{slug}_previous_{stamp}/")
    _PROV_PATH.append(run_dir / "results" / "provenance.json")
    shutil.rmtree(run_dir / "scripts", ignore_errors=True)
    (run_dir / "scripts").mkdir(parents=True, exist_ok=True)
    for s in PIPELINE_SCRIPTS:
        shutil.copy2(SCRIPTS / s, run_dir / "scripts" / s)
    cpdb_dir = run_dir / "data" / "cellphonedb"
    link_or_copy(h5ad, cpdb_dir / H5AD_NAME)
    cells = pd.read_csv(cells_tsv, sep="\t", dtype=str, keep_default_na=False)
    cells[["cell", "cell_type"]].rename(columns={"cell": "Cell"}).to_csv(cpdb_dir / META_NAME, sep="\t", index=False)
    print(f"  {len(cells)} cells, {cells['cell_type'].nunique()} cell types")
    if SHARED_CPDB_DB.exists():
        shutil.copytree(SHARED_CPDB_DB, cpdb_dir / "database" / "v5.0.0", dirs_exist_ok=True)
        print("  CellPhoneDB database: local copy")
    else:
        print("  CellPhoneDB database: will be downloaded (internet needed)")
    if CELLCHAT_SNAPSHOT.exists():
        dst = run_dir / "data" / "biological_context" / "cellchat"
        dst.mkdir(parents=True, exist_ok=True)
        for f in CELLCHAT_SNAPSHOT.glob("cellchatdb_*"):
            shutil.copy2(f, dst / f.name)
    (run_dir / "results" / "natmi").mkdir(parents=True, exist_ok=True)

    sc = run_dir / "scripts"
    n = [1]
    def nxt(title):
        n[0] += 1
        step(n[0], total, title)
    wsl = " (WSL)" if a.rscript and rlauncher.is_wsl(a.rscript) else ""

    if "CellPhoneDB" in tools:
        # CellPhoneDB can run under a different Python (e.g. one inside WSL), so a
        # locked-down Windows that blocks SciPy's DLLs can still use it. When no
        # tools-python is given, it uses the pipeline Python, exactly as before.
        cpy = a.tools_python or py
        cwsl = " (WSL)" if a.tools_python and pylauncher.is_wsl(a.tools_python) else ""
        nxt("CellPhoneDB: running" + cwsl);      run(pylauncher.command(cpy, [str(sc / "run_cellphonedb.py")], cwd=run_dir), run_dir, "CellPhoneDB")
        nxt("CellPhoneDB: standardising");       run(pylauncher.command(cpy, [str(sc / "standardize_cellphonedb.py")], cwd=run_dir), run_dir, "CellPhoneDB standardisation")
        nxt("CellPhoneDB: resolving complexes"); run(pylauncher.command(cpy, [str(sc / "enrich_cellphonedb_standardized_complexes.py")], cwd=run_dir), run_dir, "CellPhoneDB complexes")

    if "NATMI" in tools:
        nxt("NATMI: preparing input");           run([py, str(sc / "prepare_natmi_input.py")], run_dir, "NATMI input")
        natmi_in = run_dir / "data" / "natmi"
        nxt("NATMI: running")
        run([py, str(NATMI), "--interDB", "lrc2p", "--species", "human",
             "--emFile", str(natmi_in / "endometrium_expression.csv"),
             "--annFile", str(natmi_in / "endometrium_annotations.csv"),
             "--coreNum", str(a.cores), "--out", str(run_dir / "results" / "natmi")],
            NATMI.parent, "NATMI")
        PROV["parameters"]["natmi"] = {"interaction_db": "lrc2p", "species": "human", "cores": a.cores,
                                       "note": "bundled NATMI, see external/NATMI/PATCHES_CCC_PLATFORM.md"}
        nxt("NATMI: standardising");             run([py, str(sc / "standardize_natmi.py")], run_dir, "NATMI standardisation")

    if R_TOOLS & set(tools):
        # Normalised expression and labels, shared by the R tools.
        nxt("Preparing the R input (" + ", ".join(t for t in tools if t in R_TOOLS) + ")")
        run([py, str(sc / "prepare_cellchat_input.py")], run_dir, "R input")

    if "CellChat" in tools:
        nxt("CellChat: running in R" + wsl)
        run(rlauncher.command(a.rscript, ["scripts/run_cellchat.R"], cwd=run_dir), run_dir, "CellChat")
        nxt("CellChat: standardising");          run([py, str(sc / "standardize_cellchat.py")], run_dir, "CellChat standardisation")
        nxt("CellChat: resolving complexes");    run([py, str(sc / "enrich_cellchat_raw_complexes.py")], run_dir, "CellChat complexes")

    if "SingleCellSignalR" in tools:
        nxt("SingleCellSignalR: running in R" + wsl)
        run(rlauncher.command(a.rscript, ["scripts/run_singlecellsignalr.R"], cwd=run_dir), run_dir, "SingleCellSignalR")
        nxt("SingleCellSignalR: standardising"); run([py, str(sc / "standardize_singlecellsignalr.py")], run_dir, "SingleCellSignalR standardisation")

    if "ICELLNET" in tools:
        if ICELLNET_DB.exists():
            (run_dir / "data" / "icellnet").mkdir(parents=True, exist_ok=True)
            shutil.copy2(ICELLNET_DB, run_dir / "data" / "icellnet" / "ICELLNETdb.tsv")
        nxt("ICELLNET: running in R" + wsl)
        run(rlauncher.command(a.rscript, ["scripts/run_icellnet.R"], cwd=run_dir), run_dir, "ICELLNET")
        nxt("ICELLNET: standardising");          run([py, str(sc / "standardize_icellnet.py")], run_dir, "ICELLNET standardisation")
        PROV["databases"]["icellnet"] = file_record(run_dir / "data" / "icellnet" / "ICELLNETdb.tsv")

    nxt("Combining the tools’ results")
    run([py, str(sc / "build_complex_aware_candidate_consensus.py"), "--write", "--dataset-id", a.id,
         "--methods", ",".join(tools)], run_dir, "Consensus")

    nxt("Registering the results")
    consensus = run_dir / consensus_rel(tools)
    if not consensus.exists():
        fail(f"The consensus file was not produced: {consensus}")
    rel = lambda p: str(p.relative_to(PROJECT_ROOT)).replace("\\", "/")
    entry = {"consensus": rel(consensus), "tools": tools, "created": time.strftime("%Y-%m-%d %H:%M:%S")}
    if "CellChat" in tools:
        entry["cellchat"] = rel(run_dir / "results" / "cellchat" / "standardized_interactions.csv")
        entry["cellchat_run_info"] = rel(run_dir / "results" / "cellchat" / "run_info.csv")
    # Every set of tools keeps its own results; the latest run becomes the one shown.
    runs = meta.setdefault("runs", {})
    if "results" in meta and not runs and meta["results"].get("consensus"):
        old = dict(meta["results"]); runs[old.get("run", "earlier")] = old      # keep results made before runs existed
    entry["run"] = slug
    runs[slug] = entry
    meta["active_run"] = slug
    meta["results"] = dict(entry)
    meta["results_created"] = entry["created"]

    # ---- provenance: parameters, databases, R side
    import re
    m = re.search(r"\[(Non Statistical Method|Statistical Method)\] Threshold:([\d.]+)", _OUTPUT.get("CellPhoneDB", ""))
    PROV["parameters"]["cellphonedb"] = ({"analysis": m.group(1).lower(), "expression_threshold": float(m.group(2))}
                                         if m else {"note": "parameters as set in scripts/run_cellphonedb.py"})
    run_info = run_dir / "results" / "cellchat" / "run_info.csv"
    if run_info.exists():
        import csv as _csv
        with open(run_info, newline="", encoding="utf-8") as fh:
            PROV["parameters"]["cellchat"] = {r["field"]: r["value"] for r in _csv.DictReader(fh)}
    for _tool, _dir in (("icellnet", "icellnet"),):
        _info = run_dir / "results" / _dir / "run_info.csv"
        if _info.exists():
            import csv as _csv3
            with open(_info, newline="", encoding="utf-8") as fh:
                PROV["parameters"][_tool] = {r["field"]: r["value"] for r in _csv3.DictReader(fh)}
    scsr_info = run_dir / "results" / "singlecellsignalr" / "run_info.csv"
    if scsr_info.exists():
        import csv as _csv2
        with open(scsr_info, newline="", encoding="utf-8") as fh:
            PROV["parameters"]["singlecellsignalr"] = {r["field"]: r["value"] for r in _csv2.DictReader(fh)}
    cpdb_db = cpdb_dir / "database" / "v5.0.0"
    if "CellPhoneDB" in tools:
        PROV["databases"]["cellphonedb"] = {"version": "v5.0.0", **file_record(cpdb_db / "cellphonedb.zip")}
    PROV["databases"]["cellchatdb"] = {"name": PROV["parameters"].get("cellchat", {}).get("database")}
    if "NATMI" in tools:
        PROV["databases"]["natmi_lrc2p"] = file_record(NATMI.parent / "lrdbs" / "lrc2p.csv")
    if "SingleCellSignalR" in tools:
        PROV["databases"]["singlecellsignalr_lrdb"] = {"name": "LRdb bundled with SingleCellSignalR v1"}
    kb = PROJECT_ROOT / "data" / "biological_context"
    for name, rel_path in (("uniprot", "uniprot_human_reviewed.tsv"), ("reactome", "reactome/reactome_local.sqlite"),
                           ("signor", "signor/signor_local.sqlite")):
        PROV["knowledge_bases"][name] = file_record(kb / rel_path)
    PROV["software"]["r"] = r_versions(a.rscript)
    PROV["status"] = "finished"
    PROV["results"] = dict(meta["results"])
    write_provenance()
    meta["provenance"] = rel(_PROV_PATH[0])
    meta["runs"][slug]["provenance"] = meta["provenance"]
    meta["results"]["provenance"] = meta["provenance"]
    ds_file.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    if not SHARED_CPDB_DB.exists() and (cpdb_dir / "database" / "v5.0.0").exists():
        shutil.copytree(cpdb_dir / "database" / "v5.0.0", SHARED_CPDB_DB)   # reuse for next datasets
    print(RESULT_MARK + json.dumps({"consensus": meta["results"]["consensus"]}), flush=True)
    print("\nAll tools finished. The results are available in the web app.", flush=True)


if __name__ == "__main__":
    main()
