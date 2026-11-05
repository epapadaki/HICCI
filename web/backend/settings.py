"""
Local settings of the web app (web/workspace/settings.json).

analysis_python: the Python interpreter of the analysis environment, the one
that has anndata/scanpy and the CCI methods installed. The web app itself
stays light and runs heavy steps there, as separate processes.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any

from . import config

import sys as _sys
_sys.path.insert(0, str(config.WEB_ROOT / "tools"))
import rlauncher  # noqa: E402
import pylauncher  # noqa: E402

SETTINGS_FILE = config.WORKSPACE_DIR / "settings.json"


def load() -> dict[str, Any]:
    if SETTINGS_FILE.exists():
        try:
            # utf-8-sig tolerates a UTF-8 BOM (editors on Windows add one); a BOM
            # under plain utf-8 would break json parsing and silently reset every
            # setting to its default — including ephemeral_mode, which would then
            # wipe all datasets on the next launch.
            return json.loads(SETTINGS_FILE.read_text(encoding="utf-8-sig"))
        except (OSError, json.JSONDecodeError):
            pass
    return {}


def save(values: dict[str, Any]) -> dict[str, Any]:
    data = {**load(), **values}
    SETTINGS_FILE.write_text(json.dumps(data, indent=2), encoding="utf-8")
    # A changed Python/Rscript invalidates cached tool checks.
    if any(k in values for k in ("analysis_python", "rscript", "tools_python")):
        try:
            _CHECK_CACHE.clear()
        except NameError:
            pass
    return data


def _candidates() -> list[Path]:
    home = Path.home()
    names = ["ccc-platform"]
    out = []
    for base in ("miniconda3", "anaconda3", "miniforge3", "mambaforge"):
        for name in names:
            out.append(home / base / "envs" / name / "python.exe")
            out.append(home / base / "envs" / name / "bin" / "python")
    return out


def analysis_python() -> str | None:
    """Configured path, else CCC_ANALYSIS_PYTHON, else a conda env named ccc-platform."""
    # The launcher (start.sh) sets CCC_ANALYSIS_PYTHON; it wins over saved settings,
    # which may belong to another way of running the platform.
    configured = os.environ.get("CCC_ANALYSIS_PYTHON") or load().get("analysis_python")
    if configured:
        return configured
    import importlib.util, sys
    if importlib.util.find_spec("anndata") is not None:
        return sys.executable          # the app's own environment has everything
    for c in _candidates():
        if c.exists():
            return str(c)
    return None


def clean_path(value: str) -> str:
    """Accept paths pasted with quotes or with forward/back slashes."""
    return str(value or "").strip().strip('"').strip("'").strip()


def test_analysis_python(path: str | None = None) -> dict[str, Any]:
    exe = clean_path(path or analysis_python() or "")
    if not exe:
        return {"ok": False, "message": "No analysis environment is set."}
    if not Path(exe).exists():
        return {"ok": False, "path": exe, "message": f"File not found: {exe}"}
    code = ("import sys, anndata, numpy, pandas; "
            "print(sys.version.split()[0], anndata.__version__)")
    try:
        res = subprocess.run([exe, "-c", code], capture_output=True, text=True, timeout=180)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "path": exe, "message": f"Could not start it: {exc}"}
    if res.returncode != 0:
        last = (res.stderr or res.stdout).strip().splitlines()[-1:] or ["unknown error"]
        return {"ok": False, "path": exe,
                "message": "It started, but anndata could not be loaded: " + last[0]}
    py, adv = (res.stdout.strip().split() + ["?", "?"])[:2]
    return {"ok": True, "path": exe, "message": f"Python {py}, anndata {adv}"}


def tools_python() -> str | None:
    """Optional separate Python for the Python CCI tools (CellPhoneDB). A
    ``wsl:<linux path>`` value runs CellPhoneDB inside WSL, which lets a
    locked-down Windows (one that blocks SciPy's DLLs) still use it. When unset,
    CellPhoneDB uses the ordinary analysis Python."""
    return os.environ.get("CCC_TOOLS_PYTHON") or load().get("tools_python") or None


def test_tools_python(path: str | None = None) -> dict[str, Any]:
    """Can CellPhoneDB be imported by the given Python? Supports a WSL Python."""
    exe = clean_path(path if path is not None else (tools_python() or ""))
    if not exe:
        return {"ok": False, "message": "No CellPhoneDB Python is set."}
    if not pylauncher.exists(exe):
        return {"ok": False, "path": exe, "message": f"Not found: {exe}"}
    code = "import cellphonedb, importlib.metadata as m; print(m.version('cellphonedb'))"
    try:
        res = subprocess.run(pylauncher.command(exe, ["-c", code]),
                             capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "path": exe, "message": f"Could not start it: {exc}"}
    if res.returncode != 0:
        last = (res.stderr or res.stdout).strip().splitlines()[-1:] or ["unknown error"]
        return {"ok": False, "path": exe,
                "message": "It started, but CellPhoneDB could not be imported: " + last[0]}
    return {"ok": True, "path": exe, "message": f"CellPhoneDB {res.stdout.strip()} ({'WSL' if pylauncher.is_wsl(exe) else 'native'})"}


def rscript() -> str | None:
    """Configured Rscript, else Rscript on PATH, else the newest R under Program Files."""
    configured = os.environ.get("CCC_RSCRIPT") or load().get("rscript")
    if configured:
        return configured
    import shutil
    found = shutil.which("Rscript")
    if found:
        return found
    for base in (Path("C:/Program Files/R"), Path("C:/Program Files (x86)/R")):
        if base.exists():
            versions = sorted(base.glob("R-*"), reverse=True)
            for v in versions:
                exe = v / "bin" / "Rscript.exe"
                if exe.exists():
                    return str(exe)
    return None


def _run(cmd: list[str], timeout: int = 240) -> tuple[bool, str]:
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, str(exc)
    text = ((res.stdout or "") + "\n" + (res.stderr or "")).strip()
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if res.returncode != 0:
        errors = [l for l in lines if "error" in l.lower() and "execution halted" not in l.lower()]
        return False, (errors[0] if errors else (lines[-1] if lines else "unknown error"))
    return True, (lines[-1] if lines else "")


_CELLCHAT_TEST = "suppressPackageStartupMessages(library(CellChat)); cat(as.character(packageVersion('CellChat')))"


def _check_cellchat(r: str | None) -> dict[str, Any]:
    if not r or not rlauncher.exists(r):
        return {"ok": False, "message": "Rscript not found. Give its path below."}
    ok, last = _run(rlauncher.command(r, ["-e", _CELLCHAT_TEST]), timeout=300)
    where = " (in WSL)" if rlauncher.is_wsl(r) else ""
    return {"ok": ok, "message": f"CellChat {last}{where}" if ok else "R found, but CellChat could not be loaded: " + last}


def _wsl_rscripts() -> list[str]:
    import shutil
    if not shutil.which("wsl"):
        return []
    try:
        res = subprocess.run(["wsl", "bash", "-lc",
                              "ls -d ~/miniforge3/envs/*/bin/Rscript ~/miniconda3/envs/*/bin/Rscript "
                              "~/anaconda3/envs/*/bin/Rscript ~/mambaforge/envs/*/bin/Rscript 2>/dev/null; "
                              "command -v Rscript"],
                             capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return []
    paths = [l.strip() for l in res.stdout.splitlines() if l.strip().startswith("/")]
    paths.sort(key=lambda p: ("cellchat" not in p.lower(), p))      # likely candidates first
    return [rlauncher.WSL_PREFIX + p for p in dict.fromkeys(paths)]


_TOOL_KEYS = {"CellPhoneDB": "cellphonedb", "NATMI": "natmi", "CellChat": "cellchat",
              "SingleCellSignalR": "singlecellsignalr", "ICELLNET": "icellnet"}

# In-memory cache of per-tool check results, so a tool is probed once per session
# (loading an R package can take 30-60s). Cleared when the environment changes.
_CHECK_CACHE: dict[str, Any] = {}


def clear_check_cache() -> None:
    _CHECK_CACHE.clear()


def _check_python_tool(key: str) -> dict[str, Any]:
    if key == "cellphonedb":
        # CellPhoneDB may use a separate Python (e.g. one inside WSL), so a
        # locked-down Windows that blocks SciPy's DLLs can still run it.
        tp = tools_python()
        if tp:
            if not pylauncher.exists(tp):
                return {"ok": False, "message": f"The CellPhoneDB Python was not found: {tp}"}
            ok, last = _run(pylauncher.command(tp, ["-c", "import cellphonedb, importlib.metadata as m; print(m.version('cellphonedb'))"]))
            if ok:
                return {"ok": True, "message": f"CellPhoneDB {last} ({'WSL' if pylauncher.is_wsl(tp) else 'native'})"}
            return {"ok": False, "message": "cellphonedb could not be imported by the configured CellPhoneDB Python: " + last}
    py = analysis_python()
    if not py or not Path(py).exists():
        return {"ok": False, "message": "Set the analysis environment first."}
    if key == "cellphonedb":
        ok, last = _run([py, "-c", "import cellphonedb, importlib.metadata as m; print(m.version('cellphonedb'))"])
        return {"ok": ok, "message": f"CellPhoneDB {last}" if ok else
                "cellphonedb is not installed in the analysis environment (on a locked-down Windows, set a WSL Python for CellPhoneDB)."}
    natmi = config.PROJECT_ROOT / "external" / "NATMI" / "ExtractEdges.py"
    if not natmi.exists():
        return {"ok": False, "message": f"NATMI not found at {natmi}"}
    ok, last = _run([py, "-c", "import pandas, numpy, xlsxwriter, openpyxl; print('ok')"])
    return {"ok": ok, "message": "NATMI (bundled with the platform)" if ok else
            "NATMI needs pandas, numpy, xlsxwriter and openpyxl in the analysis environment: " + last}


def _resolved_rscript() -> str:
    """The Rscript to use, resolved (and cached) once: the configured/PATH R if it
    can load CellChat, otherwise a WSL R that can. Also caches the CellChat result."""
    if "rscript" in _CHECK_CACHE:
        return _CHECK_CACHE["rscript"]
    r = rscript()
    res = _check_cellchat(r)
    resolved = r or ""
    if not res["ok"] and not load().get("rscript") and os.name == "nt":
        for cand in _wsl_rscripts():
            cand_res = _check_cellchat(cand)
            if cand_res["ok"]:
                save({"rscript": cand})      # persist the discovered R (this clears the cache)
                resolved, res = cand, cand_res
                break
    _CHECK_CACHE["rscript"] = resolved
    _CHECK_CACHE["cellchat"] = res
    return resolved


def _check_r_tool(key: str) -> dict[str, Any]:
    if key == "cellchat":
        _resolved_rscript()
        return _CHECK_CACHE.get("cellchat", {"ok": False, "message": "R not checked."})
    rr = _resolved_rscript()
    if not rr or not rlauncher.exists(rr):
        return {"ok": False, "message": "Rscript not found."}
    if key == "singlecellsignalr":
        ok, last = _run(rlauncher.command(rr, ["-e", "suppressPackageStartupMessages(library(SingleCellSignalR)); "
                                               "v <- as.character(packageVersion('SingleCellSignalR')); "
                                               "if (!exists('cell_signaling', envir = asNamespace('SingleCellSignalR'))) "
                                               "stop(paste('version', v, 'is not supported; version 1 is needed')); cat(v)"]), timeout=300)
        return {"ok": ok, "message": f"SingleCellSignalR {last}" if ok else
                "SingleCellSignalR is not installed in this R (or is version 2): " + last}
    ok, last = _run(rlauncher.command(rr, ["-e", "suppressPackageStartupMessages(library(icellnet)); "
                                           "cat(as.character(packageVersion('icellnet')))"]), timeout=300)
    return {"ok": ok, "message": f"ICELLNET {last}" if ok else "ICELLNET is not installed in this R: " + last}


def check_one(tool: str, force: bool = False) -> dict[str, Any]:
    """Check (and load) a single tool, caching the result for the session."""
    key = _TOOL_KEYS.get(tool, tool.lower())
    if not force and key in _CHECK_CACHE and key != "rscript":
        return _CHECK_CACHE[key]
    res = _check_python_tool(key) if key in ("cellphonedb", "natmi") else _check_r_tool(key)
    _CHECK_CACHE[key] = res
    return res


def method_status() -> dict[str, Any]:
    """Fast: configured env + whatever tool results are already cached (no probing)."""
    tools = {}
    for name, key in _TOOL_KEYS.items():
        tools[name] = _CHECK_CACHE.get(key) if key in _CHECK_CACHE else {"ok": False, "checked": False, "message": "Select to load."}
        if tools[name].get("checked") is None:
            tools[name] = {**tools[name], "checked": True}
    flat = {key: tools[name] for name, key in _TOOL_KEYS.items()}
    return {"python": analysis_python() or "", "rscript": _CHECK_CACHE.get("rscript") or rscript() or "",
            "tools_python": tools_python() or "", "tools": tools,
            "any_ok": any(t.get("ok") for t in tools.values()), **flat}


def check_methods(force: bool = False) -> dict[str, Any]:
    """Can each CCI method run on this computer? Checks every tool (used by Setup)."""
    if force:
        clear_check_cache()
    tools = {name: check_one(name) for name in _TOOL_KEYS}
    return {"python": analysis_python() or "", "rscript": _CHECK_CACHE.get("rscript") or rscript() or "",
            "tools_python": tools_python() or "",
            "cellphonedb": tools["CellPhoneDB"], "natmi": tools["NATMI"], "cellchat": tools["CellChat"],
            "singlecellsignalr": tools["SingleCellSignalR"], "icellnet": tools["ICELLNET"],
            "tools": tools, "all_ok": all(v["ok"] for v in tools.values()), "any_ok": any(v["ok"] for v in tools.values())}


def ephemeral_mode() -> bool:
    """Demo mode: wipe all datasets and assessments on every launch.

    Defaults to True (the user asked for a clean slate each time). Turn it off
    in Setup once you want datasets and assessments to persist.
    """
    val = load().get("ephemeral_mode", True)
    return bool(val)
