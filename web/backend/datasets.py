"""
Datasets.

Every dataset is a folder ``<project>/datasets/<id>/`` with a
``dataset.json`` describing it:

    {
      "id": "melanoma_tme",
      "label": "Melanoma tumour microenvironment (Tirosh 2016)",
      "organism": "Human",
      "tissue": "Melanoma tumour microenvironment",
      "cell_type_column": "cell_type",
      "results": {
        "consensus": "results/.../consensus_cellphonedb__natmi__cellchat.csv",
        "cellchat": "results/cellchat/standardized_interactions.csv",
        "cellchat_run_info": "results/cellchat/run_info.csv"
      },
      "expression": "expression"          (optional; written by prepare_dataset.py)
    }

Paths in "results" are relative to the project folder; "expression" is
relative to the dataset folder. Everything the web app writes for a dataset
(results index, experiment description) goes to ``web/workspace/<id>/``.
"""

from __future__ import annotations

import json
import shutil
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import config

DATASETS_DIR = config.PROJECT_ROOT / "datasets"
ACTIVE_FILE = config.WORKSPACE_DIR / "active_dataset.txt"
_lock = threading.Lock()


@dataclass(frozen=True)
class Dataset:
    id: str
    label: str
    folder: Path
    meta: dict[str, Any]

    def _result_path(self, key: str) -> Path | None:
        rel = (self.meta.get("results") or {}).get(key)
        return (config.PROJECT_ROOT / rel) if rel else None

    @property
    def consensus_file(self) -> Path | None:
        return self._result_path("consensus")

    @property
    def cellchat_file(self) -> Path | None:
        return self._result_path("cellchat")

    @property
    def cellchat_run_info(self) -> Path | None:
        return self._result_path("cellchat_run_info")

    @property
    def expression_dir(self) -> Path | None:
        rel = self.meta.get("expression")
        path = (self.folder / rel) if rel else None
        return path if path and (path / "cells.tsv").exists() else None

    @property
    def workspace(self) -> Path:
        path = config.WORKSPACE_DIR / self.id
        path.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def index_db(self) -> Path:
        run = self.meta.get("active_run")
        return self.workspace / (f"interactions_index__{run}.sqlite" if run else "interactions_index.sqlite")

    @property
    def runs(self) -> list[dict]:
        """Results available for this dataset, one entry per set of tools that was run."""
        out = []
        runs = self.meta.get("runs") or {}
        if not runs and (self.meta.get("results") or {}).get("consensus"):
            r = self.meta["results"]
            runs = {"earlier": {**r, "tools": r.get("tools") or []}}
        for slug, r in runs.items():
            path = config.PROJECT_ROOT / r["consensus"] if r.get("consensus") else None
            out.append({"run": slug, "tools": r.get("tools") or [], "created": r.get("created"),
                        "available": bool(path and path.exists()),
                        "active": slug == (self.meta.get("active_run") or ("earlier" if not self.meta.get("runs") else None))})
        return out

    @property
    def context_file(self) -> Path:
        return self.workspace / "project_context.json"

    def summary(self) -> dict[str, Any]:
        report = self.meta.get("validation") or {}
        return {
            "id": self.id,
            "label": self.label,
            "organism": self.meta.get("organism"),
            "tissue": self.meta.get("tissue"),
            "has_results": bool(self.consensus_file and self.consensus_file.exists()),
            "has_expression": self.expression_dir is not None,
            "n_cells": report.get("n_cells"),
            "n_cell_types": report.get("n_cell_types"),
            "groups": report.get("groups") or [],
            "has_provenance": bool(self.meta.get("provenance")),
            "runs": self.runs,
            "active_run": self.meta.get("active_run"),
            "results_created": self.meta.get("results_created"),
        }


def _load(folder: Path) -> Dataset | None:
    f = folder / "dataset.json"
    if not f.exists():
        return None
    try:
        meta = json.loads(f.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    ds_id = str(meta.get("id") or folder.name)
    return Dataset(id=ds_id, label=str(meta.get("label") or ds_id), folder=folder, meta=meta)


def all_datasets() -> list[Dataset]:
    if not DATASETS_DIR.exists():
        return []
    out = [d for d in (_load(p) for p in sorted(DATASETS_DIR.iterdir()) if p.is_dir()) if d]
    return sorted(out, key=lambda d: d.label.lower())


def get(ds_id: str) -> Dataset | None:
    return next((d for d in all_datasets() if d.id == ds_id), None)


_NONE = "__none__"


class NoDatasetSelected(RuntimeError):
    """Raised when an operation needs a dataset but the user has not loaded one."""


def active_or_none() -> Dataset | None:
    """The dataset the user has loaded, or None if they have not chosen one.

    Unlike the old behaviour, nothing is auto-loaded: the user explicitly picks
    a dataset (upload their own, or load the example) from the start page.
    """
    with _lock:
        wanted = ACTIVE_FILE.read_text(encoding="utf-8").strip() if ACTIVE_FILE.exists() else ""
    if not wanted or wanted == _NONE:
        return None
    return get(wanted)  # None if that id no longer exists


def active() -> Dataset:
    ds = active_or_none()
    if ds is None:
        raise NoDatasetSelected("No dataset is loaded.")
    return ds


def set_active(ds_id: str) -> Dataset:
    ds = get(ds_id)
    if not ds:
        raise LookupError(f"Unknown dataset: {ds_id}")
    with _lock:
        ACTIVE_FILE.write_text(ds.id, encoding="utf-8")
    return ds


def clear_active() -> None:
    """Return to the start state with no dataset loaded."""
    with _lock:
        ACTIVE_FILE.write_text(_NONE, encoding="utf-8")


#: Datasets that ship with the platform and must not be deleted from the UI.
PROTECTED_IDS = {"melanoma_tme"}


def wipe_user_datasets() -> None:
    """Delete every dataset the user added (keep the shipped example)."""
    for d in all_datasets():
        if d.id not in PROTECTED_IDS:
            try:
                delete(d.id)
            except Exception:
                pass


def delete(ds_id: str) -> None:
    """Remove a dataset the user added: its folder and its web workspace.

    The shipped example is protected. Expert assessments in the database are
    keyed by dataset id and are left untouched (harmless if orphaned).
    """
    ds_id = str(ds_id).strip()
    if ds_id in PROTECTED_IDS:
        raise PermissionError("The example dataset cannot be removed.")
    ds = get(ds_id)
    if ds is None:
        raise LookupError(f"Unknown dataset: {ds_id}")
    if active_or_none() and active().id == ds_id:
        clear_active()
    folder = DATASETS_DIR / ds_id
    if folder.exists():
        shutil.rmtree(folder, ignore_errors=True)
    ws = config.WORKSPACE_DIR / ds_id
    if ws.exists():
        shutil.rmtree(ws, ignore_errors=True)


def migrate_legacy_workspace() -> None:
    """Earlier versions kept one index and one experiment description in web/workspace/."""
    try:
        ds = active()
    except RuntimeError:
        return
    old_ctx = config.WORKSPACE_DIR / "project_context.json"
    if old_ctx.exists() and not ds.context_file.exists():
        shutil.move(str(old_ctx), str(ds.context_file))
    for old in ("interactions_index.sqlite", "interactions_index.building"):
        p = config.WORKSPACE_DIR / old
        if p.exists():
            try:
                p.unlink()
            except OSError:
                pass


def set_run(ds_id: str, run: str) -> "Dataset":
    """Show the results of another set of tools for this dataset."""
    import json as _json
    ds = get(ds_id)
    if ds is None:
        raise ValueError("Dataset not found.")
    runs = ds.meta.get("runs") or {}
    if run not in runs:
        raise ValueError("These results do not exist for this dataset.")
    meta = dict(ds.meta)
    meta["active_run"] = run
    meta["results"] = dict(runs[run])
    if runs[run].get("provenance"):
        meta["provenance"] = runs[run]["provenance"]
    (ds.folder / "dataset.json").write_text(_json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    return get(ds_id)
