"""
Expression evidence in the user's own dataset (section D).

Reads the matrix written by tools/prepare_dataset.py with numpy only (no
anndata/scipy needed in the web environment). Counts are stored gene by gene,
so looking up any gene is a slice.

Reported per gene and population:
  detected   share of cells with a count > 0
  mean       mean of log1p(counts per 10,000) over all cells of the population
  rank       position of the population among all populations by detection
For a multi-subunit receptor, the share of receiver cells in which ALL required
subunits are detected in the same cell is reported as well.
"""

from __future__ import annotations

import threading
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import datasets, index_store

_lock = threading.Lock()


class ExpressionStore:
    def __init__(self, folder: Path, raw_counts: bool = True):
        self.folder = folder
        # Raw counts are normalised here (log1p of counts per 10,000). Values that
        # are already normalised are used as given, never normalised twice.
        self.raw_counts = raw_counts
        self.genes = (folder / "genes.txt").read_text(encoding="utf-8").split()
        self.gene_index: dict[str, int] = {}
        for i, g in enumerate(self.genes):
            self.gene_index.setdefault(g.upper(), i)
        self.cells = pd.read_csv(folder / "cells.tsv", sep="\t", dtype={"sample": str, "group": str}, keep_default_na=False)
        self.indptr = np.load(folder / "csc_indptr.npy", mmap_mode="r")
        self.indices = np.load(folder / "csc_indices.npy", mmap_mode="r")
        self.data = np.load(folder / "csc_data.npy", mmap_mode="r")
        self.library = np.array(self.cells["library_size"], dtype=float)  # writable copy
        self.library[self.library <= 0] = 1.0
        self.cell_type = self.cells["cell_type"].astype(str).to_numpy()
        self.types = sorted(set(self.cell_type))
        self.type_rows = {t: np.flatnonzero(self.cell_type == t) for t in self.types}

    def column(self, gene: str) -> tuple[np.ndarray, np.ndarray] | None:
        i = self.gene_index.get(gene.upper())
        if i is None:
            return None
        a, b = int(self.indptr[i]), int(self.indptr[i + 1])
        rows = np.asarray(self.indices[a:b])
        vals = np.asarray(self.data[a:b], dtype=float)
        keep = vals > 0
        return rows[keep], vals[keep]

    def gene_stats(self, gene: str) -> dict[str, Any] | None:
        col = self.column(gene)
        if col is None:
            return None
        rows, vals = col
        lognorm = np.log1p(vals / self.library[rows] * 1e4) if self.raw_counts else vals
        detected = np.zeros(len(self.cell_type), dtype=bool)
        detected[rows] = True
        expr = np.zeros(len(self.cell_type))
        expr[rows] = lognorm
        per_type = {}
        for t, r in self.type_rows.items():
            per_type[t] = {
                "n_cells": int(len(r)),
                "detected": float(detected[r].mean()) if len(r) else 0.0,
                "mean": float(expr[r].mean()) if len(r) else 0.0,
            }
        return {"detected_mask": detected, "per_type": per_type}


@lru_cache(maxsize=4)
def _store(folder: str, _stamp: float, raw_counts: bool) -> ExpressionStore:
    return ExpressionStore(Path(folder), raw_counts=raw_counts)


def _is_raw(ds) -> bool:
    counts = (ds.meta.get("validation") or {}).get("counts") or {}
    return bool(counts.get("looks_like_raw_counts", True))


def store_for_active() -> ExpressionStore | None:
    ds = datasets.active()
    folder = ds.expression_dir
    if folder is None:
        return None
    with _lock:
        return _store(str(folder), (folder / "cells.tsv").stat().st_mtime, _is_raw(ds))


def _side(store: ExpressionStore, entity: dict[str, Any], population: str) -> dict[str, Any]:
    genes = entity.get("subunits") or []
    out = {"label": entity.get("label"), "kind": entity.get("kind"), "population": population,
           "n_cells": int(len(store.type_rows.get(population, []))), "subunits": [], "joint": None,
           "note": None}
    if population not in store.type_rows:
        out["note"] = f"The population “{population}” is not in this dataset’s cell annotation."
        return out
    if not genes:
        out["note"] = ("This partner is not resolved to verified genes (small molecule or "
                       "unresolved CellPhoneDB entity), so its expression is not assessed.")
        return out

    masks = []
    for g in genes:
        st = store.gene_stats(g)
        if st is None:
            out["subunits"].append({"gene": g, "found": False})
            continue
        pt = st["per_type"]
        order = sorted(pt, key=lambda t: pt[t]["detected"], reverse=True)
        top = [{"population": t, "detected": round(pt[t]["detected"], 4)} for t in order[:3]]
        out["subunits"].append({
            "gene": g, "found": True,
            "detected": round(pt[population]["detected"], 4),
            "mean": round(pt[population]["mean"], 3),
            "rank": order.index(population) + 1,
            "n_populations": len(order),
            "top": top,
        })
        masks.append(st["detected_mask"])

    found = [s for s in out["subunits"] if s.get("found")]
    if len(genes) > 1 and len(masks) == len(genes):
        r = store.type_rows[population]
        joint = np.logical_and.reduce([m[r] for m in masks]) if len(r) else np.array([])
        limiting = min(found, key=lambda s: s["detected"])
        out["joint"] = {
            "detected": round(float(joint.mean()) if joint.size else 0.0, 4),
            "limiting_subunit": limiting["gene"],
        }
    missing = [s["gene"] for s in out["subunits"] if not s.get("found")]
    if missing:
        out["note"] = "Not measured in this dataset: " + ", ".join(missing) + "."
    return out


def dataset_evidence(result_id: int) -> dict[str, Any] | None:
    rec = index_store.get_result(result_id)
    if not rec:
        return None
    ds = datasets.active()
    store = store_for_active()
    if store is None:
        return {
            "available": False,
            "message": (
                f"The expression matrix of “{ds.label}” has not been prepared yet. "
                "Run web/tools/prepare_dataset.py once for this dataset (see web/README.md)."
            ),
        }

    lig = index_store.entity_display(rec["ligand_entity_key"], rec.get("m_ligand"))
    rcp = index_store.entity_display(rec["receptor_entity_key"], rec.get("m_receptor"))
    groups = sorted({g for g in store.cells["group"] if g})
    samples = sorted({s for s in store.cells["sample"] if s})

    if len(groups) >= 2 and len(samples) >= 2:
        pseudobulk = {"available": False, "reason": (
            "Experimental groups are defined. Condition comparison in the receiver "
            "(pseudobulk) is the next step to be added to the web app.")}
    else:
        pseudobulk = {"available": False, "reason": (
            "This dataset has no experimental groups with several samples each, so "
            "expression cannot be compared between conditions. Define --sample-column "
            "and --group-column when preparing a dataset to enable this.")}

    if store.raw_counts:
        mean_def = "Mean of log1p(counts per 10,000) over all cells of the population, including cells with zero."
        data_note = None
    else:
        mean_def = ("Mean of the values stored in the dataset over all cells of the population, "
                    "including cells with zero. The dataset contains normalised values, not raw "
                    "counts, so they are used as given.")
        data_note = ("This dataset contains normalised values, not raw counts. Detection is "
                     "reliable; mean expression is on the dataset’s own scale.")

    return {
        "available": True,
        "dataset": ds.label,
        "data_note": data_note,
        "receptor": _side(store, rcp, rec["receiver"]),
        "ligand": _side(store, lig, rec["sender"]),
        "pseudobulk": pseudobulk,
        "definitions": {
            "detected": "Share of cells in the population with at least one count.",
            "mean": mean_def,
            "rank": "Position of this population among all populations, ordered by detection.",
            "joint": "Share of receiver cells in which every required subunit is detected in the same cell.",
        },
        "caution": (
            "RNA detected in the receiver is compatible with the interaction but does not "
            "show that the receptor protein is on the cell surface, that it binds this "
            "ligand, or that downstream signalling is active."
        ),
    }
