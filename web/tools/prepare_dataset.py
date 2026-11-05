"""
Prepare a single-cell dataset for the CCC Platform.

Run this ONCE per dataset, in the analysis environment (the one that has
anndata/scanpy installed). It checks the dataset, and writes the expression
matrix in a simple format that the web app reads with numpy only:

    datasets/<id>/dataset.json          description + validation report
    datasets/<id>/expression/genes.txt  gene symbols
    datasets/<id>/expression/cells.tsv  cell, cell_type, sample, group, library_size
    datasets/<id>/expression/csc_*.npy  counts, stored gene by gene (CSC)

Example (endometrium use case, from the project folder):

    python web/tools/prepare_dataset.py ^
        --id endometrium_example ^
        --h5ad data/cellphonedb/endometrium_example_counts.h5ad ^
        --metadata data/cellphonedb/endometrium_example_meta.tsv ^
        --cell-type-column cell_type

Options for other datasets: --label, --organism, --tissue, --layer (use a
counts layer instead of X), --gene-symbol-column (var column with symbols),
--sample-column and --group-column (needed for condition comparisons).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULT_MARK = "@@RESULT@@ "
MIN_CELLS_WARNING = 20


def _to_csc(matrix):
    from scipy import sparse
    if sparse.issparse(matrix):
        return sparse.csc_matrix(matrix)
    return sparse.csc_matrix(np.asarray(matrix))


def _check_counts(csc) -> dict[str, Any]:
    data = np.asarray(csc.data)
    sample = data[: min(len(data), 200_000)]
    if sample.size == 0:
        return {"looks_like_raw_counts": False, "note": "The matrix has no non-zero values."}
    integer_like = bool(np.all(np.abs(sample - np.round(sample)) < 1e-6))
    negative = bool((sample < 0).any())
    raw = integer_like and not negative and float(sample.max()) >= 1
    note = None if raw else (
        "Values are not whole numbers, so they look normalised or transformed. "
        "Detection (share of cells with expression > 0) is still valid; mean "
        "expression is computed on the values as given; condition comparisons "
        "(pseudobulk) need raw counts. Use --layer to point to a counts layer if one exists."
    )
    return {"looks_like_raw_counts": raw, "min": float(sample.min()), "max": float(sample.max()), "note": note}


def _check_genes(genes: list[str]) -> dict[str, Any]:
    ensembl = sum(bool(re.match(r"^ENS[A-Z]*G\d+", g)) for g in genes)
    share = ensembl / max(1, len(genes))
    return {
        "n_genes": len(genes),
        "ensembl_share": round(share, 3),
        "note": ("Most gene names look like Ensembl IDs. Ligand and receptor genes are "
                 "matched by symbol; use --gene-symbol-column.") if share > 0.5 else None,
    }


_SUGGEST = {
    "cell_type": ("cell_type", "celltype", "cell_types", "annotation", "cell_annotation", "cluster_label", "cell_label", "labels"),
    "sample": ("sample", "sample_id", "donor", "donor_id", "patient", "patient_id", "individual", "orig.ident", "subject", "batch"),
    "group": ("condition", "group", "disease", "status", "treatment", "state", "diagnosis", "phenotype"),
}


def inspect(adata, metadata: pd.DataFrame | None = None) -> dict[str, Any]:
    """Describe a dataset so the user can choose its columns."""
    obs = adata.obs.copy()
    obs.index = obs.index.astype(str)
    matched = None
    if metadata is not None:
        meta = metadata.copy()
        meta.index = meta.index.astype(str)
        matched = int(obs.index.isin(meta.index).sum())
        common = obs.index[obs.index.isin(meta.index)]
        for col in meta.columns:
            obs.loc[common, col] = meta.loc[common, col].values
    columns = []
    for col in obs.columns:
        values = obs[col].dropna().astype(str)
        n_unique = int(values.nunique())
        suggested = None
        low = str(col).lower()
        for role, names in _SUGGEST.items():
            if low in names or any(low.startswith(n) for n in names):
                suggested = role
                break
        if n_unique <= 1 or n_unique > max(500, len(obs) // 2):
            suggested = None     # constant, or looks like a per-cell value
        columns.append({
            "name": str(col),
            "n_unique": n_unique,
            "examples": [str(v) for v in values.value_counts().head(6).index],
            "suggested": suggested,
        })
    try:
        sample = adata.X[: min(adata.n_obs, 500)]
        counts = _check_counts(_to_csc(sample))
    except Exception as exc:          # backed files with unusual storage
        counts = {"looks_like_raw_counts": None, "note": f"Could not read values: {exc}"}
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "columns": columns,
        "layers": [str(k) for k in getattr(adata, "layers", {}).keys()],
        "has_raw": getattr(adata, "raw", None) is not None,
        "counts": counts,
        "genes": _check_genes([str(g) for g in adata.var_names[:5000]]),
        "var_columns": [str(c) for c in adata.var.columns],
        "metadata_cells_matched": matched,
    }


def prepare(
    adata,
    *,
    ds_id: str,
    cell_type_column: str,
    label: str | None = None,
    organism: str | None = None,
    tissue: str | None = None,
    metadata: pd.DataFrame | None = None,
    layer: str | None = None,
    use_raw: bool = False,
    gene_symbol_column: str | None = None,
    sample_column: str | None = None,
    group_column: str | None = None,
    source: dict[str, str] | None = None,
    datasets_dir: Path | None = None,
) -> dict[str, Any]:
    datasets_dir = datasets_dir or (PROJECT_ROOT / "datasets")
    folder = datasets_dir / ds_id
    out = folder / "expression"
    out.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    obs = adata.obs.copy()
    obs.index = obs.index.astype(str)
    if metadata is not None:
        meta = metadata.copy()
        meta.index = meta.index.astype(str)
        missing = int((~obs.index.isin(meta.index)).sum())
        if missing:
            warnings.append(f"{missing} cells are not in the metadata file and are left out.")
        keep = obs.index.isin(meta.index)
        adata = adata[keep]
        obs = obs.loc[keep]
        for col in meta.columns:
            obs[col] = meta.loc[obs.index, col].values

    for col, what in ((cell_type_column, "cell type"), (sample_column, "sample"), (group_column, "group")):
        if col and col not in obs.columns:
            raise SystemExit(f"The {what} column '{col}' was not found. Available columns: {list(obs.columns)}")

    available_layers = [str(k) for k in getattr(adata, "layers", {}).keys()]
    has_raw = getattr(adata, "raw", None) is not None
    if use_raw and not has_raw:
        raise SystemExit("--use-raw was given but the file has no .raw matrix.")
    matrix = adata.raw.X if use_raw else (adata.layers[layer] if layer else adata.X)
    csc = _to_csc(matrix)
    counts_check = _check_counts(csc)
    if counts_check.get("note"):
        other = [f"layers: {', '.join(available_layers)}" if available_layers else None,
                 "a .raw matrix (use --use-raw)" if has_raw and not use_raw else None]
        other = [o for o in other if o]
        warnings.append(counts_check["note"] + (" The file also contains " + " and ".join(other) + "."
                                                if other else " The file contains no other matrix (no layers, no .raw)."))

    var = adata.raw.var if use_raw else adata.var
    if gene_symbol_column:
        genes = [str(g) for g in var[gene_symbol_column]]
    else:
        genes = [str(g) for g in var.index]
    genes_check = _check_genes(genes)
    if genes_check.get("note"):
        warnings.append(genes_check["note"])

    cell_types = obs[cell_type_column].astype(str)
    per_type = cell_types.value_counts().sort_index()
    small = per_type[per_type < MIN_CELLS_WARNING]
    if len(small):
        warnings.append(
            "Populations with fewer than %d cells give unstable expression estimates: %s."
            % (MIN_CELLS_WARNING, ", ".join(f"{k} ({v})" for k, v in small.items()))
        )

    library = np.asarray(csc.sum(axis=1)).ravel()
    cells = pd.DataFrame({
        "cell": obs.index,
        "cell_type": cell_types.values,
        "sample": obs[sample_column].astype(str).values if sample_column else "",
        "group": obs[group_column].astype(str).values if group_column else "",
        "library_size": library,
    })
    cells.to_csv(out / "cells.tsv", sep="\t", index=False)
    (out / "genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    np.save(out / "csc_indptr.npy", csc.indptr.astype(np.int64))
    np.save(out / "csc_indices.npy", csc.indices.astype(np.int32))
    np.save(out / "csc_data.npy", csc.data.astype(np.float32))

    groups = sorted(set(cells["group"])) if group_column else []
    samples = sorted(set(cells["sample"])) if sample_column else []
    report = {
        "n_cells": int(len(cells)),
        "n_genes": len(genes),
        "n_cell_types": int(len(per_type)),
        "cells_per_type": {str(k): int(v) for k, v in per_type.items()},
        "samples": samples,
        "groups": groups,
        "counts": counts_check,
        "genes": genes_check,
        "warnings": warnings,
        "matrix_source": "raw.X" if use_raw else (f"layers['{layer}']" if layer else "X"),
        "available_layers": available_layers,
        "has_raw": has_raw,
    }

    ds_file = folder / "dataset.json"
    meta_json = json.loads(ds_file.read_text(encoding="utf-8")) if ds_file.exists() else {"id": ds_id}
    meta_json.update({k: v for k, v in {
        "label": label, "organism": organism, "tissue": tissue,
        "cell_type_column": cell_type_column, "sample_column": sample_column,
        "group_column": group_column,
    }.items() if v})
    meta_json.setdefault("label", ds_id)
    if source:
        meta_json["source"] = {**meta_json.get("source", {}), **source}
    meta_json["expression"] = "expression"
    meta_json["validation"] = report
    ds_file.write_text(json.dumps(meta_json, indent=2, ensure_ascii=False), encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--id", help="Short identifier, e.g. endometrium_example")
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--cell-type-column")
    ap.add_argument("--inspect", action="store_true", help="Only describe the file (columns, layers) as JSON")
    ap.add_argument("--json", action="store_true", help="Print a machine-readable result line")
    ap.add_argument("--metadata", help="Optional TSV/CSV: first column = cell id, other columns = annotations")
    ap.add_argument("--label"); ap.add_argument("--organism"); ap.add_argument("--tissue")
    ap.add_argument("--layer"); ap.add_argument("--gene-symbol-column")
    ap.add_argument("--use-raw", action="store_true", help="Use adata.raw (often holds the raw counts)")
    ap.add_argument("--sample-column"); ap.add_argument("--group-column")
    a = ap.parse_args(argv)

    try:
        import anndata as ad
    except ImportError:
        sys.exit("anndata is not installed in this environment. Run this script in the analysis environment.")

    h5ad = Path(a.h5ad)
    if not h5ad.exists():
        sys.exit(f"File not found: {h5ad}")
    metadata = None
    if a.metadata:
        if not Path(a.metadata).exists():
            sys.exit(f"File not found: {a.metadata}")
        sep = "\t" if a.metadata.endswith((".tsv", ".txt")) else ","
        metadata = pd.read_csv(a.metadata, sep=sep, index_col=0)

    if a.inspect:
        try:
            adata = ad.read_h5ad(h5ad, backed="r")     # fast: does not load the matrix
        except Exception:
            adata = ad.read_h5ad(h5ad)
        print(RESULT_MARK + json.dumps(inspect(adata, metadata)))
        return

    if not a.id or not a.cell_type_column:
        sys.exit("--id and --cell-type-column are required.")
    print(f"Reading {h5ad} …", flush=True)
    adata = ad.read_h5ad(h5ad)
    print("Checking and converting the matrix …", flush=True)

    def rel(p):
        try:
            return str(Path(p).resolve().relative_to(PROJECT_ROOT)).replace("\\", "/")
        except ValueError:
            return str(Path(p).resolve())

    report = prepare(
        adata, ds_id=a.id, cell_type_column=a.cell_type_column, label=a.label,
        organism=a.organism, tissue=a.tissue, metadata=metadata, layer=a.layer, use_raw=a.use_raw,
        gene_symbol_column=a.gene_symbol_column, sample_column=a.sample_column,
        group_column=a.group_column,
        source={"counts": rel(h5ad), **({"metadata": rel(a.metadata)} if a.metadata else {})},
    )
    print(f"\nDataset '{a.id}': {report['n_cells']} cells, {report['n_genes']} genes, "
          f"{report['n_cell_types']} cell types.")
    for w in report["warnings"]:
        print("  Note:", w)
    if a.json:
        print(RESULT_MARK + json.dumps(report))
    print("\nDone.")


if __name__ == "__main__":
    main()
