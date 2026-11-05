from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


# --------------------------------------------------
# Project paths
# --------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Real dataset id: the scripts run inside datasets/<id>/runs/<slug>/, so the
# run folder's grandparent is the dataset id. Outside a run folder (the
# original development layout) the example id is kept.
def _dataset_id_from_run_folder(root):
    try:
        if root.parent.name == "runs" and root.parent.parent.name:
            return root.parent.parent.name
    except Exception:
        pass
    return "endometrium_example"


DATASET_ID = _dataset_id_from_run_folder(PROJECT_ROOT)

H5AD_FILE = (
    PROJECT_ROOT
    / "data"
    / "cellphonedb"
    / "endometrium_example_counts.h5ad"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "cellchat"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# --------------------------------------------------
# Output files
# --------------------------------------------------
MATRIX_FILE = (
    OUTPUT_DIR
    / "endometrium_normalized_expression.mtx"
)

GENES_FILE = (
    OUTPUT_DIR
    / "endometrium_genes.tsv"
)

CELLS_FILE = (
    OUTPUT_DIR
    / "endometrium_cells.tsv"
)

METADATA_FILE = (
    OUTPUT_DIR
    / "endometrium_metadata.tsv"
)

PROVENANCE_FILE = (
    OUTPUT_DIR
    / "input_provenance.csv"
)


# --------------------------------------------------
# Load AnnData
# --------------------------------------------------
print("Loading AnnData:")
print(H5AD_FILE)

adata = ad.read_h5ad(
    H5AD_FILE
)


print(
    f"\nShape: "
    f"{adata.n_obs:,} cells × "
    f"{adata.n_vars:,} genes"
)


# --------------------------------------------------
# Validate cell-type annotation
# --------------------------------------------------
CELL_TYPE_COLUMN = "cell_type"


if (
    CELL_TYPE_COLUMN
    not in adata.obs.columns
):
    raise ValueError(
        f"Required annotation column "
        f"'{CELL_TYPE_COLUMN}' was not found."
    )


if (
    adata.obs[
        CELL_TYPE_COLUMN
    ].isna().any()
):
    raise ValueError(
        "Missing cell-type annotations were found."
    )


# --------------------------------------------------
# Cell identifiers
# --------------------------------------------------
cell_ids = (
    adata.obs_names
    .astype(str)
    .tolist()
)


if len(cell_ids) != len(set(cell_ids)):
    raise ValueError(
        "Duplicate cell identifiers were found."
    )


# --------------------------------------------------
# Gene identifiers
# --------------------------------------------------
gene_names = (
    adata.var_names
    .astype(str)
    .tolist()
)


duplicate_gene_count = int(
    pd.Series(
        gene_names
    )
    .duplicated()
    .sum()
)


if duplicate_gene_count > 0:
    raise ValueError(
        f"Found {duplicate_gene_count} duplicate "
        "gene names. Resolve duplicated biological "
        "gene symbols before running CellChat."
    )


# --------------------------------------------------
# Expression matrix
#
# AnnData:
#     cells × genes
#
# CellChat:
#     genes × cells
#
# IMPORTANT:
# adata.X is already normalized / transformed.
# Do NOT normalize it again.
# --------------------------------------------------
X = adata.X


if sparse.issparse(X):
    X = X.tocsr()

else:
    X = sparse.csr_matrix(
        np.asarray(X)
    )


if X.nnz == 0:
    raise ValueError(
        "Expression matrix contains no "
        "non-zero values."
    )


minimum_nonzero = float(
    X.data.min()
)

maximum_value = float(
    X.data.max()
)


if minimum_nonzero < 0:
    raise ValueError(
        "Negative values were found in adata.X."
    )


# --------------------------------------------------
# Diagnostic
# --------------------------------------------------
sample_values = X.data[
    : min(
        X.nnz,
        500000,
    )
]


integer_like = bool(
    np.allclose(
        sample_values,
        np.round(
            sample_values
        ),
        atol=1e-6,
    )
)


print(
    "\n=== INPUT DIAGNOSTICS ==="
)

print(
    "Minimum non-zero expression:",
    minimum_nonzero,
)

print(
    "Maximum expression:",
    maximum_value,
)

print(
    "Integer-like:",
    integer_like,
)


if integer_like:
    raise ValueError(
        "adata.X appears integer-like. "
        "This script is configured for "
        "normalized expression input."
    )


print(
    "Interpretation: using adata.X directly "
    "as normalized CellChat expression."
)


# --------------------------------------------------
# Transpose to genes × cells
# --------------------------------------------------
expression = (
    X.T
    .tocsc()
)


# --------------------------------------------------
# Cell metadata
# --------------------------------------------------
metadata = pd.DataFrame(
    {
        "cell_id":
            cell_ids,

        "labels":
            (
                adata.obs[
                    CELL_TYPE_COLUMN
                ]
                .astype(str)
                .values
            ),
    }
)


# --------------------------------------------------
# Alignment checks
# --------------------------------------------------
if (
    expression.shape[0]
    != len(gene_names)
):
    raise ValueError(
        "Gene dimension mismatch."
    )


if (
    expression.shape[1]
    != len(cell_ids)
):
    raise ValueError(
        "Cell dimension mismatch."
    )


if len(metadata) != len(cell_ids):
    raise ValueError(
        "Metadata dimension mismatch."
    )


# --------------------------------------------------
# Write expression matrix
# --------------------------------------------------
print(
    "\nWriting Matrix Market expression matrix..."
)


mmwrite(
    MATRIX_FILE,
    expression,
)


# --------------------------------------------------
# Write genes
# --------------------------------------------------
pd.DataFrame(
    {
        "gene":
            gene_names,
    }
).to_csv(
    GENES_FILE,
    sep="\t",
    index=False,
)


# --------------------------------------------------
# Write cells
# --------------------------------------------------
pd.DataFrame(
    {
        "cell_id":
            cell_ids,
    }
).to_csv(
    CELLS_FILE,
    sep="\t",
    index=False,
)


# --------------------------------------------------
# Write metadata
# --------------------------------------------------
metadata.to_csv(
    METADATA_FILE,
    sep="\t",
    index=False,
)


# --------------------------------------------------
# Provenance
# --------------------------------------------------
provenance = pd.DataFrame(
    [
        {
            "dataset_id":
                DATASET_ID,

            "source_file":
                str(H5AD_FILE),

            "expression_source":
                "adata.X",

            "expression_state":
                "pre-normalized / non-integer expression",

            "additional_normalization":
                False,

            "log_transform_in_cellchat_script":
                False,

            "cells":
                adata.n_obs,

            "genes":
                adata.n_vars,

            "cell_type_column":
                CELL_TYPE_COLUMN,

            "cell_types":
                metadata[
                    "labels"
                ].nunique(),

            "nonzero_values":
                expression.nnz,

            "minimum_nonzero":
                minimum_nonzero,

            "maximum_value":
                maximum_value,

            "integer_like":
                integer_like,
        }
    ]
)


provenance.to_csv(
    PROVENANCE_FILE,
    index=False,
)


# --------------------------------------------------
# Summary
# --------------------------------------------------
print(
    "\n=== CELLCHAT INPUT SUMMARY ==="
)


print(
    "Genes:",
    f"{expression.shape[0]:,}",
)

print(
    "Cells:",
    f"{expression.shape[1]:,}",
)

print(
    "Cell types:",
    metadata[
        "labels"
    ].nunique(),
)

print(
    "Non-zero expression values:",
    f"{expression.nnz:,}",
)


print(
    "\nCell-type counts:"
)


print(
    metadata[
        "labels"
    ]
    .value_counts()
    .to_string()
)


print(
    "\n=== SAVED FILES ==="
)


for path in [
    MATRIX_FILE,
    GENES_FILE,
    CELLS_FILE,
    METADATA_FILE,
    PROVENANCE_FILE,
]:
    print(path)