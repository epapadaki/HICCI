from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad  # only read_h5ad is needed; scanpy would also load numba
from scipy import sparse


# --------------------------------------------------
# Paths
# --------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

COUNTS_FILE = (
    PROJECT_ROOT
    / "data"
    / "cellphonedb"
    / "endometrium_example_counts.h5ad"
)

META_FILE = (
    PROJECT_ROOT
    / "data"
    / "cellphonedb"
    / "endometrium_example_meta.tsv"
)

NATMI_DATA_DIR = (
    PROJECT_ROOT
    / "data"
    / "natmi"
)

EXPRESSION_FILE = (
    NATMI_DATA_DIR
    / "endometrium_expression.csv"
)

ANNOTATION_FILE = (
    NATMI_DATA_DIR
    / "endometrium_annotations.csv"
)

NATMI_DATA_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# --------------------------------------------------
# Load AnnData and metadata
# --------------------------------------------------
print("Loading AnnData...")

adata = ad.read_h5ad(
    COUNTS_FILE
)

print(
    f"AnnData shape: "
    f"{adata.n_obs:,} cells × "
    f"{adata.n_vars:,} genes"
)


print("\nLoading metadata...")

metadata = pd.read_csv(
    META_FILE,
    sep="\t",
)


# --------------------------------------------------
# Validate metadata
# --------------------------------------------------
required_columns = {
    "Cell",
    "cell_type",
}

missing_columns = (
    required_columns
    - set(metadata.columns)
)

if missing_columns:
    raise ValueError(
        "Missing required metadata columns: "
        + ", ".join(
            sorted(missing_columns)
        )
    )


metadata["Cell"] = (
    metadata["Cell"]
    .astype(str)
)

metadata["cell_type"] = (
    metadata["cell_type"]
    .astype(str)
)


adata.obs_names = (
    adata.obs_names
    .astype(str)
)


# --------------------------------------------------
# Ensure metadata matches AnnData cells
# --------------------------------------------------
metadata = (
    metadata
    .drop_duplicates(
        subset="Cell"
    )
    .set_index(
        "Cell"
    )
)


missing_cells = [
    cell
    for cell in adata.obs_names
    if cell not in metadata.index
]


if missing_cells:

    raise ValueError(
        f"{len(missing_cells)} AnnData cells "
        "do not have matching metadata."
    )


metadata = metadata.loc[
    adata.obs_names
]


print(
    f"Matched cells: "
    f"{len(metadata):,}"
)


# --------------------------------------------------
# Inspect expression matrix
# --------------------------------------------------
X = adata.X


if sparse.issparse(X):

    data_values = X.data

else:

    data_values = np.asarray(
        X
    ).ravel()


data_values = data_values[
    np.isfinite(
        data_values
    )
]


if len(data_values) == 0:

    raise ValueError(
        "Expression matrix contains no "
        "finite expression values."
    )


sample_size = min(
    len(data_values),
    100000,
)


sample = data_values[
    :sample_size
]


integer_fraction = np.mean(
    np.isclose(
        sample,
        np.round(sample),
        atol=1e-8,
    )
)


print("\nExpression diagnostics")
print("----------------------")

print(
    f"Minimum value: "
    f"{float(np.min(data_values)):.4f}"
)

print(
    f"Maximum value: "
    f"{float(np.max(data_values)):.4f}"
)

print(
    f"Mean value: "
    f"{float(np.mean(data_values)):.4f}"
)

print(
    f"Fraction approximately integer: "
    f"{integer_fraction:.3f}"
)


if "log1p" in adata.uns:

    print(
        "AnnData contains log1p metadata."
    )

else:

    print(
        "No log1p metadata found in adata.uns."
    )


# --------------------------------------------------
# NATMI requires normalized expression.
#
# For this official CellPhoneDB example dataset,
# we preserve the expression matrix used in the
# CellPhoneDB workflow rather than applying a second
# normalization step.
#
# This ensures both CCC methods start from the same
# expression representation.
# --------------------------------------------------
print(
    "\nUsing expression values as stored "
    "in the official CellPhoneDB example AnnData."
)


# --------------------------------------------------
# Convert AnnData:
#
# AnnData:
#     rows    = cells
#     columns = genes
#
# NATMI:
#     rows    = genes
#     columns = cells
# --------------------------------------------------
print(
    "\nPreparing NATMI expression matrix..."
)


if sparse.issparse(
    adata.X
):

    expression_df = pd.DataFrame.sparse.from_spmatrix(
        adata.X.T,
        index=adata.var_names.astype(str),
        columns=adata.obs_names.astype(str),
    )

else:

    expression_df = pd.DataFrame(
        np.asarray(
            adata.X
        ).T,
        index=adata.var_names.astype(str),
        columns=adata.obs_names.astype(str),
    )


expression_df.index.name = (
    "Gene"
)


# --------------------------------------------------
# Remove duplicated gene symbols
#
# NATMI works with official gene symbols.
# Duplicate rows are summed so one gene symbol
# corresponds to one row.
# --------------------------------------------------
if expression_df.index.duplicated().any():

    duplicate_count = int(
        expression_df.index
        .duplicated()
        .sum()
    )

    print(
        f"Found {duplicate_count:,} duplicated "
        "gene identifiers."
    )

    expression_df = (
        expression_df
        .groupby(
            level=0
        )
        .sum()
    )

else:

    print(
        "No duplicated gene identifiers found."
    )


# --------------------------------------------------
# Prepare NATMI annotation file
#
# Required structure:
#
# Cell,Annotation
# cell_1,Myeloid
# cell_2,Fibroblast
# --------------------------------------------------
annotation_df = pd.DataFrame(
    {
        "Cell":
            metadata.index.astype(str),

        "Annotation":
            metadata[
                "cell_type"
            ].astype(str).values,
    }
)


# --------------------------------------------------
# Final validation
# --------------------------------------------------
if list(
    expression_df.columns
) != list(
    annotation_df["Cell"]
):

    raise ValueError(
        "Cell order mismatch between "
        "expression and annotation files."
    )


print(
    f"\nNATMI expression shape: "
    f"{expression_df.shape[0]:,} genes × "
    f"{expression_df.shape[1]:,} cells"
)


print(
    f"NATMI annotations: "
    f"{len(annotation_df):,} cells"
)


print(
    f"Cell types: "
    f"{annotation_df['Annotation'].nunique()}"
)


# --------------------------------------------------
# Save NATMI inputs
# --------------------------------------------------
print(
    "\nSaving NATMI expression matrix..."
)

expression_df.to_csv(
    EXPRESSION_FILE,
)


print(
    "Saving NATMI annotation file..."
)

annotation_df.to_csv(
    ANNOTATION_FILE,
    index=False,
)


# --------------------------------------------------
# Preview
# --------------------------------------------------
print(
    "\nAnnotation preview:"
)

print(
    annotation_df
    .head(10)
    .to_string(
        index=False
    )
)


print(
    "\nFirst 10 gene identifiers:"
)

print(
    expression_df
    .index[:10]
    .tolist()
)


print(
    "\nFiles created:"
)

print(
    EXPRESSION_FILE
)

print(
    ANNOTATION_FILE
)

print(
    "\nNATMI input preparation completed successfully."
)