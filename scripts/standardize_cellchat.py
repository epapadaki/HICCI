from pathlib import Path

import numpy as np
import pandas as pd


# --------------------------------------------------
# Project paths
# --------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_FILE = (
    PROJECT_ROOT
    / "results"
    / "cellchat"
    / "cellchat_interactions_raw.csv"
)

RUN_INFO_FILE = (
    PROJECT_ROOT
    / "results"
    / "cellchat"
    / "run_info.csv"
)

GENES_FILE = (
    PROJECT_ROOT
    / "data"
    / "cellchat"
    / "endometrium_genes.tsv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "cellchat"
)

STANDARDIZED_FILE = (
    OUTPUT_DIR
    / "standardized_interactions.csv"
)

COMPARISON_FILE = (
    OUTPUT_DIR
    / "standardized_interactions_comparison_ready.csv"
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "standardization_summary.csv"
)


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
METHOD_NAME = "CellChat"


# --------------------------------------------------
# Helpers
# --------------------------------------------------
def clean_text(series):
    return (
        series
        .astype("string")
        .str.strip()
    )


def canonical_gene(series):
    return (
        clean_text(series)
        .str.upper()
    )


def create_interaction_key(df):
    return (
        df["sender_cell_type"].astype(str)
        + "||"
        + df["receiver_cell_type"].astype(str)
        + "||"
        + df["ligand"].astype(str)
        + "||"
        + df["receptor"].astype(str)
    )


def add_percentiles(df):
    """
    Higher CellChat communication probability is better.

    Percentiles are method-relative ranks.
    They are not probabilities of biological truth.
    """

    df = df.copy()

    df["method_global_percentile"] = (
        df["score"]
        .rank(
            method="average",
            ascending=True,
            pct=True,
        )
        * 100.0
    )

    df["method_cell_pair_percentile"] = (
        df.groupby(
            [
                "sender_cell_type",
                "receiver_cell_type",
            ]
        )["score"]
        .rank(
            method="average",
            ascending=True,
            pct=True,
        )
        * 100.0
    )

    return df


# --------------------------------------------------
# Validate inputs
# --------------------------------------------------
for path in [
    RAW_FILE,
    GENES_FILE,
]:
    if not path.exists():
        raise FileNotFoundError(
            f"Required file does not exist: {path}"
        )


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# --------------------------------------------------
# Load CellChat raw interactions
# --------------------------------------------------
print("Loading CellChat raw output:")
print(RAW_FILE)

raw = pd.read_csv(
    RAW_FILE
)


print(
    f"\nRaw shape: "
    f"{raw.shape[0]:,} rows × "
    f"{raw.shape[1]:,} columns"
)


required_columns = [
    "source",
    "target",
    "ligand",
    "receptor",
    "prob",
    "pval",
]


missing_columns = [
    column
    for column in required_columns
    if column not in raw.columns
]


if missing_columns:
    raise ValueError(
        "Missing required CellChat columns: "
        + ", ".join(missing_columns)
    )


# --------------------------------------------------
# Load actual dataset gene symbols
#
# We use these to distinguish individual genes
# from CellChat complex/composite molecular names.
# --------------------------------------------------
genes = pd.read_csv(
    GENES_FILE,
    sep="\t",
)


if "gene" not in genes.columns:
    raise ValueError(
        "Expected a 'gene' column in "
        "endometrium_genes.tsv."
    )


dataset_gene_set = set(
    genes["gene"]
    .astype(str)
    .str.strip()
    .str.upper()
)


print(
    "Dataset gene symbols:",
    f"{len(dataset_gene_set):,}",
)


# --------------------------------------------------
# Optional CellChat version provenance
# --------------------------------------------------
cellchat_version = "unknown"


if RUN_INFO_FILE.exists():

    run_info = pd.read_csv(
        RUN_INFO_FILE
    )

    if {
        "field",
        "value",
    }.issubset(
        run_info.columns
    ):

        version_rows = run_info[
            run_info["field"]
            == "CellChat_version"
        ]

        if not version_rows.empty:
            cellchat_version = str(
                version_rows.iloc[0][
                    "value"
                ]
            )


print(
    "CellChat version:",
    cellchat_version,
)


# --------------------------------------------------
# Standardize basic fields
# --------------------------------------------------
df = raw.copy()


df["sender_cell_type"] = clean_text(
    df["source"]
)

df["receiver_cell_type"] = clean_text(
    df["target"]
)


# Preserve native CellChat molecular labels
df["ligand_raw"] = clean_text(
    df["ligand"]
)

df["receptor_raw"] = clean_text(
    df["receptor"]
)


# Canonical comparison representation
df["ligand"] = canonical_gene(
    df["ligand_raw"]
)

df["receptor"] = canonical_gene(
    df["receptor_raw"]
)


df["prob"] = pd.to_numeric(
    df["prob"],
    errors="coerce",
)

df["pval"] = pd.to_numeric(
    df["pval"],
    errors="coerce",
)


# --------------------------------------------------
# Remove structurally invalid rows
# --------------------------------------------------
valid_mask = (
    df["sender_cell_type"].notna()
    & df["receiver_cell_type"].notna()
    & df["ligand"].notna()
    & df["receptor"].notna()
    & df["prob"].notna()
    & df["pval"].notna()
)


df = (
    df.loc[
        valid_mask
    ]
    .copy()
)


# --------------------------------------------------
# Identify gene-gene versus complex/composite events
#
# A molecule is treated as an individual gene only
# if its CellChat label maps exactly to a gene symbol
# present in the expression dataset.
#
# We deliberately do NOT explode CellChat complexes
# into artificial single-subunit interactions.
# --------------------------------------------------
df["ligand_is_dataset_gene"] = (
    df["ligand"].isin(
        dataset_gene_set
    )
)


df["receptor_is_dataset_gene"] = (
    df["receptor"].isin(
        dataset_gene_set
    )
)


df["interaction_type"] = np.where(
    (
        df["ligand_is_dataset_gene"]
        & df["receptor_is_dataset_gene"]
    ),
    "ligand-receptor",
    "complex-or-composite",
)


# --------------------------------------------------
# Method metadata
# --------------------------------------------------
df["dataset_id"] = DATASET_ID
df["method"] = METHOD_NAME
df["method_version"] = cellchat_version

df["resource"] = (
    "CellChatDB v2 human"
)

df["score"] = df["prob"]

df["score_name"] = (
    "CellChat communication probability"
)

df["score_higher_is_better"] = True


# --------------------------------------------------
# Human-readable molecular pair
# --------------------------------------------------
df["interacting_pair"] = (
    df["ligand"].astype(str)
    + " → "
    + df["receptor"].astype(str)
)


# --------------------------------------------------
# Canonical biological event
# --------------------------------------------------
df["interaction_key"] = (
    create_interaction_key(
        df
    )
)


# --------------------------------------------------
# Check duplicate canonical events
# --------------------------------------------------
duplicate_count = int(
    df.duplicated(
        subset=[
            "interaction_key",
        ]
    ).sum()
)


print(
    "\nDuplicate canonical keys "
    "before defensive collapse:",
    f"{duplicate_count:,}",
)


# --------------------------------------------------
# Defensive collapse
#
# CellChat output is expected to have one row per
# sender/receiver/ligand/receptor event.
#
# If duplicates appear, retain:
#   1. highest communication probability
#   2. lowest p-value
# --------------------------------------------------
df = (
    df.sort_values(
        [
            "score",
            "pval",
        ],
        ascending=[
            False,
            True,
        ],
    )
    .drop_duplicates(
        subset=[
            "interaction_key",
        ],
        keep="first",
    )
    .reset_index(
        drop=True
    )
)


# --------------------------------------------------
# Full standardized percentiles
#
# These rankings are within all significant CellChat
# events exported by subsetCommunication.
# --------------------------------------------------
df = add_percentiles(
    df
)


# --------------------------------------------------
# Reorder columns
# --------------------------------------------------
preferred_columns = [
    "dataset_id",
    "method",
    "method_version",
    "resource",

    "interaction_key",
    "interaction_type",

    "sender_cell_type",
    "receiver_cell_type",

    "ligand",
    "receptor",

    "ligand_raw",
    "receptor_raw",

    "interacting_pair",

    "score",
    "score_name",
    "score_higher_is_better",

    "method_global_percentile",
    "method_cell_pair_percentile",

    "pval",

    "ligand_is_dataset_gene",
    "receptor_is_dataset_gene",

    "interaction_name",
    "interaction_name_2",
    "pathway_name",
    "annotation",
    "evidence",
]


preferred_columns = [
    column
    for column in preferred_columns
    if column in df.columns
]


remaining_columns = [
    column
    for column in df.columns
    if column not in preferred_columns
]


df = df[
    preferred_columns
    + remaining_columns
]


# --------------------------------------------------
# Save full standardized CellChat output
# --------------------------------------------------
df.to_csv(
    STANDARDIZED_FILE,
    index=False,
)


# --------------------------------------------------
# Comparison-ready subset
#
# Requirements:
#   - significant CellChat inference
#   - positive communication score
#   - individual ligand gene
#   - individual receptor gene
#
# Complex/composite CellChat entities remain in the
# full standardized file but are not forced into the
# exact gene-level consensus universe.
# --------------------------------------------------
comparison = df[
    (
        df["interaction_type"]
        == "ligand-receptor"
    )
    & (
        df["score"] > 0
    )
    & (
        df["pval"] < 0.05
    )
].copy()


# --------------------------------------------------
# Recalculate percentiles within the exact universe
# used for cross-method comparison.
#
# This mirrors the NATMI comparison-ready strategy:
# method-relative ranks are defined inside the actual
# set entering consensus.
# --------------------------------------------------
comparison = add_percentiles(
    comparison
)


comparison = (
    comparison
    .sort_values(
        "score",
        ascending=False,
    )
    .reset_index(
        drop=True
    )
)


comparison.to_csv(
    COMPARISON_FILE,
    index=False,
)


# --------------------------------------------------
# Summary
# --------------------------------------------------
full_gene_gene = int(
    (
        df["interaction_type"]
        == "ligand-receptor"
    ).sum()
)


full_complex = int(
    (
        df["interaction_type"]
        == "complex-or-composite"
    ).sum()
)


unique_lr_full = int(
    df[
        [
            "ligand",
            "receptor",
        ]
    ]
    .drop_duplicates()
    .shape[0]
)


unique_lr_comparison = int(
    comparison[
        [
            "ligand",
            "receptor",
        ]
    ]
    .drop_duplicates()
    .shape[0]
)


unique_cell_pairs = int(
    df[
        [
            "sender_cell_type",
            "receiver_cell_type",
        ]
    ]
    .drop_duplicates()
    .shape[0]
)


summary = pd.DataFrame(
    [
        {
            "metric":
                "Raw significant CellChat rows",

            "value":
                len(raw),
        },
        {
            "metric":
                "Full standardized canonical events",

            "value":
                len(df),
        },
        {
            "metric":
                "Individual gene-gene events",

            "value":
                full_gene_gene,
        },
        {
            "metric":
                "Complex/composite events preserved",

            "value":
                full_complex,
        },
        {
            "metric":
                "Comparison-ready canonical events",

            "value":
                len(comparison),
        },
        {
            "metric":
                "Full unique molecular pairs",

            "value":
                unique_lr_full,
        },
        {
            "metric":
                "Comparison-ready unique LR pairs",

            "value":
                unique_lr_comparison,
        },
        {
            "metric":
                "Unique sender-receiver pairs",

            "value":
                unique_cell_pairs,
        },
        {
            "metric":
                "Duplicate canonical keys before collapse",

            "value":
                duplicate_count,
        },
    ]
)


summary.to_csv(
    SUMMARY_FILE,
    index=False,
)


# --------------------------------------------------
# Terminal report
# --------------------------------------------------
print(
    "\n========================================"
)

print(
    "=== CELLCHAT STANDARDIZATION COMPLETE ==="
)

print(
    "========================================"
)


print(
    "\nRaw significant rows:",
    f"{len(raw):,}",
)


print(
    "Full standardized canonical events:",
    f"{len(df):,}",
)


print(
    "Individual gene-gene events:",
    f"{full_gene_gene:,}",
)


print(
    "Complex/composite events preserved:",
    f"{full_complex:,}",
)


print(
    "Comparison-ready canonical events:",
    f"{len(comparison):,}",
)


print(
    "Unique comparison-ready LR pairs:",
    f"{unique_lr_comparison:,}",
)


print(
    "Unique sender-receiver pairs:",
    f"{unique_cell_pairs:,}",
)


print(
    "\nScore range:",
    f"{df['score'].min():.8g}",
    "to",
    f"{df['score'].max():.8g}",
)


print(
    "\nTop 10 comparison-ready events:"
)


top_columns = [
    "sender_cell_type",
    "receiver_cell_type",
    "ligand",
    "receptor",
    "score",
    "pval",
    "method_global_percentile",
]


if len(comparison) > 0:

    print(
        comparison[
            top_columns
        ]
        .head(10)
        .to_string(
            index=False
        )
    )


print(
    "\nSaved:"
)

print(
    STANDARDIZED_FILE
)

print(
    COMPARISON_FILE
)

print(
    SUMMARY_FILE
)