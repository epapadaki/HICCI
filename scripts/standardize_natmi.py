from pathlib import Path

import pandas as pd


# --------------------------------------------------
# Project paths
# --------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "results"
    / "natmi"
    / "Edges_lrc2p.csv"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "results"
    / "natmi"
    / "standardized_interactions.csv"
)

COMPARISON_FILE = (
    PROJECT_ROOT
    / "results"
    / "natmi"
    / "standardized_interactions_comparison_ready.csv"
)


# --------------------------------------------------
# Settings
# --------------------------------------------------
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

METHOD_NAME = "NATMI"

COMMON_DETECTION_THRESHOLD = 0.10


# --------------------------------------------------
# Validate input
# --------------------------------------------------
if not INPUT_FILE.exists():

    raise FileNotFoundError(
        f"NATMI result file not found:\n{INPUT_FILE}"
    )


# --------------------------------------------------
# Load NATMI output
# --------------------------------------------------
print("Loading NATMI results...")

df = pd.read_csv(
    INPUT_FILE
)


print(
    f"Original shape: "
    f"{df.shape[0]:,} rows × "
    f"{df.shape[1]:,} columns"
)


# --------------------------------------------------
# Required NATMI columns
# --------------------------------------------------
required_columns = [
    "Sending cluster",
    "Ligand symbol",
    "Receptor symbol",
    "Target cluster",

    "Ligand-expressing cells",
    "Ligand detection rate",
    "Ligand average expression value",
    "Ligand total expression value",

    "Receptor-expressing cells",
    "Receptor detection rate",
    "Receptor average expression value",
    "Receptor total expression value",

    "Edge average expression weight",
    "Edge average expression derived specificity",

    "Edge total expression weight",
    "Edge total expression derived specificity",
]


missing_columns = [
    column
    for column in required_columns
    if column not in df.columns
]


if missing_columns:

    raise ValueError(
        "Missing required NATMI columns:\n"
        + "\n".join(missing_columns)
    )


# --------------------------------------------------
# Rename NATMI fields into platform-friendly names
# --------------------------------------------------
rename_map = {

    "Sending cluster":
        "sender_cell_type",

    "Target cluster":
        "receiver_cell_type",

    "Ligand symbol":
        "ligand",

    "Receptor symbol":
        "receptor",

    "Ligand-expressing cells":
        "ligand_expressing_cells",

    "Ligand detection rate":
        "ligand_detection_rate",

    "Ligand average expression value":
        "ligand_average_expression",

    "Ligand total expression value":
        "ligand_total_expression",

    "Ligand derived specificity of average expression value":
        "ligand_average_expression_specificity",

    "Ligand derived specificity of total expression value":
        "ligand_total_expression_specificity",

    "Receptor-expressing cells":
        "receptor_expressing_cells",

    "Receptor detection rate":
        "receptor_detection_rate",

    "Receptor average expression value":
        "receptor_average_expression",

    "Receptor total expression value":
        "receptor_total_expression",

    "Receptor derived specificity of average expression value":
        "receptor_average_expression_specificity",

    "Receptor derived specificity of total expression value":
        "receptor_total_expression_specificity",

    "Edge average expression weight":
        "score",

    "Edge average expression derived specificity":
        "edge_average_expression_specificity",

    "Edge total expression weight":
        "edge_total_expression_weight",

    "Edge total expression derived specificity":
        "edge_total_expression_specificity",
}


df = df.rename(
    columns=rename_map
)


# --------------------------------------------------
# Clean textual fields
# --------------------------------------------------
text_columns = [
    "sender_cell_type",
    "receiver_cell_type",
    "ligand",
    "receptor",
]


for column in text_columns:

    df[column] = (
        df[column]
        .astype("string")
        .str.strip()
    )


# --------------------------------------------------
# Convert numerical fields
# --------------------------------------------------
numeric_columns = [
    "ligand_expressing_cells",
    "ligand_detection_rate",
    "ligand_average_expression",
    "ligand_total_expression",
    "ligand_average_expression_specificity",
    "ligand_total_expression_specificity",

    "receptor_expressing_cells",
    "receptor_detection_rate",
    "receptor_average_expression",
    "receptor_total_expression",
    "receptor_average_expression_specificity",
    "receptor_total_expression_specificity",

    "score",
    "edge_average_expression_specificity",

    "edge_total_expression_weight",
    "edge_total_expression_specificity",
]


for column in numeric_columns:

    if column in df.columns:

        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )


# --------------------------------------------------
# Remove invalid rows
# --------------------------------------------------
before_cleaning = len(df)


df = df[
    df["sender_cell_type"].notna()
    & df["receiver_cell_type"].notna()
    & df["ligand"].notna()
    & df["receptor"].notna()
    & df["score"].notna()
].copy()


df = df[
    df["score"] > 0
].copy()


print(
    f"Valid positive-score interactions: "
    f"{len(df):,}"
)


removed = (
    before_cleaning
    - len(df)
)


print(
    f"Removed invalid/non-positive rows: "
    f"{removed:,}"
)


# --------------------------------------------------
# Add common platform fields
# --------------------------------------------------
df["dataset_id"] = DATASET_ID

df["method"] = METHOD_NAME

df["interaction_type"] = (
    "ligand-receptor"
)


# --------------------------------------------------
# Canonical molecular interaction label
# --------------------------------------------------
df["interacting_pair"] = (
    df["ligand"].astype(str)
    + "_"
    + df["receptor"].astype(str)
)


# --------------------------------------------------
# Canonical interaction key
#
# This will later allow:
#
# CellPhoneDB vs NATMI
#
# same sender
# same receiver
# same ligand
# same receptor
# --------------------------------------------------
df["interaction_key"] = (
    df["sender_cell_type"].astype(str)
    + "||"
    + df["receiver_cell_type"].astype(str)
    + "||"
    + df["ligand"].astype(str)
    + "||"
    + df["receptor"].astype(str)
)


# --------------------------------------------------
# Score metadata
#
# Raw scores from different CCC methods must NOT
# be compared directly because they have different
# mathematical definitions and scales.
# --------------------------------------------------
df["score_name"] = (
    "edge_average_expression_weight"
)

df["score_higher_is_better"] = True


# --------------------------------------------------
# Common detection threshold
#
# CellPhoneDB was run using threshold = 0.10.
#
# We preserve ALL NATMI interactions, but mark
# interactions that satisfy the same broad
# expression prevalence criterion for future
# cross-method comparison.
# --------------------------------------------------
df[
    "passes_common_detection_threshold"
] = (
    (
        df["ligand_detection_rate"]
        >= COMMON_DETECTION_THRESHOLD
    )
    &
    (
        df["receptor_detection_rate"]
        >= COMMON_DETECTION_THRESHOLD
    )
)


# --------------------------------------------------
# Check canonical uniqueness
# --------------------------------------------------
duplicate_mask = df[
    "interaction_key"
].duplicated(
    keep=False
)


duplicate_count = int(
    duplicate_mask.sum()
)


print(
    f"Duplicate canonical interaction keys: "
    f"{duplicate_count:,}"
)


if duplicate_count > 0:

    print(
        "\nWARNING:"
        "\nDuplicate sender/receiver/ligand/receptor "
        "combinations were found."
    )

    print(
        df.loc[
            duplicate_mask,
            [
                "interaction_key",
                "score",
            ],
        ]
        .head(20)
        .to_string(index=False)
    )

    raise ValueError(
        "Canonical NATMI interactions are not unique. "
        "Inspect duplicates before continuing."
    )


# --------------------------------------------------
# Calculate NATMI-specific computational percentiles
#
# These are useful because NATMI raw scores cannot
# later be directly compared to CellPhoneDB raw scores.
# --------------------------------------------------
df[
    "method_global_percentile"
] = (
    df["score"]
    .rank(
        method="average",
        pct=True,
    )
    * 100
)


df[
    "method_cell_pair_percentile"
] = (
    df
    .groupby(
        [
            "sender_cell_type",
            "receiver_cell_type",
        ]
    )["score"]
    .rank(
        method="average",
        pct=True,
    )
    * 100
)


# --------------------------------------------------
# Sort strongest NATMI interactions first
# --------------------------------------------------
df = (
    df
    .sort_values(
        "score",
        ascending=False,
    )
    .reset_index(
        drop=True
    )
)


# --------------------------------------------------
# Final standardized column order
# --------------------------------------------------
final_columns = [

    # General platform fields
    "dataset_id",
    "method",

    "interaction_key",
    "interacting_pair",

    "interaction_type",

    # Direction
    "sender_cell_type",
    "receiver_cell_type",

    # Molecular pair
    "ligand",
    "receptor",

    # Generic method score
    "score",
    "score_name",
    "score_higher_is_better",

    # Method-normalized ranks
    "method_global_percentile",
    "method_cell_pair_percentile",

    # Harmonized comparison flag
    "passes_common_detection_threshold",

    # NATMI ligand evidence
    "ligand_expressing_cells",
    "ligand_detection_rate",
    "ligand_average_expression",
    "ligand_total_expression",
    "ligand_average_expression_specificity",
    "ligand_total_expression_specificity",

    # NATMI receptor evidence
    "receptor_expressing_cells",
    "receptor_detection_rate",
    "receptor_average_expression",
    "receptor_total_expression",
    "receptor_average_expression_specificity",
    "receptor_total_expression_specificity",

    # NATMI edge evidence
    "edge_average_expression_specificity",
    "edge_total_expression_weight",
    "edge_total_expression_specificity",
]


final_columns = [
    column
    for column in final_columns
    if column in df.columns
]


standardized_df = df[
    final_columns
].copy()


# --------------------------------------------------
# Comparison-ready subset
# --------------------------------------------------
comparison_df = standardized_df[
    standardized_df[
        "passes_common_detection_threshold"
    ]
].copy()


# --------------------------------------------------
# Save files
# --------------------------------------------------
standardized_df.to_csv(
    OUTPUT_FILE,
    index=False,
)


comparison_df.to_csv(
    COMPARISON_FILE,
    index=False,
)


# --------------------------------------------------
# Summary
# --------------------------------------------------
print(
    "\nNATMI standardization completed."
)

print(
    f"\nAll standardized NATMI interactions: "
    f"{len(standardized_df):,}"
)

print(
    f"Comparison-ready interactions "
    f"(ligand & receptor detection >= "
    f"{COMMON_DETECTION_THRESHOLD:.2f}): "
    f"{len(comparison_df):,}"
)


print(
    "\nTop 10 NATMI interactions:"
)


preview_columns = [
    "sender_cell_type",
    "receiver_cell_type",
    "ligand",
    "receptor",
    "score",
    "method_global_percentile",
    "ligand_detection_rate",
    "receptor_detection_rate",
    "passes_common_detection_threshold",
]


print(
    standardized_df[
        preview_columns
    ]
    .head(10)
    .to_string(
        index=False
    )
)


print(
    "\nFiles saved:"
)

print(
    OUTPUT_FILE
)

print(
    COMPARISON_FILE
)