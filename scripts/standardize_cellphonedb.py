from pathlib import Path

import pandas as pd


# --------------------------------------------------
# Paths
# --------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "results"
    / "cellphonedb"
    / "simple_analysis_means_result_prototype.txt"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "results"
    / "cellphonedb"
    / "standardized_interactions.csv"
)

OUTPUT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True,
)


# --------------------------------------------------
# Load CellPhoneDB results
# --------------------------------------------------
df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
)

print("Original shape:", df.shape)


# --------------------------------------------------
# Identify cell-pair columns
#
# CellPhoneDB stores cell-type combinations as:
#
# CellTypeA|CellTypeB
# --------------------------------------------------
pair_columns = [
    column
    for column in df.columns
    if "|" in column
]

print("Cell-pair columns:", len(pair_columns))


# --------------------------------------------------
# CellPhoneDB metadata columns to preserve
# --------------------------------------------------
metadata_columns = [
    "id_cp_interaction",
    "interacting_pair",
    "partner_a",
    "partner_b",
    "gene_a",
    "gene_b",
    "receptor_a",
    "receptor_b",
    "annotation_strategy",
    "is_integrin",
]

metadata_columns = [
    column
    for column in metadata_columns
    if column in df.columns
]


# --------------------------------------------------
# Convert wide CellPhoneDB result to long format
#
# Wide:
# interaction | CellA|CellB | CellA|CellC ...
#
# Long:
# interaction | CellA | CellB | score
# --------------------------------------------------
long_df = df[
    metadata_columns + pair_columns
].melt(
    id_vars=metadata_columns,
    value_vars=pair_columns,
    var_name="cell_pair",
    value_name="score",
)


# --------------------------------------------------
# Separate cell pair into two cell types
# --------------------------------------------------
long_df[
    ["cell_type_a", "cell_type_b"]
] = long_df["cell_pair"].str.split(
    "|",
    n=1,
    expand=True,
    regex=False,
)


# --------------------------------------------------
# Convert interaction value to numeric
# --------------------------------------------------
long_df["score"] = pd.to_numeric(
    long_df["score"],
    errors="coerce",
)


# --------------------------------------------------
# Keep only expressed interactions
# --------------------------------------------------
long_df = long_df[
    long_df["score"].notna()
    & (long_df["score"] > 0)
].copy()


print(
    "Expressed interaction × cell-pair rows:",
    len(long_df),
)


# --------------------------------------------------
# Helper functions
# --------------------------------------------------
def clean_value(value):
    """
    Convert missing values to None and clean strings.
    """

    if pd.isna(value):
        return None

    value = str(value).strip()

    if value == "":
        return None

    return value


def to_bool(value):
    """
    Safely convert CellPhoneDB receptor flags to bool.

    Handles:
        True / False
        "True" / "False"
        1 / 0
    """

    if pd.isna(value):
        return False

    if isinstance(value, bool):
        return value

    if isinstance(value, (int, float)):
        return bool(value)

    return str(value).strip().lower() in {
        "true",
        "1",
        "yes",
        "t",
    }


def preferred_name(gene, partner):
    """
    Prefer human-readable gene symbols.

    Example:

        gene:
            APP

        CellPhoneDB partner:
            simple:P05067

        platform display:
            APP

    If a gene symbol is unavailable, as may happen
    for complexes, retain the original CellPhoneDB
    partner identifier.
    """

    gene = clean_value(gene)
    partner = clean_value(partner)

    if gene is not None:
        return gene

    return partner


# --------------------------------------------------
# Determine biological direction
#
# CellPhoneDB contains receptor_a and receptor_b.
#
# If B is receptor:
#
#     Cell A = sender
#     Partner A = ligand
#
#     Cell B = receiver
#     Partner B = receptor
#
#
# If A is receptor:
#
#     Cell B = sender
#     Partner B = ligand
#
#     Cell A = receiver
#     Partner A = receptor
#
#
# Other cases are not considered safely directional.
# --------------------------------------------------
def assign_direction(row):

    receptor_a = to_bool(
        row.get("receptor_a", False)
    )

    receptor_b = to_bool(
        row.get("receptor_b", False)
    )


    gene_a = row.get("gene_a")
    gene_b = row.get("gene_b")

    partner_a = row.get("partner_a")
    partner_b = row.get("partner_b")


    # --------------------------------------------------
    # B is receptor
    #
    # A -> B
    # --------------------------------------------------
    if receptor_b and not receptor_a:

        return pd.Series(
            {
                "sender_cell_type":
                    row["cell_type_a"],

                "receiver_cell_type":
                    row["cell_type_b"],

                "ligand":
                    preferred_name(
                        gene_a,
                        partner_a,
                    ),

                "receptor":
                    preferred_name(
                        gene_b,
                        partner_b,
                    ),

                "ligand_partner_id":
                    clean_value(
                        partner_a
                    ),

                "receptor_partner_id":
                    clean_value(
                        partner_b
                    ),

                "interaction_type":
                    "ligand-receptor",
            }
        )


    # --------------------------------------------------
    # A is receptor
    #
    # B -> A
    # --------------------------------------------------
    if receptor_a and not receptor_b:

        return pd.Series(
            {
                "sender_cell_type":
                    row["cell_type_b"],

                "receiver_cell_type":
                    row["cell_type_a"],

                "ligand":
                    preferred_name(
                        gene_b,
                        partner_b,
                    ),

                "receptor":
                    preferred_name(
                        gene_a,
                        partner_a,
                    ),

                "ligand_partner_id":
                    clean_value(
                        partner_b
                    ),

                "receptor_partner_id":
                    clean_value(
                        partner_a
                    ),

                "interaction_type":
                    "ligand-receptor",
            }
        )


    # --------------------------------------------------
    # No safe ligand -> receptor direction
    # --------------------------------------------------
    return pd.Series(
        {
            "sender_cell_type": None,
            "receiver_cell_type": None,

            "ligand": None,
            "receptor": None,

            "ligand_partner_id": None,
            "receptor_partner_id": None,

            "interaction_type":
                "non-directional-or-complex",
        }
    )


# --------------------------------------------------
# Apply biological direction
# --------------------------------------------------
print(
    "Assigning ligand–receptor direction..."
)

direction_df = long_df.apply(
    assign_direction,
    axis=1,
)

long_df = pd.concat(
    [
        long_df,
        direction_df,
    ],
    axis=1,
)


# --------------------------------------------------
# Add platform metadata
# --------------------------------------------------
long_df["method"] = "CellPhoneDB"

# Stamp the real dataset id. The scripts run inside the run folder
# datasets/<id>/runs/<slug>/, so PROJECT_ROOT is that run folder and its
# grandparent is the dataset id. (Previously this was hard-coded to the
# original development dataset "endometrium_example", which made the consensus
# step drop every row for any other dataset.)
try:
    _dataset_id = PROJECT_ROOT.parent.parent.name or "dataset"
except Exception:
    _dataset_id = "dataset"
long_df["dataset_id"] = _dataset_id


# --------------------------------------------------
# Final common platform schema
# --------------------------------------------------
final_columns = [
    "dataset_id",

    "method",

    "id_cp_interaction",
    "interacting_pair",

    "cell_type_a",
    "cell_type_b",

    "partner_a",
    "partner_b",

    "gene_a",
    "gene_b",

    "receptor_a",
    "receptor_b",

    "sender_cell_type",
    "receiver_cell_type",

    "ligand",
    "receptor",

    "ligand_partner_id",
    "receptor_partner_id",

    "interaction_type",

    "annotation_strategy",
    "is_integrin",

    "score",
]


final_columns = [
    column
    for column in final_columns
    if column in long_df.columns
]


standardized_df = long_df[
    final_columns
].copy()


# --------------------------------------------------
# Sort by computational evidence
# --------------------------------------------------
standardized_df = (
    standardized_df
    .sort_values(
        "score",
        ascending=False,
    )
    .reset_index(
        drop=True
    )
)


# --------------------------------------------------
# Save standardized result
# --------------------------------------------------
standardized_df.to_csv(
    OUTPUT_FILE,
    index=False,
)


# --------------------------------------------------
# Summary
# --------------------------------------------------
print(
    "\nStandardization completed."
)

print(
    "Final shape:",
    standardized_df.shape,
)


print(
    "\nInteraction types:"
)

print(
    standardized_df[
        "interaction_type"
    ].value_counts()
)


# --------------------------------------------------
# Directional interaction summary
# --------------------------------------------------
directional_df = standardized_df[
    standardized_df[
        "interaction_type"
    ]
    == "ligand-receptor"
].copy()


print(
    "\nDirectional ligand-receptor rows:",
    len(directional_df),
)


print(
    "\nTop 10 directional interactions:"
)


preview_columns = [
    "sender_cell_type",
    "receiver_cell_type",
    "ligand",
    "receptor",
    "ligand_partner_id",
    "receptor_partner_id",
    "interacting_pair",
    "score",
]


preview_columns = [
    column
    for column in preview_columns
    if column in directional_df.columns
]


print(
    directional_df[
        preview_columns
    ]
    .head(10)
    .to_string(
        index=False
    )
)


print(
    f"\nSaved to:\n{OUTPUT_FILE}"
)