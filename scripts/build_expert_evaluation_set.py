from pathlib import Path
import argparse
import hashlib
import sqlite3

import numpy as np
import pandas as pd


# --------------------------------------------------
# Project paths
# --------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_CONSENSUS_FILE = (
    PROJECT_ROOT
    / "results"
    / "consensus"
    / "consensus_cellphonedb__natmi__cellchat.csv"
)

DEFAULT_DATABASE_FILE = (
    PROJECT_ROOT
    / "database"
    / "ccc_platform.db"
)

DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "evaluation"
)


# --------------------------------------------------
# Operational stratum thresholds
#
# IMPORTANT:
# These are dataset-calibrated operational criteria.
# They are NOT claimed to be biologically optimal.
# --------------------------------------------------
ROBUST_CONSENSUS_MIN_PERCENTILE = 95.0
ROBUST_CONSENSUS_MAX_DISAGREEMENT = 10.0

EXTREME_DISAGREEMENT_MIN = 80.0

EXTREME_SINGLE_MIN_PERCENTILE = 99.9

CONTROL_CONSENSUS_MIN = 40.0
CONTROL_CONSENSUS_MAX = 80.0
CONTROL_DISAGREEMENT_MAX = 40.0


# --------------------------------------------------
# Canonical event helper
# --------------------------------------------------
def canonical_key(
    sender,
    receiver,
    ligand,
    receptor,
):
    return (
        str(sender)
        + "||"
        + str(receiver)
        + "||"
        + str(ligand).upper()
        + "||"
        + str(receptor).upper()
    )


# --------------------------------------------------
# Stable anonymized ID
# --------------------------------------------------
def make_event_hash(
    interaction_key,
):
    digest = hashlib.sha256(
        str(interaction_key).encode("utf-8")
    ).hexdigest()

    return digest[:12].upper()


# --------------------------------------------------
# Load reviewed canonical events
# --------------------------------------------------
def load_reviewed_event_keys(
    database_file,
):
    if not database_file.exists():
        print(
            "WARNING: database does not exist:",
            database_file,
        )
        return set()

    connection = sqlite3.connect(
        database_file
    )

    query = """
        SELECT
            sender_cell_type,
            receiver_cell_type,
            ligand,
            receptor
        FROM expert_feedback
        WHERE
            sender_cell_type IS NOT NULL
            AND receiver_cell_type IS NOT NULL
            AND ligand IS NOT NULL
            AND receptor IS NOT NULL
    """

    feedback = pd.read_sql_query(
        query,
        connection,
    )

    connection.close()

    if feedback.empty:
        return set()

    feedback[
        "interaction_key"
    ] = feedback.apply(
        lambda row:
            canonical_key(
                row[
                    "sender_cell_type"
                ],
                row[
                    "receiver_cell_type"
                ],
                row[
                    "ligand"
                ],
                row[
                    "receptor"
                ],
            ),
        axis=1,
    )

    return set(
        feedback[
            "interaction_key"
        ].astype(str)
    )


# --------------------------------------------------
# Validate consensus input
# --------------------------------------------------
def validate_consensus(
    consensus,
):
    required_columns = {
        "interaction_key",
        "sender_cell_type",
        "receiver_cell_type",
        "ligand",
        "receptor",
        "method_count",
        "total_selected_methods",
        "support_fraction",
        "method_support",
        "consensus_percentile",
        "rank_disagreement",
        "max_method_percentile",
    }

    missing = (
        required_columns
        - set(
            consensus.columns
        )
    )

    if missing:
        raise ValueError(
            "Consensus file is missing required "
            f"columns: {sorted(missing)}"
        )

    duplicate_keys = (
        consensus[
            "interaction_key"
        ]
        .duplicated()
        .sum()
    )

    if duplicate_keys:
        raise ValueError(
            "Consensus file contains "
            f"{duplicate_keys} duplicate "
            "canonical interaction keys."
        )


# --------------------------------------------------
# Define mutually exclusive evaluation strata
# --------------------------------------------------
def assign_strata(
    consensus,
):
    work = consensus.copy()

    # Numerical safety
    numeric_columns = [
        "method_count",
        "total_selected_methods",
        "support_fraction",
        "consensus_percentile",
        "rank_disagreement",
        "max_method_percentile",
    ]

    for column in numeric_columns:
        work[column] = pd.to_numeric(
            work[column],
            errors="coerce",
        )

    work[
        "evaluation_stratum"
    ] = pd.NA

    # --------------------------------------------------
    # 1. Robust full consensus
    # --------------------------------------------------
    robust_mask = (
        (
            work[
                "method_count"
            ]
            == work[
                "total_selected_methods"
            ]
        )
        & (
            work[
                "consensus_percentile"
            ]
            >= ROBUST_CONSENSUS_MIN_PERCENTILE
        )
        & (
            work[
                "rank_disagreement"
            ]
            <= ROBUST_CONSENSUS_MAX_DISAGREEMENT
        )
    )

    work.loc[
        robust_mask,
        "evaluation_stratum",
    ] = "robust_full_consensus"

    # --------------------------------------------------
    # 2. Extreme cross-method disagreement
    #
    # Requires detection by at least two methods.
    # This is presence agreement combined with
    # major rank disagreement.
    # --------------------------------------------------
    disagreement_mask = (
        work[
            "evaluation_stratum"
        ].isna()
        & (
            work[
                "method_count"
            ]
            >= 2
        )
        & (
            work[
                "rank_disagreement"
            ]
            >= EXTREME_DISAGREEMENT_MIN
        )
    )

    work.loc[
        disagreement_mask,
        "evaluation_stratum",
    ] = "extreme_cross_method_disagreement"

    # --------------------------------------------------
    # 3. Extreme single-method discovery
    # --------------------------------------------------
    single_mask = (
        work[
            "evaluation_stratum"
        ].isna()
        & (
            work[
                "method_count"
            ]
            == 1
        )
        & (
            work[
                "max_method_percentile"
            ]
            >= EXTREME_SINGLE_MIN_PERCENTILE
        )
    )

    work.loc[
        single_mask,
        "evaluation_stratum",
    ] = "extreme_single_method_discovery"

    # --------------------------------------------------
    # 4. Intermediate shared-event control
    #
    # This deliberately avoids both:
    # - very high consensus
    # - extreme disagreement
    #
    # It is not a biological negative control.
    # It is a computationally intermediate control.
    # --------------------------------------------------
    control_mask = (
        work[
            "evaluation_stratum"
        ].isna()
        & (
            work[
                "method_count"
            ]
            >= 2
        )
        & (
            work[
                "consensus_percentile"
            ]
            >= CONTROL_CONSENSUS_MIN
        )
        & (
            work[
                "consensus_percentile"
            ]
            <= CONTROL_CONSENSUS_MAX
        )
        & (
            work[
                "rank_disagreement"
            ]
            < CONTROL_DISAGREEMENT_MAX
        )
    )

    work.loc[
        control_mask,
        "evaluation_stratum",
    ] = "intermediate_shared_control"

    return work


# --------------------------------------------------
# Stratified sampling
# --------------------------------------------------
def stratified_sample(
    candidates,
    n_per_stratum,
    seed,
):
    strata = [
        "robust_full_consensus",
        "extreme_cross_method_disagreement",
        "extreme_single_method_discovery",
        "intermediate_shared_control",
    ]

    sampled_frames = []

    for index, stratum in enumerate(
        strata
    ):
        subset = candidates[
            candidates[
                "evaluation_stratum"
            ]
            == stratum
        ].copy()

        available = len(
            subset
        )

        requested = int(
            n_per_stratum
        )

        take = min(
            available,
            requested,
        )

        if take == 0:
            print(
                f"WARNING: no candidates for {stratum}"
            )
            continue

        if take < requested:
            print(
                "WARNING:",
                stratum,
                "contains only",
                available,
                "eligible events; requested",
                requested,
            )

        sampled = subset.sample(
            n=take,
            replace=False,
            random_state=(
                seed
                + index
            ),
        )

        sampled_frames.append(
            sampled
        )

    if not sampled_frames:
        raise RuntimeError(
            "No events could be sampled."
        )

    evaluation = pd.concat(
        sampled_frames,
        ignore_index=True,
    )

    # Global randomized order.
    evaluation = (
        evaluation
        .sample(
            frac=1.0,
            random_state=seed + 1000,
        )
        .reset_index(
            drop=True
        )
    )

    return evaluation


# --------------------------------------------------
# Add blinded evaluation identifiers
# --------------------------------------------------
def add_evaluation_ids(
    evaluation,
):
    result = (
        evaluation.copy()
    )

    result[
        "event_hash"
    ] = (
        result[
            "interaction_key"
        ]
        .map(
            make_event_hash
        )
    )

    result[
        "evaluation_id"
    ] = [
        f"EVAL-{number:04d}"
        for number
        in range(
            1,
            len(result) + 1,
        )
    ]

    return result


# --------------------------------------------------
# Build Stage 1 blinded file
# --------------------------------------------------
def build_stage1_file(
    evaluation,
):
    columns = [
        "evaluation_id",
        "sender_cell_type",
        "receiver_cell_type",
        "ligand",
        "receptor",
    ]

    stage1 = (
        evaluation[
            columns
        ]
        .copy()
    )

    stage1[
        "pre_evidence_decision"
    ] = ""

    stage1[
        "pre_evidence_biological_relevance"
    ] = ""

    stage1[
        "pre_evidence_confidence"
    ] = ""

    stage1[
        "pre_evidence_comment"
    ] = ""

    return stage1


# --------------------------------------------------
# Build Stage 2 evidence-reveal file
# --------------------------------------------------
def build_stage2_file(
    evaluation,
):
    base_columns = [
        "evaluation_id",
        "sender_cell_type",
        "receiver_cell_type",
        "ligand",
        "receptor",
        "method_count",
        "total_selected_methods",
        "support_fraction",
        "method_support",
        "consensus_percentile",
        "rank_disagreement",
    ]

    method_columns = []

    for method_slug in [
        "cellphonedb",
        "natmi",
        "cellchat",
    ]:
        for suffix in [
            "present",
            "score",
            "score_name",
            "global_percentile",
            "cell_pair_percentile",
        ]:
            column = (
                f"{method_slug}_{suffix}"
            )

            if column in evaluation.columns:
                method_columns.append(
                    column
                )

    selected_columns = [
        column
        for column
        in (
            base_columns
            + method_columns
        )
        if column
        in evaluation.columns
    ]

    stage2 = (
        evaluation[
            selected_columns
        ]
        .copy()
    )

    stage2[
        "post_evidence_decision"
    ] = ""

    stage2[
        "post_evidence_biological_relevance"
    ] = ""

    stage2[
        "post_evidence_confidence"
    ] = ""

    stage2[
        "post_evidence_comment"
    ] = ""

    return stage2


# --------------------------------------------------
# Sampling summary
# --------------------------------------------------
def build_summary(
    all_events,
    eligible_candidates,
    evaluation,
    reviewed_count,
    requested_per_stratum,
    seed,
):
    strata = [
        "robust_full_consensus",
        "extreme_cross_method_disagreement",
        "extreme_single_method_discovery",
        "intermediate_shared_control",
    ]

    rows = []

    for stratum in strata:
        candidate_count = int(
            (
                eligible_candidates[
                    "evaluation_stratum"
                ]
                == stratum
            ).sum()
        )

        sampled_count = int(
            (
                evaluation[
                    "evaluation_stratum"
                ]
                == stratum
            ).sum()
        )

        rows.append(
            {
                "stratum":
                    stratum,

                "eligible_candidates":
                    candidate_count,

                "requested_sample_size":
                    requested_per_stratum,

                "sampled_events":
                    sampled_count,

                "random_seed":
                    seed,
            }
        )

    rows.append(
        {
            "stratum":
                "__overall__",

            "eligible_candidates":
                int(
                    len(
                        eligible_candidates
                    )
                ),

            "requested_sample_size":
                (
                    requested_per_stratum
                    * 4
                ),

            "sampled_events":
                int(
                    len(
                        evaluation
                    )
                ),

            "random_seed":
                seed,
        }
    )

    summary = pd.DataFrame(
        rows
    )

    summary[
        "total_consensus_events"
    ] = len(
        all_events
    )

    summary[
        "previously_reviewed_events_excluded"
    ] = reviewed_count

    return summary


# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a reproducible stratified "
            "expert-evaluation set from the "
            "3-method CCC consensus."
        )
    )

    parser.add_argument(
        "--consensus-file",
        type=Path,
        default=DEFAULT_CONSENSUS_FILE,
    )

    parser.add_argument(
        "--database-file",
        type=Path,
        default=DEFAULT_DATABASE_FILE,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
    )

    parser.add_argument(
        "--n-per-stratum",
        type=int,
        default=25,
        help=(
            "Number of events sampled from each "
            "computational stratum. Default: 25."
        ),
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20260914,
        help=(
            "Random seed for reproducible sampling."
        ),
    )

    args = parser.parse_args()

    if args.n_per_stratum <= 0:
        raise ValueError(
            "--n-per-stratum must be > 0."
        )

    if not args.consensus_file.exists():
        raise FileNotFoundError(
            f"Consensus file not found: "
            f"{args.consensus_file}"
        )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------
    # Load consensus
    # --------------------------------------------------
    consensus = pd.read_csv(
        args.consensus_file
    )

    validate_consensus(
        consensus
    )

    print(
        "Loaded consensus events:",
        len(consensus),
    )

    # --------------------------------------------------
    # Assign computational strata
    # --------------------------------------------------
    stratified = assign_strata(
        consensus
    )

    # --------------------------------------------------
    # Exclude any previously reviewed event
    # --------------------------------------------------
    reviewed_keys = (
        load_reviewed_event_keys(
            args.database_file
        )
    )

    reviewed_mask = (
        stratified[
            "interaction_key"
        ]
        .astype(str)
        .isin(
            reviewed_keys
        )
    )

    reviewed_count = int(
        reviewed_mask.sum()
    )

    eligible = stratified[
        ~reviewed_mask
    ].copy()

    # Keep only the four evaluation strata.
    eligible = eligible[
        eligible[
            "evaluation_stratum"
        ].notna()
    ].copy()

    print(
        "Previously reviewed events excluded:",
        reviewed_count,
    )

    print(
        "\nEligible candidate pools:"
    )

    print(
        eligible[
            "evaluation_stratum"
        ]
        .value_counts()
        .to_string()
    )

    # --------------------------------------------------
    # Sample
    # --------------------------------------------------
    evaluation = stratified_sample(
        candidates=eligible,
        n_per_stratum=(
            args.n_per_stratum
        ),
        seed=args.seed,
    )

    evaluation = add_evaluation_ids(
        evaluation
    )

    # --------------------------------------------------
    # Safety checks
    # --------------------------------------------------
    if (
        evaluation[
            "interaction_key"
        ]
        .duplicated()
        .any()
    ):
        raise RuntimeError(
            "Duplicate canonical events appeared "
            "in the sampled evaluation set."
        )

    if (
        set(
            evaluation[
                "interaction_key"
            ]
        )
        & reviewed_keys
    ):
        raise RuntimeError(
            "A previously reviewed event leaked "
            "into the evaluation set."
        )

    # --------------------------------------------------
    # Full internal researcher file
    # --------------------------------------------------
    preferred_internal_columns = [
        "evaluation_id",
        "event_hash",
        "evaluation_stratum",
        "interaction_key",
        "sender_cell_type",
        "receiver_cell_type",
        "ligand",
        "receptor",
        "method_count",
        "total_selected_methods",
        "support_fraction",
        "support_label",
        "method_support",
        "support_category",
        "consensus_percentile",
        "rank_disagreement",
        "max_method_percentile",
        "min_method_percentile",
        "mean_available_percentile",
    ]

    method_columns = [
        column
        for column
        in evaluation.columns
        if (
            column.startswith(
                "cellphonedb_"
            )
            or column.startswith(
                "natmi_"
            )
            or column.startswith(
                "cellchat_"
            )
        )
    ]

    internal_columns = [
        column
        for column
        in (
            preferred_internal_columns
            + method_columns
        )
        if column
        in evaluation.columns
    ]

    full_internal = (
        evaluation[
            internal_columns
        ]
        .copy()
    )

    # --------------------------------------------------
    # Stage-specific files
    # --------------------------------------------------
    stage1 = build_stage1_file(
        evaluation
    )

    stage2 = build_stage2_file(
        evaluation
    )

    summary = build_summary(
        all_events=stratified,
        eligible_candidates=eligible,
        evaluation=evaluation,
        reviewed_count=reviewed_count,
        requested_per_stratum=(
            args.n_per_stratum
        ),
        seed=args.seed,
    )

    # --------------------------------------------------
    # Save
    # --------------------------------------------------
    full_file = (
        args.output_dir
        / "expert_evaluation_set_full.csv"
    )

    stage1_file = (
        args.output_dir
        / "expert_evaluation_stage1_blinded.csv"
    )

    stage2_file = (
        args.output_dir
        / "expert_evaluation_stage2_evidence.csv"
    )

    summary_file = (
        args.output_dir
        / "expert_evaluation_sampling_summary.csv"
    )

    full_internal.to_csv(
        full_file,
        index=False,
    )

    stage1.to_csv(
        stage1_file,
        index=False,
    )

    stage2.to_csv(
        stage2_file,
        index=False,
    )

    summary.to_csv(
        summary_file,
        index=False,
    )

    print(
        "\nSampled evaluation strata:"
    )

    print(
        evaluation[
            "evaluation_stratum"
        ]
        .value_counts()
        .to_string()
    )

    print(
        "\nTotal evaluation events:",
        len(evaluation),
    )

    print(
        "\nSaved:"
    )

    print(
        " ",
        full_file,
    )

    print(
        " ",
        stage1_file,
    )

    print(
        " ",
        stage2_file,
    )

    print(
        " ",
        summary_file,
    )

    print(
        "\nIMPORTANT:"
    )

    print(
        "Stage 1 is intended to be blinded to "
        "method support, scores and stratum."
    )

    print(
        "Do not begin a real expert study using "
        "the current Feedback UI until a blinded "
        "review mode has been implemented."
    )


if __name__ == "__main__":
    main()