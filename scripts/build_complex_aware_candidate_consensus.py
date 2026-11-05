from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

CANDIDATE_OUTPUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "consensus_candidate_complex_aware"
)

METHOD_CONFIG = {
    "CellPhoneDB": {
        "path": (
            PROJECT_ROOT
            / "results"
            / "cellphonedb"
            / "standardized_interactions_complex_aware.csv"
        ),
        "interaction_key_column": (
            "complex_aware_event_key_candidate"
        ),
    },
    "NATMI": {
        "path": (
            PROJECT_ROOT
            / "results"
            / "natmi"
            / "standardized_interactions_comparison_ready.csv"
        ),
        "interaction_key_column": None,
    },
    "CellChat": {
        "path": (
            PROJECT_ROOT
            / "results"
            / "cellchat"
            / "cellchat_interactions_complex_aware.csv"
        ),
        "interaction_key_column": (
            "complex_aware_event_key_candidate"
        ),
    },
    # ICELLNET (v2): curated pairs with complexes (keys built by standardize_icellnet.py).
    "ICELLNET": {
        "path": (
            PROJECT_ROOT
            / "results"
            / "icellnet"
            / "standardized_interactions.csv"
        ),
        "interaction_key_column": (
            "complex_aware_event_key_candidate"
        ),
    },
    # SingleCellSignalR (v1): gene-level ligand-receptor pairs, LRscore.
    "SingleCellSignalR": {
        "path": (
            PROJECT_ROOT
            / "results"
            / "singlecellsignalr"
            / "standardized_interactions.csv"
        ),
        "interaction_key_column": None,
    },
}


def clean_text(series: pd.Series) -> pd.Series:
    out = (
        series
        .astype("string")
        .str.strip()
    )

    out = out.replace(
        {
            "": pd.NA,
            "nan": pd.NA,
            "None": pd.NA,
            "<NA>": pd.NA,
        }
    )

    return out


def parse_bool(value, default=True):
    if pd.isna(value):
        return default

    if isinstance(value, str):
        return (
            value
            .strip()
            .lower()
            in {
                "true",
                "1",
                "yes",
            }
        )

    return bool(value)


def method_slug(method_name: str) -> str:
    return (
        method_name
        .strip()
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
    )


def create_gene_key(df: pd.DataFrame) -> pd.Series:
    ligand_entity = (
        "GENE:"
        + df["ligand"].astype(str).str.upper()
    )

    receptor_entity = (
        "GENE:"
        + df["receptor"].astype(str).str.upper()
    )

    return (
        df["sender_cell_type"].astype(str)
        + "||"
        + df["receiver_cell_type"].astype(str)
        + "||"
        + ligand_entity
        + "||"
        + receptor_entity
    )


def split_candidate_key(
    interaction_key: pd.Series,
) -> pd.DataFrame:
    columns = [
        "key_sender",
        "key_receiver",
        "ligand_entity_key",
        "receptor_entity_key",
    ]

    # No events for this method (e.g. everything was filtered out upstream):
    # return an empty, correctly-shaped frame instead of failing on the split.
    if len(interaction_key) == 0:
        return pd.DataFrame(columns=columns)

    parts = (
        interaction_key
        .astype("string")
        .str.split(
            "||",
            n=3,
            expand=True,
            regex=False,
        )
    )

    if parts.shape[1] != 4:
        example = next(
            (str(v) for v in interaction_key.tolist() if v is not None and str(v) != "<NA>"),
            "<none>",
        )
        raise ValueError(
            "Candidate interaction keys must contain exactly four logical "
            "fields separated by '||' (sender||receiver||ligand_entity||"
            "receptor_entity). The keys for this method do not — they are "
            "empty or malformed, so standardisation produced no usable "
            f"candidate events. Example key seen: {example!r}"
        )

    parts.columns = [
        "key_sender",
        "key_receiver",
        "ligand_entity_key",
        "receptor_entity_key",
    ]

    return parts


def prepare_candidate_method_events(
    method_name: str,
    *,
    dataset_id: str | None = None,
) -> pd.DataFrame:
    config = METHOD_CONFIG[
        method_name
    ]

    file_path = Path(
        config["path"]
    )

    if not file_path.exists():
        raise FileNotFoundError(
            f"Candidate input not found for {method_name}: "
            f"{file_path}"
        )

    if file_path.stat().st_size == 0:
        raise ValueError(
            f"Candidate input is empty for {method_name}: "
            f"{file_path}"
        )

    df = pd.read_csv(
        file_path,
        low_memory=False,
    )

    required_columns = [
        "sender_cell_type",
        "receiver_cell_type",
        "ligand",
        "receptor",
        "score",
    ]

    missing = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing:
        raise ValueError(
            f"{method_name} missing required columns: "
            + ", ".join(
                missing
            )
        )

    if (
        dataset_id is not None
        and "dataset_id" in df.columns
    ):
        filtered = df[
            df["dataset_id"].astype(str)
            == str(dataset_id)
        ].copy()
        # Each per-run method file contains exactly one dataset's results, so
        # when the stamped dataset_id does not match the requested one (e.g. a
        # placeholder left by standardisation), keep all rows rather than
        # dropping every candidate and failing later with an empty key table.
        df = filtered if len(filtered) else df.copy()

    if "interaction_type" in df.columns:
        interaction_type = (
            df["interaction_type"]
            .astype("string")
            .str.strip()
            .str.lower()
        )

        # Match active consensus semantics when the
        # schema actually provides directional LR labels.
        if interaction_type.notna().any():
            df = df[
                interaction_type
                == "ligand-receptor"
            ].copy()

    for column in [
        "sender_cell_type",
        "receiver_cell_type",
        "ligand",
        "receptor",
    ]:
        df[column] = clean_text(
            df[column]
        )

    df["ligand"] = (
        df["ligand"]
        .str.upper()
    )

    df["receptor"] = (
        df["receptor"]
        .str.upper()
    )

    df["score"] = pd.to_numeric(
        df["score"],
        errors="coerce",
    )

    df = df[
        df["sender_cell_type"].notna()
        &
        df["receiver_cell_type"].notna()
        &
        df["ligand"].notna()
        &
        df["receptor"].notna()
        &
        df["score"].notna()
    ].copy()

    higher_is_better = True

    if (
        "score_higher_is_better"
        in df.columns
        and not df[
            "score_higher_is_better"
        ].dropna().empty
    ):
        higher_is_better = parse_bool(
            df[
                "score_higher_is_better"
            ]
            .dropna()
            .iloc[0],
            default=True,
        )

    key_column = config[
        "interaction_key_column"
    ]

    if key_column is None:
        df["interaction_key"] = (
            create_gene_key(
                df
            )
        )
    else:
        if key_column not in df.columns:
            raise ValueError(
                f"{method_name} is missing candidate key "
                f"column: {key_column}"
            )

        df["interaction_key"] = clean_text(
            df[key_column]
        )

        df = df[
            df[
                "interaction_key"
            ].notna()
        ].copy()

    # Strongest method-specific representation wins,
    # exactly matching active consensus semantics.
    df = df.sort_values(
        "score",
        ascending=not higher_is_better,
    )

    raw_valid_rows = len(
        df
    )

    unique_events = (
        df
        .drop_duplicates(
            "interaction_key",
            keep="first",
        )
        .copy()
        .reset_index(
            drop=True
        )
    )

    key_parts = split_candidate_key(
        unique_events[
            "interaction_key"
        ]
    )

    # Defensive consistency check: sender and receiver
    # encoded in the candidate key must agree with the row.
    sender_mismatch = (
        key_parts[
            "key_sender"
        ].astype(str)
        !=
        unique_events[
            "sender_cell_type"
        ].astype(str)
    )

    receiver_mismatch = (
        key_parts[
            "key_receiver"
        ].astype(str)
        !=
        unique_events[
            "receiver_cell_type"
        ].astype(str)
    )

    if sender_mismatch.any():
        raise ValueError(
            f"{method_name}: sender encoded in "
            "candidate key disagrees with row."
        )

    if receiver_mismatch.any():
        raise ValueError(
            f"{method_name}: receiver encoded in "
            "candidate key disagrees with row."
        )

    unique_events[
        "ligand_entity_key"
    ] = key_parts[
        "ligand_entity_key"
    ].to_numpy()

    unique_events[
        "receptor_entity_key"
    ] = key_parts[
        "receptor_entity_key"
    ].to_numpy()

    unique_events[
        "global_percentile"
    ] = (
        unique_events[
            "score"
        ]
        .rank(
            method="average",
            pct=True,
            ascending=higher_is_better,
        )
        * 100
    )

    unique_events[
        "cell_pair_percentile"
    ] = (
        unique_events
        .groupby(
            [
                "sender_cell_type",
                "receiver_cell_type",
            ]
        )[
            "score"
        ]
        .rank(
            method="average",
            pct=True,
            ascending=higher_is_better,
        )
        * 100
    )

    if (
        "score_name"
        in unique_events.columns
    ):
        score_names = (
            unique_events[
                "score_name"
            ]
            .dropna()
            .astype(str)
        )

        score_name = (
            score_names.iloc[0]
            if not score_names.empty
            else "method-specific score"
        )

    elif method_name == "CellPhoneDB":
        score_name = (
            "mean interaction expression"
        )

    elif method_name == "NATMI":
        score_name = (
            "edge average expression weight"
        )

    elif method_name == "CellChat":
        score_name = (
            "communication probability"
        )

    else:
        score_name = (
            "method-specific score"
        )

    unique_events[
        "method"
    ] = method_name

    unique_events[
        "score_name"
    ] = score_name

    unique_events[
        "score_higher_is_better"
    ] = higher_is_better

    unique_events.attrs[
        "raw_valid_rows"
    ] = raw_valid_rows

    unique_events.attrs[
        "unique_events"
    ] = len(
        unique_events
    )

    unique_events.attrs[
        "key_collapses"
    ] = (
        raw_valid_rows
        -
        len(
            unique_events
        )
    )

    return unique_events[
        [
            "interaction_key",
            "sender_cell_type",
            "receiver_cell_type",
            "ligand",
            "receptor",
            "ligand_entity_key",
            "receptor_entity_key",
            "method",
            "score",
            "score_name",
            "score_higher_is_better",
            "global_percentile",
            "cell_pair_percentile",
        ]
    ].copy()


def build_candidate_consensus(
    selected_methods: list[str],
    *,
    dataset_id: str | None = None,
):
    selected_methods = list(
        dict.fromkeys(
            selected_methods
        )
    )

    if not selected_methods:
        raise ValueError(
            "Select at least one method."
        )

    method_tables = {}

    method_diagnostics = {}

    for method_name in selected_methods:
        table = (
            prepare_candidate_method_events(
                method_name,
                dataset_id=dataset_id,
            )
        )

        method_tables[
            method_name
        ] = table

        method_diagnostics[
            method_name
        ] = {
            "input_file": str(
                METHOD_CONFIG[
                    method_name
                ]["path"]
            ),
            "valid_rows": int(
                table.attrs.get(
                    "raw_valid_rows",
                    len(table),
                )
            ),
            "unique_events": int(
                len(
                    table
                )
            ),
            "within_method_key_collapses": int(
                table.attrs.get(
                    "key_collapses",
                    0,
                )
            ),
        }

    biological_frames = []

    for table in method_tables.values():
        biological_frames.append(
            table[
                [
                    "interaction_key",
                    "sender_cell_type",
                    "receiver_cell_type",
                    "ligand_entity_key",
                    "receptor_entity_key",
                ]
            ]
        )

    consensus = (
        pd.concat(
            biological_frames,
            ignore_index=True,
        )
        .drop_duplicates(
            "interaction_key"
        )
        .reset_index(
            drop=True
        )
    )

    percentile_columns = []
    present_columns = []

    for method_name in selected_methods:
        slug = method_slug(
            method_name
        )

        table = (
            method_tables[
                method_name
            ]
            .copy()
        )

        table[
            f"{slug}_present"
        ] = True

        table = table.rename(
            columns={
                "score":
                    f"{slug}_score",
                "score_name":
                    f"{slug}_score_name",
                "score_higher_is_better":
                    f"{slug}_score_higher_is_better",
                "global_percentile":
                    f"{slug}_global_percentile",
                "cell_pair_percentile":
                    f"{slug}_cell_pair_percentile",
                "ligand":
                    f"{slug}_ligand_label",
                "receptor":
                    f"{slug}_receptor_label",
            }
        )

        keep_columns = [
            "interaction_key",
            f"{slug}_present",
            f"{slug}_score",
            f"{slug}_score_name",
            f"{slug}_score_higher_is_better",
            f"{slug}_global_percentile",
            f"{slug}_cell_pair_percentile",
            f"{slug}_ligand_label",
            f"{slug}_receptor_label",
        ]

        consensus = consensus.merge(
            table[
                keep_columns
            ],
            on="interaction_key",
            how="left",
            validate="one_to_one",
        )

        percentile_columns.append(
            f"{slug}_global_percentile"
        )

        present_columns.append(
            f"{slug}_present"
        )

    for column in present_columns:
        consensus[column] = (
            consensus[column]
            .fillna(False)
            .astype(bool)
        )

    consensus[
        "method_count"
    ] = (
        consensus[
            present_columns
        ]
        .sum(
            axis=1
        )
        .astype(int)
    )

    total_methods = len(
        selected_methods
    )

    consensus[
        "support_fraction"
    ] = (
        consensus[
            "method_count"
        ]
        /
        total_methods
    )

    consensus[
        "support_label"
    ] = (
        consensus[
            "method_count"
        ].astype(str)
        + "/"
        + str(
            total_methods
        )
    )

    percentile_matrix = (
        consensus[
            percentile_columns
        ]
    )

    consensus[
        "mean_available_percentile"
    ] = percentile_matrix.mean(
        axis=1,
        skipna=True,
    )

    consensus[
        "max_method_percentile"
    ] = percentile_matrix.max(
        axis=1,
        skipna=True,
    )

    consensus[
        "min_method_percentile"
    ] = percentile_matrix.min(
        axis=1,
        skipna=True,
    )

    # Same semantics as active engine:
    # a consensus percentile exists only when at least
    # two methods support the event.
    consensus[
        "consensus_percentile"
    ] = np.where(
        consensus[
            "method_count"
        ]
        >= 2,
        consensus[
            "mean_available_percentile"
        ],
        np.nan,
    )

    consensus[
        "rank_disagreement"
    ] = np.where(
        consensus[
            "method_count"
        ]
        >= 2,
        (
            consensus[
                "max_method_percentile"
            ]
            -
            consensus[
                "min_method_percentile"
            ]
        ),
        np.nan,
    )

    consensus[
        "percentile_gap"
    ] = consensus[
        "rank_disagreement"
    ]

    def supporting_methods(row):
        supported = []

        for method_name in selected_methods:
            slug = method_slug(
                method_name
            )

            if bool(
                row[
                    f"{slug}_present"
                ]
            ):
                supported.append(
                    method_name
                )

        return " | ".join(
            supported
        )

    consensus[
        "supporting_methods"
    ] = consensus.apply(
        supporting_methods,
        axis=1,
    )

    # Human-readable display labels are deliberately
    # separate from canonical entity identity.
    #
    # We prefer the first available method label, but
    # matching itself NEVER depends on these labels.
    ligand_label_columns = [
        f"{method_slug(m)}_ligand_label"
        for m in selected_methods
    ]

    receptor_label_columns = [
        f"{method_slug(m)}_receptor_label"
        for m in selected_methods
    ]

    consensus[
        "ligand"
    ] = (
        consensus[
            ligand_label_columns
        ]
        .bfill(
            axis=1
        )
        .iloc[
            :,
            0
        ]
    )

    consensus[
        "receptor"
    ] = (
        consensus[
            receptor_label_columns
        ]
        .bfill(
            axis=1
        )
        .iloc[
            :,
            0
        ]
    )

    # Sort is descriptive only.
    consensus = consensus.sort_values(
        by=[
            "method_count",
            "consensus_percentile",
            "mean_available_percentile",
        ],
        ascending=[
            False,
            False,
            False,
        ],
        na_position="last",
    ).reset_index(
        drop=True
    )

    support_counts = (
        consensus[
            "method_count"
        ]
        .value_counts()
        .sort_index()
        .to_dict()
    )

    exact_patterns = {}

    for _, row in consensus.iterrows():
        pattern = tuple(
            method_name
            for method_name in selected_methods
            if bool(
                row[
                    f"{method_slug(method_name)}_present"
                ]
            )
        )

        label = (
            " + ".join(
                pattern
            )
            if pattern
            else "None"
        )

        exact_patterns[
            label
        ] = (
            exact_patterns.get(
                label,
                0,
            )
            + 1
        )

    diagnostics = {
        "runtime_network_required": False,
        "active_registry_changed": False,
        "active_consensus_changed": False,
        "selected_methods": selected_methods,
        "dataset_id": dataset_id,
        "candidate_union": int(
            len(
                consensus
            )
        ),
        "support_counts": {
            str(k): int(v)
            for k, v in support_counts.items()
        },
        "exact_support_patterns": {
            str(k): int(v)
            for k, v in exact_patterns.items()
        },
        "methods": method_diagnostics,
    }

    return (
        consensus,
        diagnostics,
    )


def validate_candidate(
    consensus: pd.DataFrame,
    diagnostics: dict,
    selected_methods: list[str],
):
    if consensus[
        "interaction_key"
    ].duplicated().any():
        raise AssertionError(
            "Candidate consensus contains duplicate "
            "interaction_key values."
        )

    if (
        consensus[
            "method_count"
        ]
        .lt(1)
        .any()
    ):
        raise AssertionError(
            "Candidate consensus contains method_count < 1."
        )

    if (
        consensus[
            "method_count"
        ]
        .gt(
            len(
                selected_methods
            )
        )
        .any()
    ):
        raise AssertionError(
            "Candidate consensus contains impossible "
            "method_count values."
        )

    single_method = (
        consensus[
            "method_count"
        ]
        == 1
    )

    if (
        consensus.loc[
            single_method,
            "consensus_percentile",
        ]
        .notna()
        .any()
    ):
        raise AssertionError(
            "Single-method events must not receive "
            "consensus_percentile."
        )

    if (
        consensus.loc[
            single_method,
            "rank_disagreement",
        ]
        .notna()
        .any()
    ):
        raise AssertionError(
            "Single-method events must not receive "
            "rank_disagreement."
        )

    multi_method = (
        consensus[
            "method_count"
        ]
        >= 2
    )

    if (
        consensus.loc[
            multi_method,
            "consensus_percentile",
        ]
        .isna()
        .any()
    ):
        raise AssertionError(
            "Multi-method events are missing "
            "consensus_percentile."
        )

    if (
        consensus.loc[
            multi_method,
            "rank_disagreement",
        ]
        .isna()
        .any()
    ):
        raise AssertionError(
            "Multi-method events are missing "
            "rank_disagreement."
        )

    for method_name in selected_methods:
        method_info = (
            diagnostics[
                "methods"
            ][
                method_name
            ]
        )

        if (
            method_info[
                "within_method_key_collapses"
            ]
            < 0
        ):
            raise AssertionError(
                "Negative key collapse count."
            )


def print_report(
    consensus: pd.DataFrame,
    diagnostics: dict,
):
    print(
        "=== COMPLEX-AWARE CANDIDATE CONSENSUS ==="
    )
    print(
        "Runtime network required: False"
    )
    print(
        "Active registry changed: False"
    )
    print(
        "Active consensus changed: False"
    )

    print(
        "\n=== METHOD INPUTS ==="
    )

    for (
        method_name,
        info,
    ) in diagnostics[
        "methods"
    ].items():
        print(
            f"{method_name}: "
            f"valid rows={info['valid_rows']:,}, "
            f"unique events={info['unique_events']:,}, "
            f"within-method key collapses="
            f"{info['within_method_key_collapses']:,}"
        )

    print(
        "\nCandidate union:",
        f"{diagnostics['candidate_union']:,}",
    )

    print(
        "\nCandidate support counts:"
    )

    for key, value in diagnostics[
        "support_counts"
    ].items():
        print(
            f"{key}: {value:,}"
        )

    print(
        "\nCandidate exact support patterns:"
    )

    for key, value in diagnostics[
        "exact_support_patterns"
    ].items():
        print(
            f"{key}: {value:,}"
        )

    print(
        "\n=== PERCENTILE SEMANTICS ==="
    )
    print(
        "Percentiles are recalculated within each "
        "candidate method universe."
    )
    print(
        "Raw method scores are never averaged "
        "across methods."
    )
    print(
        "Missing-method percentiles remain missing; "
        "they are never replaced with zero."
    )
    print(
        "consensus_percentile is the mean available "
        "within-method percentile for events supported "
        "by at least two methods."
    )
    print(
        "rank_disagreement is max minus min available "
        "within-method percentile for events supported "
        "by at least two methods."
    )

    print(
        "\n=== STRUCTURAL VALIDATION ==="
    )
    print(
        "Unique interaction keys:",
        not consensus[
            "interaction_key"
        ].duplicated().any(),
    )
    print(
        "Single-method consensus percentile is NaN:",
        consensus.loc[
            consensus[
                "method_count"
            ]
            == 1,
            "consensus_percentile",
        ].isna().all(),
    )
    print(
        "Single-method rank disagreement is NaN:",
        consensus.loc[
            consensus[
                "method_count"
            ]
            == 1,
            "rank_disagreement",
        ].isna().all(),
    )

    print(
        "\n=== CANDIDATE BUILD COMPLETE ==="
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a complex-aware CCC consensus candidate "
            "without changing the active registry, active UI, "
            "or active consensus."
        )
    )

    parser.add_argument(
        "--dataset-id",
        default=None,
    )

    parser.add_argument(
        "--write",
        action="store_true",
        help=(
            "Write candidate-only outputs to "
            "results/consensus_candidate_complex_aware/. "
            "Without --write, the script is dry-run only."
        ),
    )

    parser.add_argument(
        "--methods",
        default="CellPhoneDB,NATMI,CellChat",
        help=(
            "Comma-separated tools to include "
            f"(available: {', '.join(METHOD_CONFIG)})."
        ),
    )

    args = parser.parse_args()

    selected_methods = [
        m.strip()
        for m in args.methods.split(",")
        if m.strip()
    ]
    unknown = [
        m for m in selected_methods
        if m not in METHOD_CONFIG
    ]
    if unknown:
        raise SystemExit(
            "Unknown tool(s): " + ", ".join(unknown)
        )

    consensus, diagnostics = (
        build_candidate_consensus(
            selected_methods,
            dataset_id=args.dataset_id,
        )
    )

    validate_candidate(
        consensus,
        diagnostics,
        selected_methods,
    )

    print_report(
        consensus,
        diagnostics,
    )

    if not args.write:
        print(
            "\nDry run only: no files written."
        )
        print(
            "Active consensus changed: False"
        )
        return

    CANDIDATE_OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    methods_slug = "__".join(
        method_slug(
            method_name
        )
        for method_name in selected_methods
    )

    consensus_path = (
        CANDIDATE_OUTPUT_DIR
        / f"consensus_{methods_slug}.csv"
    )

    manifest_path = (
        CANDIDATE_OUTPUT_DIR
        / f"manifest_{methods_slug}.json"
    )

    consensus.to_csv(
        consensus_path,
        index=False,
    )

    with manifest_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            diagnostics,
            handle,
            indent=2,
            sort_keys=True,
        )

    print(
        "\nCandidate consensus saved:"
    )
    print(
        consensus_path
    )

    print(
        "\nCandidate manifest saved:"
    )
    print(
        manifest_path
    )

    print(
        "\nIMPORTANT:"
    )
    print(
        "The active method registry was NOT changed."
    )
    print(
        "The active consensus directory was NOT used."
    )
    print(
        "The Streamlit UI was NOT changed."
    )


if __name__ == "__main__":
    main()
