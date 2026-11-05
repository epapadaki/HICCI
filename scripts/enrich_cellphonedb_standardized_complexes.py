from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(
        0,
        str(PROJECT_ROOT),
    )


from utils.cellphonedb_standardized_enrichment import (
    enrich_standardized_cellphonedb,
    save_enrichment_summary,
)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Add local CellPhoneDB complex/proxy semantics and "
            "verified component genes to the standardized output. "
            "The original standardized file is not overwritten."
        )
    )

    parser.add_argument(
        "--standardized",
        default=(
            "results/cellphonedb/"
            "standardized_interactions.csv"
        ),
    )

    parser.add_argument(
        "--complex-input",
        default=(
            "data/cellphonedb/database/"
            "v5.0.0/complex_input.csv"
        ),
    )

    parser.add_argument(
        "--deconvoluted",
        default=(
            "results/cellphonedb/"
            "simple_analysis_deconvoluted_percents_prototype.txt"
        ),
    )

    parser.add_argument(
        "--output",
        default=(
            "results/cellphonedb/"
            "standardized_interactions_complex_aware.csv"
        ),
    )

    parser.add_argument(
        "--summary",
        default=(
            "results/cellphonedb/"
            "standardized_interactions_complex_aware_summary.json"
        ),
    )

    args = parser.parse_args()

    def project_path(
        value: str,
    ) -> Path:
        path = Path(
            value
        )

        if not path.is_absolute():
            path = (
                PROJECT_ROOT
                / path
            )

        return path.resolve()

    standardized_path = project_path(
        args.standardized
    )

    complex_input_path = project_path(
        args.complex_input
    )

    deconvoluted_path = project_path(
        args.deconvoluted
    )

    output_path = project_path(
        args.output
    )

    summary_path = project_path(
        args.summary
    )

    print(
        "Runtime network required: False"
    )

    print(
        "\nLoading standardized CellPhoneDB:"
    )

    print(
        standardized_path
    )

    frame = pd.read_csv(
        standardized_path,
        low_memory=False,
    )

    original_columns = list(
        frame.columns
    )

    enriched, summary = (
        enrich_standardized_cellphonedb(
            frame,
            complex_input_path=(
                complex_input_path
            ),
            deconvoluted_path=(
                deconvoluted_path
            ),
        )
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    enriched.to_csv(
        output_path,
        index=False,
    )

    save_enrichment_summary(
        summary,
        summary_path,
    )

    print(
        "\n=== CELLPHONEDB COMPLEX-AWARE ENRICHMENT ==="
    )

    for key, value in (
        summary.items()
    ):
        print(
            f"{key}: {value}"
        )

    print(
        "\nOriginal columns preserved:",
        all(
            column
            in enriched.columns
            for column
            in original_columns
        ),
    )

    print(
        "Original row count preserved:",
        len(
            frame
        )
        == len(
            enriched
        ),
    )

    print(
        "\nExample protein-complex rows:"
    )

    protein_mask = (
        (
            enriched[
                "ligand_entity_semantic_type"
            ]
            == "protein_complex"
        )
        |
        (
            enriched[
                "receptor_entity_semantic_type"
            ]
            == "protein_complex"
        )
    )

    columns = [
        "sender_cell_type",
        "receiver_cell_type",
        "ligand",
        "receptor",
        "ligand_entity_semantic_type",
        "ligand_subunits",
        "receptor_entity_semantic_type",
        "receptor_subunits",
        "score",
    ]

    print(
        enriched.loc[
            protein_mask,
            columns,
        ]
        .head(
            15
        )
        .to_string(
            index=False
        )
    )

    print(
        "\nExample non-protein proxy rows:"
    )

    proxy_mask = (
        (
            enriched[
                "ligand_entity_semantic_type"
            ]
            == "nonprotein_ligand_proxy"
        )
        |
        (
            enriched[
                "receptor_entity_semantic_type"
            ]
            == "nonprotein_ligand_proxy"
        )
    )

    proxy_columns = [
        "sender_cell_type",
        "receiver_cell_type",
        "ligand",
        "receptor",
        "ligand_entity_semantic_type",
        "ligand_proxy_genes",
        "receptor_entity_semantic_type",
        "receptor_proxy_genes",
        "score",
    ]

    print(
        enriched.loc[
            proxy_mask,
            proxy_columns,
        ]
        .head(
            15
        )
        .to_string(
            index=False
        )
    )

    print(
        "\nSaved enriched output:"
    )

    print(
        output_path
    )

    print(
        "\nSaved summary:"
    )

    print(
        summary_path
    )

    print(
        "\nIMPORTANT:"
    )

    print(
        "The original standardized_interactions.csv was NOT overwritten."
    )

    print(
        "Existing CCC scores were NOT changed."
    )

    print(
        "Consensus logic was NOT changed."
    )


if __name__ == "__main__":
    main()
