from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(
        0,
        str(PROJECT_ROOT),
    )


from utils.cellchat_raw_complex_enrichment import (
    enrich_cellchat_raw,
    save_summary,
)


RAW_PATH = (
    PROJECT_ROOT
    / "results"
    / "cellchat"
    / "cellchat_interactions_raw.csv"
)

COMPLEX_SNAPSHOT = (
    PROJECT_ROOT
    / "data"
    / "biological_context"
    / "cellchat"
    / "cellchatdb_human_complex_local.tsv"
)

INTERACTION_SNAPSHOT = (
    PROJECT_ROOT
    / "data"
    / "biological_context"
    / "cellchat"
    / "cellchatdb_human_interaction_local.tsv"
)

OUTPUT_PATH = (
    PROJECT_ROOT
    / "results"
    / "cellchat"
    / "cellchat_interactions_complex_aware.csv"
)

SUMMARY_PATH = (
    PROJECT_ROOT
    / "results"
    / "cellchat"
    / "cellchat_interactions_complex_aware_summary.json"
)


def main():
    print(
        "Runtime network required: False"
    )

    for path in [
        RAW_PATH,
        COMPLEX_SNAPSHOT,
        INTERACTION_SNAPSHOT,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                path
            )

    raw = pd.read_csv(
        RAW_PATH,
        low_memory=False,
    )

    original_columns = list(
        raw.columns
    )

    enriched, summary = (
        enrich_cellchat_raw(
            raw,
            complex_snapshot_path=(
                COMPLEX_SNAPSHOT
            ),
            interaction_snapshot_path=(
                INTERACTION_SNAPSHOT
            ),
        )
    )

    OUTPUT_PATH.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    enriched.to_csv(
        OUTPUT_PATH,
        index=False,
    )

    save_summary(
        summary,
        SUMMARY_PATH,
    )

    print(
        "\n=== CELLCHAT COMPLEX-AWARE RAW ENRICHMENT ==="
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
            raw
        )
        == len(
            enriched
        ),
    )

    print(
        "\n=== PROTEIN COMPLEX EXAMPLES ==="
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

    print(
        enriched.loc[
            protein_mask,
            [
                "source",
                "target",
                "ligand",
                "receptor",
                "ligand_entity_semantic_type",
                "ligand_subunits",
                "receptor_entity_semantic_type",
                "receptor_subunits",
                "prob",
                "pval",
                "pathway_name",
            ],
        ]
        .head(
            20
        )
        .to_string(
            index=False
        )
    )

    print(
        "\n=== NON-PROTEIN LIGAND PROXY EXAMPLES ==="
    )

    proxy_mask = (
        enriched[
            "ligand_entity_semantic_type"
        ]
        == "nonprotein_ligand_proxy"
    )

    print(
        enriched.loc[
            proxy_mask,
            [
                "source",
                "target",
                "ligand",
                "receptor",
                "ligand_entity_semantic_type",
                "ligand_proxy_genes",
                "receptor_entity_semantic_type",
                "receptor_subunits",
                "prob",
                "pval",
                "pathway_name",
                "annotation",
            ],
        ]
        .head(
            20
        )
        .to_string(
            index=False
        )
    )

    print(
        "\nSaved enriched CellChat raw output:"
    )

    print(
        OUTPUT_PATH
    )

    print(
        "\nSaved summary:"
    )

    print(
        SUMMARY_PATH
    )

    print(
        "\nIMPORTANT:"
    )

    print(
        "The original cellchat_interactions_raw.csv was NOT overwritten."
    )

    print(
        "CellChat prob/pval values were NOT changed."
    )

    print(
        "The current comparison-ready CellChat universe was NOT changed."
    )

    print(
        "Consensus logic was NOT changed."
    )

    print(
        "complex_aware_event_key_candidate is NOT active in consensus."
    )


if __name__ == "__main__":
    main()
