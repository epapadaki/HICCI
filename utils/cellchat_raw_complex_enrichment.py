from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import json

import pandas as pd

from utils.cellchat_complex_local import (
    LocalCellChatComplexResolver,
)
from utils.cellchat_entity_semantics import (
    LocalCellChatComplexSemantics,
)


def _clean(value: Any) -> str | None:
    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except Exception:
        pass

    value = " ".join(
        str(value).split()
    ).strip()

    return value or None


def _upper(value: Any) -> str | None:
    value = _clean(
        value
    )

    return (
        value.upper()
        if value is not None
        else None
    )


def _join(values) -> str:
    return "|".join(
        str(value)
        for value in values
        if _clean(
            value
        )
    )


@dataclass(frozen=True)
class CellChatResolvedEntity:
    display_name: str
    semantic_type: str
    entity_key: str

    primary_gene: str | None

    ordered_subunits: tuple[str, ...]
    subunits: tuple[str, ...]
    subunits_resolved: bool

    proxy_genes: tuple[str, ...]
    associated_genes: tuple[str, ...]

    complex_name: str | None

    ligand_occurrences: int | None
    receptor_occurrences: int | None
    ligand_annotations: tuple[str, ...]
    receptor_annotations: tuple[str, ...]

    resolution_source: str
    resolution_note: str

    def to_columns(
        self,
        *,
        side: str,
    ) -> dict[str, Any]:
        return {
            f"{side}_entity_semantic_type":
                self.semantic_type,

            f"{side}_entity_key":
                self.entity_key,

            f"{side}_primary_gene":
                self.primary_gene,

            f"{side}_ordered_subunits":
                _join(
                    self.ordered_subunits
                ),

            f"{side}_subunits":
                _join(
                    self.subunits
                ),

            f"{side}_subunit_count":
                len(
                    self.subunits
                ),

            f"{side}_subunits_resolved":
                self.subunits_resolved,

            f"{side}_proxy_genes":
                _join(
                    self.proxy_genes
                ),

            f"{side}_proxy_gene_count":
                len(
                    self.proxy_genes
                ),

            f"{side}_associated_genes":
                _join(
                    self.associated_genes
                ),

            f"{side}_cellchat_complex_name":
                self.complex_name,

            f"{side}_cellchat_ligand_occurrences":
                self.ligand_occurrences,

            f"{side}_cellchat_receptor_occurrences":
                self.receptor_occurrences,

            f"{side}_cellchat_ligand_annotations":
                _join(
                    self.ligand_annotations
                ),

            f"{side}_cellchat_receptor_annotations":
                _join(
                    self.receptor_annotations
                ),

            f"{side}_resolution_source":
                self.resolution_source,

            f"{side}_resolution_note":
                self.resolution_note,
        }


def resolve_cellchat_entity(
    *,
    entity_name: Any,
    side: str,
    resolver: LocalCellChatComplexResolver,
    semantics: LocalCellChatComplexSemantics,
) -> CellChatResolvedEntity:
    name = _clean(
        entity_name
    )

    if name is None:
        return CellChatResolvedEntity(
            display_name="",
            semantic_type="unresolved_entity",
            entity_key="CELLCHAT_UNRESOLVED:UNKNOWN",
            primary_gene=None,
            ordered_subunits=(),
            subunits=(),
            subunits_resolved=False,
            proxy_genes=(),
            associated_genes=(),
            complex_name=None,
            ligand_occurrences=None,
            receptor_occurrences=None,
            ligand_annotations=(),
            receptor_annotations=(),
            resolution_source="unresolved",
            resolution_note=(
                "Missing CellChat ligand/receptor label."
            ),
        )

    complex_record = resolver.resolve(
        name
    )

    semantic_record = semantics.classify(
        name
    )

    # Simple gene / protein label.
    if complex_record is None:
        gene = _upper(
            name
        )

        return CellChatResolvedEntity(
            display_name=name,
            semantic_type="gene",
            entity_key=(
                "GENE:"
                + gene
            ),
            primary_gene=gene,
            ordered_subunits=(
                gene,
            ),
            subunits=(
                gene,
            ),
            subunits_resolved=True,
            proxy_genes=(),
            associated_genes=(),
            complex_name=None,
            ligand_occurrences=None,
            receptor_occurrences=None,
            ligand_annotations=(),
            receptor_annotations=(),
            resolution_source=(
                "cellchat_raw_simple_entity"
            ),
            resolution_note=(
                "Entity is not listed in the local "
                "CellChat complex table and is preserved "
                "as a simple gene/protein label."
            ),
        )

    ordered = tuple(
        complex_record.ordered_subunits
    )

    unique = tuple(
        complex_record.unique_subunits
    )

    semantic_type = (
        semantic_record.semantic_type
        if semantic_record is not None
        else "unclassified_cellchat_complex"
    )

    ligand_occurrences = (
        semantic_record.ligand_occurrences
        if semantic_record is not None
        else None
    )

    receptor_occurrences = (
        semantic_record.receptor_occurrences
        if semantic_record is not None
        else None
    )

    ligand_annotations = (
        tuple(
            semantic_record.ligand_annotations
        )
        if semantic_record is not None
        else ()
    )

    receptor_annotations = (
        tuple(
            semantic_record.receptor_annotations
        )
        if semantic_record is not None
        else ()
    )

    common = dict(
        display_name=name,
        primary_gene=None,
        ordered_subunits=ordered,
        complex_name=(
            complex_record.complex_name
        ),
        ligand_occurrences=(
            ligand_occurrences
        ),
        receptor_occurrences=(
            receptor_occurrences
        ),
        ligand_annotations=(
            ligand_annotations
        ),
        receptor_annotations=(
            receptor_annotations
        ),
    )

    if semantic_type == "protein_complex":
        return CellChatResolvedEntity(
            **common,
            semantic_type="protein_complex",
            entity_key=(
                "COMPLEX:"
                + "+".join(
                    unique
                )
                if unique
                else
                "CELLCHAT_PROTEIN_COMPLEX_UNRESOLVED:"
                + name.upper()
            ),
            subunits=unique,
            subunits_resolved=bool(
                unique
            ),
            proxy_genes=(),
            associated_genes=(),
            resolution_source=(
                "cellchat_local_complex+side_aware_semantics"
            ),
            resolution_note=(
                "CellChat complex components were resolved "
                "locally and side-aware CellChat interaction "
                "metadata supports protein-complex semantics."
            ),
        )

    if semantic_type == "nonprotein_ligand_proxy":
        return CellChatResolvedEntity(
            **common,
            semantic_type=(
                "nonprotein_ligand_proxy"
            ),
            entity_key=(
                "CELLCHAT_PROXY:"
                + name.upper()
            ),
            subunits=(),
            subunits_resolved=False,
            proxy_genes=unique,
            associated_genes=(),
            resolution_source=(
                "cellchat_local_complex+side_aware_semantics"
            ),
            resolution_note=(
                "CellChat uses this construct on the ligand "
                "side in Non-protein Signaling. Component "
                "genes are retained as proxy genes, not as "
                "molecular subunits of the non-protein ligand."
            ),
        )

    if (
        semantic_type
        == "mixed_or_ambiguous_cellchat_entity"
    ):
        return CellChatResolvedEntity(
            **common,
            semantic_type=(
                "mixed_or_ambiguous_cellchat_entity"
            ),
            entity_key=(
                "CELLCHAT_AMBIGUOUS:"
                + name.upper()
            ),
            subunits=(),
            subunits_resolved=False,
            proxy_genes=(),
            associated_genes=unique,
            resolution_source=(
                "cellchat_local_complex+side_aware_semantics"
            ),
            resolution_note=(
                "Conflicting CellChat metadata prevents "
                "safe protein-complex or proxy classification. "
                "Components are stored only as associated genes."
            ),
        )

    return CellChatResolvedEntity(
        **common,
        semantic_type=(
            "unclassified_cellchat_complex"
        ),
        entity_key=(
            "CELLCHAT_COMPLEX_UNCLASSIFIED:"
            + name.upper()
        ),
        subunits=(),
        subunits_resolved=False,
        proxy_genes=(),
        associated_genes=unique,
        resolution_source=(
            "cellchat_local_complex"
        ),
        resolution_note=(
            "CellChat complex components are available, "
            "but local side-aware interaction metadata "
            "does not support a safe semantic class."
        ),
    )


def enrich_cellchat_raw(
    raw_df: pd.DataFrame,
    *,
    complex_snapshot_path: str | Path,
    interaction_snapshot_path: str | Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Add complex-aware entity metadata to the full CellChat raw output.

    Existing CellChat columns and row order are preserved.
    `prob` and `pval` are never modified.
    No consensus or cross-method event key is activated here.
    """
    required = {
        "source",
        "target",
        "ligand",
        "receptor",
        "prob",
        "pval",
    }

    missing = (
        required
        - set(
            raw_df.columns
        )
    )

    if missing:
        raise ValueError(
            "CellChat raw output is missing: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )

    resolver = (
        LocalCellChatComplexResolver(
            complex_snapshot_path=(
                complex_snapshot_path
            )
        )
    )

    semantics = (
        LocalCellChatComplexSemantics(
            complex_snapshot_path=(
                complex_snapshot_path
            ),
            interaction_snapshot_path=(
                interaction_snapshot_path
            ),
        )
    )

    result = raw_df.copy()

    ligand_records = []
    receptor_records = []
    candidate_keys = []

    for _, row in (
        result.iterrows()
    ):
        ligand_entity = (
            resolve_cellchat_entity(
                entity_name=row.get(
                    "ligand"
                ),
                side="ligand",
                resolver=resolver,
                semantics=semantics,
            )
        )

        receptor_entity = (
            resolve_cellchat_entity(
                entity_name=row.get(
                    "receptor"
                ),
                side="receptor",
                resolver=resolver,
                semantics=semantics,
            )
        )

        ligand_records.append(
            ligand_entity.to_columns(
                side="ligand"
            )
        )

        receptor_records.append(
            receptor_entity.to_columns(
                side="receptor"
            )
        )

        source = _clean(
            row.get(
                "source"
            )
        )

        target = _clean(
            row.get(
                "target"
            )
        )

        if source and target:
            candidate_keys.append(
                "||".join(
                    [
                        source,
                        target,
                        ligand_entity.entity_key,
                        receptor_entity.entity_key,
                    ]
                )
            )
        else:
            candidate_keys.append(
                None
            )

    ligand_df = pd.DataFrame(
        ligand_records,
        index=result.index,
    )

    receptor_df = pd.DataFrame(
        receptor_records,
        index=result.index,
    )

    for column in ligand_df.columns:
        result[
            column
        ] = ligand_df[
            column
        ]

    for column in receptor_df.columns:
        result[
            column
        ] = receptor_df[
            column
        ]

    # Platform convenience aliases only.
    # Original raw CellChat fields remain unchanged.
    result[
        "method"
    ] = "CellChat"

    result[
        "sender_cell_type"
    ] = result[
        "source"
    ]

    result[
        "receiver_cell_type"
    ] = result[
        "target"
    ]

    result[
        "score"
    ] = result[
        "prob"
    ]

    result[
        "complex_aware_event_key_candidate"
    ] = candidate_keys

    summary = {
        "rows_total":
            int(
                len(
                    result
                )
            ),

        "original_row_count_preserved":
            int(
                len(
                    result
                )
            )
            == int(
                len(
                    raw_df
                )
            ),

        "prob_unchanged":
            bool(
                result[
                    "prob"
                ].equals(
                    raw_df[
                        "prob"
                    ]
                )
            ),

        "pval_unchanged":
            bool(
                result[
                    "pval"
                ].equals(
                    raw_df[
                        "pval"
                    ]
                )
            ),

        "ligand_semantic_counts":
            (
                result[
                    "ligand_entity_semantic_type"
                ]
                .value_counts(
                    dropna=False
                )
                .to_dict()
            ),

        "receptor_semantic_counts":
            (
                result[
                    "receptor_entity_semantic_type"
                ]
                .value_counts(
                    dropna=False
                )
                .to_dict()
            ),

        "rows_with_protein_complex":
            int(
                (
                    (
                        result[
                            "ligand_entity_semantic_type"
                        ]
                        == "protein_complex"
                    )
                    |
                    (
                        result[
                            "receptor_entity_semantic_type"
                        ]
                        == "protein_complex"
                    )
                ).sum()
            ),

        "rows_with_nonprotein_proxy":
            int(
                (
                    (
                        result[
                            "ligand_entity_semantic_type"
                        ]
                        == "nonprotein_ligand_proxy"
                    )
                    |
                    (
                        result[
                            "receptor_entity_semantic_type"
                        ]
                        == "nonprotein_ligand_proxy"
                    )
                ).sum()
            ),

        "rows_with_unclassified_complex":
            int(
                (
                    (
                        result[
                            "ligand_entity_semantic_type"
                        ]
                        == "unclassified_cellchat_complex"
                    )
                    |
                    (
                        result[
                            "receptor_entity_semantic_type"
                        ]
                        == "unclassified_cellchat_complex"
                    )
                ).sum()
            ),

        "rows_with_ambiguous_complex":
            int(
                (
                    (
                        result[
                            "ligand_entity_semantic_type"
                        ]
                        == "mixed_or_ambiguous_cellchat_entity"
                    )
                    |
                    (
                        result[
                            "receptor_entity_semantic_type"
                        ]
                        == "mixed_or_ambiguous_cellchat_entity"
                    )
                ).sum()
            ),

        "runtime_network_required":
            False,

        "existing_cellchat_probabilities_changed":
            False,

        "existing_cellchat_pvalues_changed":
            False,

        "consensus_changed":
            False,

        "comparison_universe_changed":
            False,

        "complex_aware_key_status":
            "candidate_only",
    }

    return (
        result,
        summary,
    )


def save_summary(
    summary: dict[str, Any],
    path: str | Path,
) -> Path:
    path = Path(
        path
    ).resolve()

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    path.write_text(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    return path
