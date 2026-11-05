from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any
import hashlib

import pandas as pd


PROTEIN_SIGNALING_ANNOTATIONS = {
    "Secreted Signaling",
    "ECM-Receptor",
    "Cell-Cell Contact",
}

NONPROTEIN_SIGNALING_ANNOTATION = (
    "Non-protein Signaling"
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open(
        "rb"
    ) as handle:
        for chunk in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(
                chunk
            )

    return digest.hexdigest()


@dataclass(frozen=True)
class CellChatComplexSemanticsRecord:
    complex_name: str
    semantic_type: str

    ligand_occurrences: int
    receptor_occurrences: int

    ligand_annotations: tuple[str, ...]
    receptor_annotations: tuple[str, ...]

    nonprotein_ligand_occurrences: int
    protein_ligand_occurrences: int

    classification_basis: str

    complex_source_path: str
    interaction_source_path: str
    complex_source_sha256: str
    interaction_source_sha256: str

    def to_dict(
        self,
    ) -> dict[str, Any]:
        result = asdict(
            self
        )

        result[
            "ligand_annotations"
        ] = list(
            self.ligand_annotations
        )

        result[
            "receptor_annotations"
        ] = list(
            self.receptor_annotations
        )

        return result


class LocalCellChatComplexSemantics:
    """
    Classify CellChat complex definitions using ONLY local CellChatDB
    metadata and the side on which a complex appears.

    Important:
    ----------
    CellChat's `annotation` is interaction-level metadata.

    Therefore:
    - ligand + Non-protein Signaling
        -> explicit non-protein ligand proxy evidence

    - receptor + Non-protein Signaling
        -> does NOT imply a non-protein receptor.
           The interaction can be non-protein because of its ligand.

    - ligand + Secreted Signaling / ECM-Receptor / Cell-Cell Contact
        -> protein-complex evidence

    - receptor occurrence
        -> protein-receptor-complex evidence unless contradicted by
           explicit ligand-side non-protein use for the same complex.

    Conflict handling:
    ------------------
    If the same complex is explicitly used as a non-protein ligand
    AND also has protein-complex evidence, classify it as
    `mixed_or_ambiguous_cellchat_entity`.

    No name heuristics are used.
    No network access is used.
    """

    def __init__(
        self,
        *,
        complex_snapshot_path: str | Path,
        interaction_snapshot_path: str | Path,
    ):
        self.complex_snapshot_path = Path(
            complex_snapshot_path
        ).resolve()

        self.interaction_snapshot_path = Path(
            interaction_snapshot_path
        ).resolve()

        for path in [
            self.complex_snapshot_path,
            self.interaction_snapshot_path,
        ]:
            if not path.exists():
                raise FileNotFoundError(
                    path
                )

        self.complex_source_sha256 = _sha256(
            self.complex_snapshot_path
        )

        self.interaction_source_sha256 = _sha256(
            self.interaction_snapshot_path
        )

        complex_df = pd.read_csv(
            self.complex_snapshot_path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )

        interaction_df = pd.read_csv(
            self.interaction_snapshot_path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            low_memory=False,
        )

        if (
            "complex_name"
            not in complex_df.columns
        ):
            raise ValueError(
                "CellChat complex snapshot must contain "
                "'complex_name'."
            )

        required_interaction = {
            "ligand",
            "receptor",
            "annotation",
        }

        missing = (
            required_interaction
            - set(
                interaction_df.columns
            )
        )

        if missing:
            raise ValueError(
                "CellChat interaction snapshot is missing: "
                + ", ".join(
                    sorted(
                        missing
                    )
                )
            )

        complex_names = [
            value
            for value
            in (
                complex_df[
                    "complex_name"
                ]
                .map(
                    _clean
                )
                .tolist()
            )
            if value
        ]

        complex_name_set = set(
            complex_names
        )

        usage: dict[
            str,
            dict[str, list[str]],
        ] = {
            name: {
                "ligand":
                    [],
                "receptor":
                    [],
            }
            for name
            in complex_names
        }

        for _, row in (
            interaction_df.iterrows()
        ):
            annotation = (
                _clean(
                    row.get(
                        "annotation"
                    )
                )
                or ""
            )

            ligand = _clean(
                row.get(
                    "ligand"
                )
            )

            receptor = _clean(
                row.get(
                    "receptor"
                )
            )

            if (
                ligand
                in complex_name_set
            ):
                usage[
                    ligand
                ][
                    "ligand"
                ].append(
                    annotation
                )

            if (
                receptor
                in complex_name_set
            ):
                usage[
                    receptor
                ][
                    "receptor"
                ].append(
                    annotation
                )

        records: dict[
            str,
            CellChatComplexSemanticsRecord,
        ] = {}

        for complex_name in (
            complex_names
        ):
            ligand_annotations_raw = (
                usage[
                    complex_name
                ][
                    "ligand"
                ]
            )

            receptor_annotations_raw = (
                usage[
                    complex_name
                ][
                    "receptor"
                ]
            )

            ligand_annotations = tuple(
                sorted(
                    {
                        value
                        for value
                        in ligand_annotations_raw
                        if value
                    }
                )
            )

            receptor_annotations = tuple(
                sorted(
                    {
                        value
                        for value
                        in receptor_annotations_raw
                        if value
                    }
                )
            )

            nonprotein_ligand_occurrences = (
                sum(
                    value
                    == NONPROTEIN_SIGNALING_ANNOTATION
                    for value
                    in ligand_annotations_raw
                )
            )

            protein_ligand_occurrences = (
                sum(
                    value
                    in PROTEIN_SIGNALING_ANNOTATIONS
                    for value
                    in ligand_annotations_raw
                )
            )

            has_nonprotein_ligand_evidence = (
                nonprotein_ligand_occurrences
                > 0
            )

            has_protein_evidence = (
                protein_ligand_occurrences
                > 0
                or
                len(
                    receptor_annotations_raw
                )
                > 0
            )

            if (
                has_nonprotein_ligand_evidence
                and
                has_protein_evidence
            ):
                semantic_type = (
                    "mixed_or_ambiguous_cellchat_entity"
                )

                basis = (
                    "The same CellChat complex has explicit "
                    "ligand-side Non-protein Signaling use and "
                    "also protein-complex evidence from a "
                    "protein-signaling ligand use or receptor use."
                )

            elif has_nonprotein_ligand_evidence:
                semantic_type = (
                    "nonprotein_ligand_proxy"
                )

                basis = (
                    "CellChat explicitly uses this complex on the "
                    "ligand side in interactions annotated "
                    "Non-protein Signaling."
                )

            elif has_protein_evidence:
                semantic_type = (
                    "protein_complex"
                )

                if (
                    len(
                        receptor_annotations_raw
                    )
                    > 0
                    and
                    protein_ligand_occurrences
                    == 0
                ):
                    basis = (
                        "The complex is used on the receptor side; "
                        "interaction-level Non-protein Signaling, "
                        "if present, does not make the receptor "
                        "non-protein."
                    )
                else:
                    basis = (
                        "CellChat uses this complex in protein "
                        "signaling contexts and there is no "
                        "ligand-side Non-protein Signaling evidence."
                    )

            else:
                semantic_type = (
                    "unclassified_cellchat_complex"
                )

                basis = (
                    "The complex definition was not referenced in "
                    "the local interaction snapshot with enough "
                    "side-aware metadata for classification."
                )

            record = (
                CellChatComplexSemanticsRecord(
                    complex_name=(
                        complex_name
                    ),
                    semantic_type=(
                        semantic_type
                    ),
                    ligand_occurrences=len(
                        ligand_annotations_raw
                    ),
                    receptor_occurrences=len(
                        receptor_annotations_raw
                    ),
                    ligand_annotations=(
                        ligand_annotations
                    ),
                    receptor_annotations=(
                        receptor_annotations
                    ),
                    nonprotein_ligand_occurrences=(
                        nonprotein_ligand_occurrences
                    ),
                    protein_ligand_occurrences=(
                        protein_ligand_occurrences
                    ),
                    classification_basis=(
                        basis
                    ),
                    complex_source_path=str(
                        self.complex_snapshot_path
                    ),
                    interaction_source_path=str(
                        self.interaction_snapshot_path
                    ),
                    complex_source_sha256=(
                        self.complex_source_sha256
                    ),
                    interaction_source_sha256=(
                        self.interaction_source_sha256
                    ),
                )
            )

            records[
                complex_name.upper()
            ] = record

        self._records = records


    def classify(
        self,
        complex_name: Any,
    ) -> CellChatComplexSemanticsRecord | None:
        value = _clean(
            complex_name
        )

        if value is None:
            return None

        return self._records.get(
            value.upper()
        )


    def summary(
        self,
    ) -> dict[str, Any]:
        counts = {}

        for record in (
            self._records.values()
        ):
            counts[
                record.semantic_type
            ] = (
                counts.get(
                    record.semantic_type,
                    0,
                )
                + 1
            )

        return {
            "complex_source_path":
                str(
                    self.complex_snapshot_path
                ),

            "interaction_source_path":
                str(
                    self.interaction_snapshot_path
                ),

            "complex_source_sha256":
                self.complex_source_sha256,

            "interaction_source_sha256":
                self.interaction_source_sha256,

            "runtime_network_required":
                False,

            "total_complex_records":
                len(
                    self._records
                ),

            "semantic_counts":
                counts,

            "consensus_changed":
                False,
        }


    def to_dataframe(
        self,
    ) -> pd.DataFrame:
        return pd.DataFrame(
            [
                record.to_dict()
                for record
                in self._records.values()
            ]
        )
