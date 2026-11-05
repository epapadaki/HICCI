from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any
import hashlib

import pandas as pd


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
    value = _clean(value)

    return (
        value.upper()
        if value is not None
        else None
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


@dataclass(frozen=True)
class CellChatComplexRecord:
    complex_name: str

    # Raw ordered component list exactly as stored by CellChat,
    # excluding empty fields.
    ordered_subunits: tuple[str, ...]

    # Unique normalized gene symbols used only for canonical comparison.
    unique_subunits: tuple[str, ...]

    raw_subunit_count: int
    unique_subunit_count: int
    contains_duplicate_subunits: bool

    source_path: str
    source_sha256: str

    def to_dict(
        self,
    ) -> dict[str, Any]:
        result = asdict(self)

        result["ordered_subunits"] = list(
            self.ordered_subunits
        )

        result["unique_subunits"] = list(
            self.unique_subunits
        )

        result["canonical_entity_key_candidate"] = (
            "COMPLEX:"
            + "+".join(
                self.unique_subunits
            )
            if self.unique_subunits
            else
            "CELLCHAT_COMPLEX_UNRESOLVED:"
            + self.complex_name.upper()
        )

        return result


class LocalCellChatComplexResolver:
    """
    Resolve CellChat complex definitions from a local snapshot exported
    from CellChatDB.human$complex.

    Important scientific constraint:
    This resolver does NOT assume that every CellChat complex definition
    represents a physical multi-protein complex. CellChat also contains
    biosynthetic / neurotransmitter / metabolite proxy constructs.

    Therefore this layer only reports:
      - exact CellChat complex identity
      - exact component fields
      - unique component gene set
      - candidate canonical key

    It does NOT alter consensus and does NOT assign a protein-complex
    semantic label.
    """

    def __init__(
        self,
        *,
        complex_snapshot_path: str | Path,
    ):
        self.complex_snapshot_path = Path(
            complex_snapshot_path
        ).resolve()

        if not self.complex_snapshot_path.exists():
            raise FileNotFoundError(
                self.complex_snapshot_path
            )

        self.source_sha256 = _sha256(
            self.complex_snapshot_path
        )

        frame = pd.read_csv(
            self.complex_snapshot_path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )

        if "complex_name" not in frame.columns:
            raise ValueError(
                "CellChat complex snapshot must contain "
                "'complex_name'."
            )

        subunit_columns = [
            column
            for column in frame.columns
            if str(column).startswith(
                "subunit_"
            )
        ]

        def subunit_index(
            column: str,
        ) -> int:
            try:
                return int(
                    str(column).split(
                        "_",
                        1,
                    )[1]
                )
            except Exception:
                return 999999

        subunit_columns = sorted(
            subunit_columns,
            key=subunit_index,
        )

        if not subunit_columns:
            raise ValueError(
                "CellChat complex snapshot contains no "
                "subunit_1, subunit_2, ... columns."
            )

        records: dict[
            str,
            CellChatComplexRecord,
        ] = {}

        for _, row in frame.iterrows():
            complex_name = _clean(
                row.get(
                    "complex_name"
                )
            )

            if complex_name is None:
                continue

            ordered = []

            for column in subunit_columns:
                value = _upper(
                    row.get(
                        column
                    )
                )

                if value:
                    ordered.append(
                        value
                    )

            ordered_tuple = tuple(
                ordered
            )

            unique_tuple = tuple(
                sorted(
                    set(
                        ordered_tuple
                    )
                )
            )

            record = (
                CellChatComplexRecord(
                    complex_name=(
                        complex_name
                    ),
                    ordered_subunits=(
                        ordered_tuple
                    ),
                    unique_subunits=(
                        unique_tuple
                    ),
                    raw_subunit_count=len(
                        ordered_tuple
                    ),
                    unique_subunit_count=len(
                        unique_tuple
                    ),
                    contains_duplicate_subunits=(
                        len(
                            ordered_tuple
                        )
                        != len(
                            unique_tuple
                        )
                    ),
                    source_path=str(
                        self.complex_snapshot_path
                    ),
                    source_sha256=(
                        self.source_sha256
                    ),
                )
            )

            records[
                complex_name.upper()
            ] = record

        self._records = records


    def resolve(
        self,
        complex_name: Any,
    ) -> CellChatComplexRecord | None:
        value = _clean(
            complex_name
        )

        if value is None:
            return None

        return self._records.get(
            value.upper()
        )


    def resolve_entity(
        self,
        complex_name: Any,
    ) -> dict[str, Any]:
        record = self.resolve(
            complex_name
        )

        if record is None:
            value = _clean(
                complex_name
            )

            return {
                "found":
                    False,

                "resolved":
                    False,

                "semantic_type":
                    "cellchat_complex_unresolved",

                "record":
                    None,

                "entity_key_candidate":
                    (
                        "CELLCHAT_COMPLEX_UNRESOLVED:"
                        + (
                            value.upper()
                            if value
                            else "UNKNOWN"
                        )
                    ),

                "reason":
                    (
                        "Complex was not found in the local "
                        "CellChat complex snapshot."
                    ),

                "runtime_network_required":
                    False,
            }

        return {
            "found":
                True,

            "resolved":
                bool(
                    record.unique_subunits
                ),

            "semantic_type":
                "cellchat_complex_definition",

            "record":
                record.to_dict(),

            "entity_key_candidate":
                (
                    "COMPLEX:"
                    + "+".join(
                        record.unique_subunits
                    )
                    if record.unique_subunits
                    else
                    "CELLCHAT_COMPLEX_UNRESOLVED:"
                    + record.complex_name.upper()
                ),

            "reason":
                (
                    "Resolved directly from local "
                    "CellChatDB.human$complex snapshot. "
                    "No biological semantic class is inferred yet."
                ),

            "runtime_network_required":
                False,
        }


    def summary(
        self,
    ) -> dict[str, Any]:
        raw_counts = [
            record.raw_subunit_count
            for record
            in self._records.values()
        ]

        duplicate_records = sum(
            record.contains_duplicate_subunits
            for record
            in self._records.values()
        )

        return {
            "source_path":
                str(
                    self.complex_snapshot_path
                ),

            "source_sha256":
                self.source_sha256,

            "runtime_network_required":
                False,

            "complex_records":
                len(
                    self._records
                ),

            "single_component_records":
                sum(
                    value == 1
                    for value
                    in raw_counts
                ),

            "multi_component_records":
                sum(
                    value >= 2
                    for value
                    in raw_counts
                ),

            "records_with_duplicate_component_fields":
                int(
                    duplicate_records
                ),

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
