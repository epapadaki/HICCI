from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any
import hashlib
import re

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


def _bool(value: Any) -> bool:
    value = _clean(
        value
    )

    if value is None:
        return False

    return (
        value.lower()
        in {
            "true",
            "1",
            "yes",
            "y",
            "t",
        }
    )


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


def _uniprot_columns(
    frame: pd.DataFrame,
) -> list[str]:
    """
    Support both common CellPhoneDB column spellings:
        uniprot_1, uniprot_2, ...
        uniprot1, uniprot2, ...
    """
    columns = []

    for column in frame.columns:
        match = re.fullmatch(
            r"uniprot_?(\d+)",
            str(
                column
            ),
            flags=re.IGNORECASE,
        )

        if match:
            columns.append(
                (
                    int(
                        match.group(
                            1
                        )
                    ),
                    column,
                )
            )

    return [
        column
        for _, column
        in sorted(
            columns
        )
    ]


@dataclass(frozen=True)
class CellPhoneDBEntitySemantics:
    complex_name: str
    semantic_type: str

    member_uniprots: tuple[str, ...]
    member_count: int

    transmembrane: bool
    peripheral: bool
    secreted: bool
    receptor: bool
    integrin: bool
    other: bool

    other_desc: str | None
    receptor_desc: str | None
    secreted_desc: str | None
    comments_complex: str | None

    classification_basis: str
    source_path: str
    source_sha256: str

    def to_dict(
        self,
    ) -> dict[str, Any]:
        result = asdict(
            self
        )

        result[
            "member_uniprots"
        ] = list(
            self.member_uniprots
        )

        return result


class LocalCellPhoneDBEntitySemantics:
    """
    Classify CellPhoneDB `complex_input.csv` rows using explicit
    CellPhoneDB metadata only.

    No name heuristics are used.

    semantic_type
    -------------
    nonprotein_ligand_proxy
        Explicit CellPhoneDB `other=TRUE` +
        `other_desc=biosynthesis_enzyme`.

    protein_complex
        A CellPhoneDB complex row with UniProt member(s) that is not
        explicitly marked as a non-protein biosynthesis proxy and is
        not marked as another `other` class.

        This includes heteromers, homomers / single-unique-gene
        protein complexes, receptor complexes and secreted complexes.

    other_cpdb_complex
        `other=TRUE`, but not the explicit biosynthesis-enzyme proxy
        category above. We preserve it conservatively rather than
        relabeling it as protein complex.

    unclassified_cpdb_complex
        Insufficient metadata for either category.

    Important
    ---------
    This class describes CellPhoneDB database semantics. It does NOT
    claim that an interaction occurs in the user's biological sample.
    """

    def __init__(
        self,
        *,
        complex_input_path: str | Path,
    ):
        self.complex_input_path = Path(
            complex_input_path
        ).resolve()

        if not self.complex_input_path.exists():
            raise FileNotFoundError(
                self.complex_input_path
            )

        self.source_sha256 = _sha256(
            self.complex_input_path
        )

        frame = pd.read_csv(
            self.complex_input_path,
            dtype=str,
            keep_default_na=False,
        )

        if "complex_name" not in frame.columns:
            raise ValueError(
                "complex_input.csv must contain complex_name."
            )

        self.uniprot_columns = (
            _uniprot_columns(
                frame
            )
        )

        if not self.uniprot_columns:
            raise ValueError(
                "complex_input.csv contains no "
                "uniprot_1/uniprot1 style columns."
            )

        records: dict[
            str,
            CellPhoneDBEntitySemantics,
        ] = {}

        for _, row in frame.iterrows():
            complex_name = _clean(
                row.get(
                    "complex_name"
                )
            )

            if complex_name is None:
                continue

            member_uniprots = []

            for column in self.uniprot_columns:
                value = _clean(
                    row.get(
                        column
                    )
                )

                if value:
                    member_uniprots.append(
                        value.upper()
                    )

            member_uniprots = tuple(
                dict.fromkeys(
                    member_uniprots
                )
            )

            other = _bool(
                row.get(
                    "other"
                )
            )

            other_desc = _clean(
                row.get(
                    "other_desc"
                )
            )

            other_desc_normalized = (
                other_desc.lower()
                if other_desc
                else None
            )

            if (
                other
                and
                other_desc_normalized
                == "biosynthesis_enzyme"
            ):
                semantic_type = (
                    "nonprotein_ligand_proxy"
                )

                basis = (
                    "CellPhoneDB metadata explicitly marks "
                    "other=TRUE and "
                    "other_desc=biosynthesis_enzyme."
                )

            elif other:
                semantic_type = (
                    "other_cpdb_complex"
                )

                basis = (
                    "CellPhoneDB metadata marks other=TRUE, "
                    "but not as biosynthesis_enzyme; kept "
                    "conservatively as an 'other' CellPhoneDB "
                    "complex class."
                )

            elif member_uniprots:
                semantic_type = (
                    "protein_complex"
                )

                basis = (
                    "CellPhoneDB complex row contains UniProt "
                    "member(s) and is not explicitly marked as "
                    "an 'other' entity."
                )

            else:
                semantic_type = (
                    "unclassified_cpdb_complex"
                )

                basis = (
                    "CellPhoneDB metadata was insufficient for "
                    "a more specific classification."
                )

            record = (
                CellPhoneDBEntitySemantics(
                    complex_name=(
                        complex_name
                    ),
                    semantic_type=(
                        semantic_type
                    ),
                    member_uniprots=(
                        member_uniprots
                    ),
                    member_count=len(
                        member_uniprots
                    ),
                    transmembrane=_bool(
                        row.get(
                            "transmembrane"
                        )
                    ),
                    peripheral=_bool(
                        row.get(
                            "peripheral"
                        )
                    ),
                    secreted=_bool(
                        row.get(
                            "secreted"
                        )
                    ),
                    receptor=_bool(
                        row.get(
                            "receptor"
                        )
                    ),
                    integrin=_bool(
                        row.get(
                            "integrin"
                        )
                    ),
                    other=(
                        other
                    ),
                    other_desc=(
                        other_desc
                    ),
                    receptor_desc=_clean(
                        row.get(
                            "receptor_desc"
                        )
                    ),
                    secreted_desc=_clean(
                        row.get(
                            "secreted_desc"
                        )
                    ),
                    comments_complex=_clean(
                        row.get(
                            "comments_complex"
                        )
                    ),
                    classification_basis=(
                        basis
                    ),
                    source_path=str(
                        self.complex_input_path
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


    def classify(
        self,
        complex_name: Any,
    ) -> CellPhoneDBEntitySemantics | None:
        value = _clean(
            complex_name
        )

        if value is None:
            return None

        if value.lower().startswith(
            "complex:"
        ):
            value = value[
                len(
                    "complex:"
                ):
            ].strip()

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
            "source_path":
                str(
                    self.complex_input_path
                ),

            "source_sha256":
                self.source_sha256,

            "runtime_network_required":
                False,

            "total_complex_rows":
                len(
                    self._records
                ),

            "semantic_counts":
                counts,
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
