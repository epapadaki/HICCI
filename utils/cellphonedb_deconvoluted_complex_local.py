from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any
import hashlib

import pandas as pd

from utils.interaction_entities import (
    normalize_molecular_entity,
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


def _read_table(path: Path) -> pd.DataFrame:
    suffix = (
        path.suffix.lower()
    )

    if suffix in {
        ".txt",
        ".tsv",
    }:
        return pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )

    return pd.read_csv(
        path,
        dtype=str,
        keep_default_na=False,
    )


@dataclass(frozen=True)
class DeconvolutedComplexRecord:
    complex_name: str
    subunits: tuple[str, ...]
    interaction_ids: tuple[str, ...]
    source_path: str
    source_sha256: str
    cpdb_is_complex: bool
    subunit_count: int

    def to_dict(
        self,
    ) -> dict[str, Any]:
        result = asdict(
            self
        )

        result[
            "subunits"
        ] = list(
            self.subunits
        )

        result[
            "interaction_ids"
        ] = list(
            self.interaction_ids
        )

        entity = normalize_molecular_entity(
            display_name=self.complex_name,
            partner_id=(
                "complex:"
                + self.complex_name
            ),
            subunits=self.subunits,
            entity_type="complex",
            source_method="CellPhoneDB",
        )

        result[
            "canonical_entity_key"
        ] = entity.canonical_key

        result[
            "entity"
        ] = entity.to_dict()

        return result


class LocalCellPhoneDBDeconvolutedResolver:
    """
    Resolve CellPhoneDB complexes from a local deconvoluted output.

    This is preferred for the current platform integration because
    it uses the exact CellPhoneDB analysis output that was generated
    for the user's dataset.

    No network access is used.
    """

    def __init__(
        self,
        *,
        source_path: str | Path,
    ):
        self.source_path = Path(
            source_path
        ).resolve()

        if not self.source_path.exists():
            raise FileNotFoundError(
                self.source_path
            )

        self.source_sha256 = _sha256(
            self.source_path
        )

        frame = _read_table(
            self.source_path
        )

        required = {
            "gene_name",
            "complex_name",
        }

        missing = (
            required
            - set(
                frame.columns
            )
        )

        if missing:
            raise ValueError(
                "CellPhoneDB deconvoluted output is missing: "
                + ", ".join(
                    sorted(
                        missing
                    )
                )
            )

        self._explicit_simple_genes: set[str] = set()
        self._explicit_simple_gene_interactions: set[
            tuple[str, str]
        ] = set()

        if "is_complex" in frame.columns:
            is_complex = (
                frame[
                    "is_complex"
                ]
                .astype(str)
                .str.strip()
                .str.lower()
                .isin(
                    {
                        "true",
                        "1",
                        "yes",
                        "y",
                        "t",
                    }
                )
            )

            simple_frame = frame.loc[
                ~is_complex
            ].copy()

            simple_frame[
                "_simple_gene_name"
            ] = (
                simple_frame[
                    "gene_name"
                ]
                .map(
                    _upper
                )
            )

            simple_interaction_column = (
                "id_cp_interaction"
                if "id_cp_interaction"
                in simple_frame.columns
                else None
            )

            for _, simple_row in (
                simple_frame.iterrows()
            ):
                simple_gene = (
                    simple_row[
                        "_simple_gene_name"
                    ]
                )

                if simple_gene is None:
                    continue

                self._explicit_simple_genes.add(
                    simple_gene
                )

                if simple_interaction_column:
                    interaction_id = _clean(
                        simple_row.get(
                            simple_interaction_column
                        )
                    )

                    if interaction_id:
                        self._explicit_simple_gene_interactions.add(
                            (
                                simple_gene,
                                interaction_id,
                            )
                        )

            frame = frame[
                is_complex
            ].copy()

        frame[
            "_complex_name"
        ] = (
            frame[
                "complex_name"
            ]
            .map(
                _clean
            )
        )

        frame[
            "_gene_name"
        ] = (
            frame[
                "gene_name"
            ]
            .map(
                _upper
            )
        )

        frame = frame[
            frame[
                "_complex_name"
            ].notna()
            &
            frame[
                "_gene_name"
            ].notna()
        ].copy()

        interaction_column = (
            "id_cp_interaction"
            if "id_cp_interaction"
            in frame.columns
            else None
        )

        records: dict[
            str,
            DeconvolutedComplexRecord,
        ] = {}

        for complex_name, group in (
            frame.groupby(
                "_complex_name",
                sort=True,
            )
        ):
            subunits = tuple(
                sorted(
                    set(
                        group[
                            "_gene_name"
                        ].dropna()
                    )
                )
            )

            if not subunits:
                continue

            if interaction_column:
                interaction_ids = tuple(
                    sorted(
                        {
                            value
                            for value
                            in (
                                group[
                                    interaction_column
                                ].map(
                                    _clean
                                )
                            )
                            if value
                        }
                    )
                )
            else:
                interaction_ids = ()

            record = (
                DeconvolutedComplexRecord(
                    complex_name=(
                        str(
                            complex_name
                        )
                    ),
                    subunits=(
                        subunits
                    ),
                    interaction_ids=(
                        interaction_ids
                    ),
                    source_path=str(
                        self.source_path
                    ),
                    source_sha256=(
                        self.source_sha256
                    ),
                    cpdb_is_complex=True,
                    subunit_count=len(
                        subunits
                    ),
                )
            )

            records[
                str(
                    complex_name
                ).upper()
            ] = record

        self._records = records


    def is_explicit_simple_gene(
        self,
        gene_name: Any,
        *,
        interaction_id: Any = None,
    ) -> bool:
        """
        Return True when the local CellPhoneDB deconvoluted output
        explicitly marks this gene as is_complex=False.

        When an interaction id is available, the evidence is matched
        within that exact CellPhoneDB interaction. This avoids using a
        global gene-level shortcut when the same gene may also occur as
        a component of a true complex in other interactions.
        """
        gene = _upper(
            gene_name
        )

        if gene is None:
            return False

        interaction = _clean(
            interaction_id
        )

        if interaction is not None:
            return (
                (
                    gene,
                    interaction,
                )
                in self._explicit_simple_gene_interactions
            )

        return (
            gene
            in self._explicit_simple_genes
        )


    def resolve(
        self,
        complex_name_or_partner_id: Any,
    ) -> DeconvolutedComplexRecord | None:
        value = _clean(
            complex_name_or_partner_id
        )

        if value is None:
            return None

        lowered = value.lower()

        if lowered.startswith(
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


    def resolve_entity(
        self,
        complex_name_or_partner_id: Any,
    ) -> dict[str, Any]:
        record = self.resolve(
            complex_name_or_partner_id
        )

        if record is None:
            return {
                "found":
                    False,

                "resolved":
                    False,

                "record":
                    None,

                "reason":
                    (
                        "Complex was not found in the local "
                        "CellPhoneDB deconvoluted output."
                    ),

                "runtime_network_required":
                    False,
            }

        return {
            "found":
                True,

            "resolved":
                True,

            "record":
                record.to_dict(),

            "reason":
                (
                    "Resolved from the local CellPhoneDB "
                    "deconvoluted analysis output."
                ),

            "runtime_network_required":
                False,
        }


    def summary(
        self,
    ) -> dict[str, Any]:
        counts = [
            record.subunit_count
            for record
            in self._records.values()
        ]

        single_subunit = sum(
            value == 1
            for value in counts
        )

        multi_subunit = sum(
            value >= 2
            for value in counts
        )

        return {
            "source_path":
                str(
                    self.source_path
                ),

            "source_sha256":
                self.source_sha256,

            "runtime_network_required":
                False,

            "complex_records":
                len(
                    self._records
                ),

            "explicit_simple_gene_records":
                len(
                    self._explicit_simple_genes
                ),

            "single_subunit_cpdb_complex_records":
                single_subunit,

            "multi_subunit_cpdb_complex_records":
                multi_subunit,
        }


    def to_dataframe(
        self,
    ) -> pd.DataFrame:
        rows = [
            record.to_dict()
            for record
            in self._records.values()
        ]

        return pd.DataFrame(
            rows
        )


def discover_local_deconvoluted_outputs(
    project_root: str | Path,
) -> list[Path]:
    project_root = Path(
        project_root
    ).resolve()

    patterns = [
        "results/cellphonedb/*deconvoluted*.txt",
        "results/cellphonedb/*deconvoluted*.tsv",
        "results/cellphonedb/*deconvoluted*.csv",
        "results/cellphonedb/**/*deconvoluted*.txt",
        "results/cellphonedb/**/*deconvoluted*.tsv",
        "results/cellphonedb/**/*deconvoluted*.csv",
    ]

    paths = []

    for pattern in patterns:
        paths.extend(
            project_root.glob(
                pattern
            )
        )

    unique = {}

    for path in paths:
        if path.is_file():
            unique[
                str(
                    path.resolve()
                )
            ] = path.resolve()

    return sorted(
        unique.values(),
        key=lambda path:
            path.stat().st_mtime,
        reverse=True,
    )
