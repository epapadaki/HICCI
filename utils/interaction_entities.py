from __future__ import annotations

from dataclasses import dataclass, asdict
import json
import re
from typing import Iterable, Any

import pandas as pd


_SIMPLE_PREFIXES = (
    "simple:",
    "protein:",
    "gene:",
)

_COMPLEX_PREFIXES = (
    "complex:",
    "heteromer:",
    "multimer:",
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


def _clean_gene_symbol(value: Any) -> str | None:
    value = _clean(value)

    if value is None:
        return None

    return value.upper()


def _normalize_subunits(
    subunits: Iterable[Any] | str | None,
) -> tuple[str, ...]:
    """
    Normalize explicitly supplied complex subunits.

    Important:
    We intentionally do NOT guess subunits from arbitrary display
    strings such as "A_B" because method-specific naming conventions
    differ and underscores can occur in identifiers.

    Accepted string separators when the caller explicitly provides
    a subunit field:
        +  |  ;  ,
    """
    if subunits is None:
        return ()

    values: list[Any]

    if isinstance(
        subunits,
        str,
    ):
        raw = _clean(
            subunits
        )

        if raw is None:
            return ()

        values = re.split(
            r"[+|;,]",
            raw,
        )

    else:
        values = list(
            subunits
        )

    normalized = []

    for value in values:
        gene = _clean_gene_symbol(
            value
        )

        if gene is not None:
            normalized.append(
                gene
            )

    # Complex membership is treated as order-invariant for the
    # canonical key, but duplicates are removed.
    return tuple(
        sorted(
            set(
                normalized
            )
        )
    )


def _partner_prefix_type(
    partner_id: Any,
) -> str | None:
    partner = _clean(
        partner_id
    )

    if partner is None:
        return None

    lowered = partner.lower()

    if lowered.startswith(
        _COMPLEX_PREFIXES
    ):
        return "complex"

    if lowered.startswith(
        _SIMPLE_PREFIXES
    ):
        return "gene"

    return None


def _strip_known_prefix(
    value: str | None,
) -> str | None:
    value = _clean(
        value
    )

    if value is None:
        return None

    lowered = value.lower()

    for prefix in (
        *_SIMPLE_PREFIXES,
        *_COMPLEX_PREFIXES,
    ):
        if lowered.startswith(
            prefix
        ):
            return value[
                len(
                    prefix
                ):
            ].strip() or None

    return value


@dataclass(frozen=True)
class MolecularEntity:
    """
    Method-agnostic representation of one ligand/receptor entity.

    entity_type:
        gene
        complex
        unknown

    `display_name` preserves what should be shown to the user.
    `partner_id` preserves the method/resource-specific identifier.
    `subunits` contains only explicitly resolved gene subunits.
    """

    display_name: str | None
    entity_type: str
    partner_id: str | None
    primary_gene: str | None
    subunits: tuple[str, ...]
    subunits_resolved: bool
    canonical_key: str
    source_method: str | None = None

    @property
    def is_complex(
        self,
    ) -> bool:
        return (
            self.entity_type
            == "complex"
        )

    @property
    def subunit_count(
        self,
    ) -> int:
        return len(
            self.subunits
        )

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
            "is_complex"
        ] = self.is_complex

        result[
            "subunit_count"
        ] = self.subunit_count

        return result


def normalize_molecular_entity(
    *,
    display_name: Any = None,
    gene: Any = None,
    partner_id: Any = None,
    subunits: Iterable[Any] | str | None = None,
    entity_type: str | None = None,
    source_method: str | None = None,
) -> MolecularEntity:
    """
    Build a conservative canonical molecular entity.

    Scientific / interoperability rules
    -------------------------------------
    1. Existing display labels are preserved.
    2. Explicit complex subunits are canonicalized order-invariantly.
    3. A method identifier such as ``complex:...`` marks an entity
       as a complex even if its subunits are not yet resolved.
    4. We never invent complex subunits by splitting arbitrary names.
    5. Unresolved complexes are NOT automatically matched across
       methods merely because their display strings look similar.
    """
    display_name_clean = _clean(
        display_name
    )

    gene_clean = _clean_gene_symbol(
        gene
    )

    partner_clean = _clean(
        partner_id
    )

    normalized_subunits = (
        _normalize_subunits(
            subunits
        )
    )

    explicit_type = _clean(
        entity_type
    )

    if explicit_type is not None:
        explicit_type = (
            explicit_type.lower()
        )

        if explicit_type not in {
            "gene",
            "complex",
            "unknown",
        }:
            raise ValueError(
                "entity_type must be one of: gene, complex, unknown."
            )

    prefix_type = (
        _partner_prefix_type(
            partner_clean
        )
    )

    if explicit_type is not None:
        resolved_type = (
            explicit_type
        )

    elif len(
        normalized_subunits
    ) > 1:
        resolved_type = "complex"

    elif prefix_type is not None:
        resolved_type = (
            prefix_type
        )

    elif gene_clean is not None:
        resolved_type = "gene"

    elif len(
        normalized_subunits
    ) == 1:
        resolved_type = "gene"

    else:
        resolved_type = "unknown"

    if (
        resolved_type
        == "gene"
        and gene_clean is None
        and len(
            normalized_subunits
        ) == 1
    ):
        gene_clean = (
            normalized_subunits[0]
        )

    if display_name_clean is None:
        if resolved_type == "gene":
            display_name_clean = (
                gene_clean
                or _strip_known_prefix(
                    partner_clean
                )
            )

        elif (
            resolved_type
            == "complex"
            and normalized_subunits
        ):
            display_name_clean = (
                " + ".join(
                    normalized_subunits
                )
            )

        else:
            display_name_clean = (
                _strip_known_prefix(
                    partner_clean
                )
            )

    if resolved_type == "gene":
        identifier = (
            gene_clean
            or _clean_gene_symbol(
                _strip_known_prefix(
                    partner_clean
                )
            )
            or _clean_gene_symbol(
                display_name_clean
            )
        )

        canonical_key = (
            "GENE:"
            + identifier
            if identifier
            else "GENE:UNRESOLVED"
        )

        resolved_subunits = (
            (identifier,)
            if identifier
            else ()
        )

        subunits_resolved = bool(
            identifier
        )

    elif resolved_type == "complex":
        if normalized_subunits:
            canonical_key = (
                "COMPLEX:"
                + "+".join(
                    normalized_subunits
                )
            )

            resolved_subunits = (
                normalized_subunits
            )

            subunits_resolved = True

        else:
            # Preserve method/resource identity. We deliberately do
            # not claim cross-method equivalence until subunits are
            # resolved.
            unresolved_id = (
                partner_clean
                or display_name_clean
                or "UNRESOLVED"
            )

            canonical_key = (
                "COMPLEX_UNRESOLVED:"
                + str(
                    unresolved_id
                ).upper()
            )

            resolved_subunits = ()

            subunits_resolved = False

    else:
        unresolved_id = (
            partner_clean
            or display_name_clean
            or "UNRESOLVED"
        )

        canonical_key = (
            "UNKNOWN:"
            + str(
                unresolved_id
            ).upper()
        )

        resolved_subunits = (
            normalized_subunits
        )

        subunits_resolved = bool(
            normalized_subunits
        )

    return MolecularEntity(
        display_name=(
            display_name_clean
        ),
        entity_type=(
            resolved_type
        ),
        partner_id=(
            partner_clean
        ),
        primary_gene=(
            gene_clean
        ),
        subunits=(
            resolved_subunits
        ),
        subunits_resolved=(
            subunits_resolved
        ),
        canonical_key=(
            canonical_key
        ),
        source_method=(
            _clean(
                source_method
            )
        ),
    )


def _entity_columns(
    entity: MolecularEntity,
    *,
    side: str,
) -> dict[str, Any]:
    return {
        f"{side}_entity_type":
            entity.entity_type,

        f"{side}_entity_key":
            entity.canonical_key,

        f"{side}_partner_id":
            entity.partner_id,

        f"{side}_primary_gene":
            entity.primary_gene,

        f"{side}_subunits":
            "|".join(
                entity.subunits
            ),

        f"{side}_subunit_count":
            entity.subunit_count,

        f"{side}_subunits_resolved":
            entity.subunits_resolved,

        f"{side}_is_complex":
            entity.is_complex,
    }


def annotate_interaction_entities(
    df: pd.DataFrame,
    *,
    ligand_column: str = "ligand",
    receptor_column: str = "receptor",
    ligand_gene_column: str | None = None,
    receptor_gene_column: str | None = None,
    ligand_partner_column: str | None = "ligand_partner_id",
    receptor_partner_column: str | None = "receptor_partner_id",
    ligand_subunits_column: str | None = None,
    receptor_subunits_column: str | None = None,
    ligand_entity_type_column: str | None = None,
    receptor_entity_type_column: str | None = None,
    method_column: str | None = "method",
) -> pd.DataFrame:
    """
    Add complex-aware entity metadata while preserving existing
    `ligand` and `receptor` columns unchanged.

    This function does NOT alter consensus calculations.
    """
    result = df.copy()

    def get_value(
        row: pd.Series,
        column: str | None,
    ):
        if (
            column is None
            or column not in row.index
        ):
            return None

        return row[
            column
        ]

    ligand_records = []
    receptor_records = []

    for _, row in result.iterrows():
        method = get_value(
            row,
            method_column,
        )

        ligand_entity = (
            normalize_molecular_entity(
                display_name=get_value(
                    row,
                    ligand_column,
                ),
                gene=get_value(
                    row,
                    ligand_gene_column,
                ),
                partner_id=get_value(
                    row,
                    ligand_partner_column,
                ),
                subunits=get_value(
                    row,
                    ligand_subunits_column,
                ),
                entity_type=get_value(
                    row,
                    ligand_entity_type_column,
                ),
                source_method=method,
            )
        )

        receptor_entity = (
            normalize_molecular_entity(
                display_name=get_value(
                    row,
                    receptor_column,
                ),
                gene=get_value(
                    row,
                    receptor_gene_column,
                ),
                partner_id=get_value(
                    row,
                    receptor_partner_column,
                ),
                subunits=get_value(
                    row,
                    receptor_subunits_column,
                ),
                entity_type=get_value(
                    row,
                    receptor_entity_type_column,
                ),
                source_method=method,
            )
        )

        ligand_records.append(
            _entity_columns(
                ligand_entity,
                side="ligand",
            )
        )

        receptor_records.append(
            _entity_columns(
                receptor_entity,
                side="receptor",
            )
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

    return result


def build_complex_aware_event_key(
    *,
    sender_cell_type: Any,
    receiver_cell_type: Any,
    ligand_entity_key: Any,
    receptor_entity_key: Any,
) -> str:
    """
    Candidate future canonical event key.

    It is provided now for testing/provenance, but existing consensus
    code should continue using its current key until complex-resolution
    coverage has been validated method by method.
    """
    sender = (
        _clean(
            sender_cell_type
        )
        or ""
    )

    receiver = (
        _clean(
            receiver_cell_type
        )
        or ""
    )

    ligand_key = (
        _clean(
            ligand_entity_key
        )
        or "UNKNOWN:UNRESOLVED"
    )

    receptor_key = (
        _clean(
            receptor_entity_key
        )
        or "UNKNOWN:UNRESOLVED"
    )

    return "||".join(
        [
            sender,
            receiver,
            ligand_key,
            receptor_key,
        ]
    )


def serialize_entity(
    entity: MolecularEntity,
) -> str:
    """
    Stable JSON representation useful for provenance/debugging.
    """
    return json.dumps(
        entity.to_dict(),
        sort_keys=True,
        separators=(
            ",",
            ":",
        ),
    )
