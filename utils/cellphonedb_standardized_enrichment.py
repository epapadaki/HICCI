from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any
import json

import pandas as pd

from utils.cellphonedb_deconvoluted_complex_local import (
    LocalCellPhoneDBDeconvolutedResolver,
)
from utils.cellphonedb_entity_semantics import (
    LocalCellPhoneDBEntitySemantics,
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


def _bool(value: Any) -> bool:
    if isinstance(
        value,
        bool,
    ):
        return value

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


def _joined(
    values,
) -> str:
    return "|".join(
        sorted(
            {
                str(value).upper()
                for value in values
                if _clean(
                    value
                )
            }
        )
    )


@dataclass(frozen=True)
class EnrichedCellPhoneDBEntity:
    display_name: str | None
    partner_id: str | None

    entity_semantic_type: str
    entity_key: str

    primary_gene: str | None

    subunits: tuple[str, ...]
    subunits_resolved: bool

    proxy_genes: tuple[str, ...]
    associated_genes: tuple[str, ...]

    cpdb_complex_name: str | None
    cpdb_receptor: bool | None
    cpdb_secreted: bool | None
    cpdb_transmembrane: bool | None
    cpdb_integrin: bool | None
    cpdb_other: bool | None
    cpdb_other_desc: str | None

    resolution_source: str
    resolution_note: str

    def to_columns(
        self,
        *,
        side: str,
    ) -> dict[str, Any]:
        return {
            f"{side}_partner_id":
                self.partner_id,

            f"{side}_entity_semantic_type":
                self.entity_semantic_type,

            f"{side}_entity_key":
                self.entity_key,

            f"{side}_primary_gene":
                self.primary_gene,

            f"{side}_subunits":
                _joined(
                    self.subunits
                ),

            f"{side}_subunit_count":
                len(
                    self.subunits
                ),

            f"{side}_subunits_resolved":
                self.subunits_resolved,

            f"{side}_proxy_genes":
                _joined(
                    self.proxy_genes
                ),

            f"{side}_proxy_gene_count":
                len(
                    self.proxy_genes
                ),

            f"{side}_associated_genes":
                _joined(
                    self.associated_genes
                ),

            f"{side}_cpdb_complex_name":
                self.cpdb_complex_name,

            f"{side}_cpdb_receptor":
                self.cpdb_receptor,

            f"{side}_cpdb_secreted":
                self.cpdb_secreted,

            f"{side}_cpdb_transmembrane":
                self.cpdb_transmembrane,

            f"{side}_cpdb_integrin":
                self.cpdb_integrin,

            f"{side}_cpdb_other":
                self.cpdb_other,

            f"{side}_cpdb_other_desc":
                self.cpdb_other_desc,

            f"{side}_resolution_source":
                self.resolution_source,

            f"{side}_resolution_note":
                self.resolution_note,
        }


def _unresolved_complex_key(
    semantic_type: str,
    complex_name: str,
) -> str:
    prefix = {
        "protein_complex":
            "CPDB_PROTEIN_COMPLEX_UNRESOLVED",

        "nonprotein_ligand_proxy":
            "CPDB_PROXY",

        "other_cpdb_complex":
            "CPDB_OTHER",

        "unclassified_cpdb_complex":
            "CPDB_COMPLEX_UNCLASSIFIED",
    }.get(
        semantic_type,
        "CPDB_COMPLEX_UNCLASSIFIED",
    )

    return (
        prefix
        + ":"
        + complex_name.upper()
    )


def resolve_cpdb_entity(
    *,
    display_name: Any,
    partner_id: Any,
    gene: Any,
    interaction_id: Any = None,
    semantics: LocalCellPhoneDBEntitySemantics,
    deconvoluted: LocalCellPhoneDBDeconvolutedResolver,
) -> EnrichedCellPhoneDBEntity:
    """
    Resolve one CellPhoneDB interaction side.

    Critical distinction:
    - protein_complex -> component genes are stored as subunits
    - nonprotein_ligand_proxy -> genes are stored as proxy_genes,
      NEVER as molecular subunits of the ligand
    - other_cpdb_complex -> genes are stored only as associated_genes
    """
    display_clean = _clean(
        display_name
    )

    partner_clean = _clean(
        partner_id
    )

    gene_clean = _upper(
        gene
    )

    # The dataset-specific CellPhoneDB deconvoluted output is the
    # strongest source for whether this partner was actually treated
    # as a simple protein/gene or as a complex in this interaction.
    #
    # This protects against name collisions where complex_input.csv
    # contains a row whose complex_name is identical to a simple gene
    # symbol (e.g. GALR2), while the actual CellPhoneDB output explicitly
    # reports is_complex=False for that interaction.
    if (
        gene_clean is not None
        and
        deconvoluted.is_explicit_simple_gene(
            gene_clean,
            interaction_id=interaction_id,
        )
    ):
        return EnrichedCellPhoneDBEntity(
            display_name=(
                display_clean
                or gene_clean
            ),
            partner_id=(
                partner_clean
            ),
            entity_semantic_type=(
                "gene"
            ),
            entity_key=(
                "GENE:"
                + gene_clean
            ),
            primary_gene=(
                gene_clean
            ),
            subunits=(
                gene_clean,
            ),
            subunits_resolved=True,
            proxy_genes=(),
            associated_genes=(),
            cpdb_complex_name=None,
            cpdb_receptor=None,
            cpdb_secreted=None,
            cpdb_transmembrane=None,
            cpdb_integrin=None,
            cpdb_other=None,
            cpdb_other_desc=None,
            resolution_source=(
                "cellphonedb_deconvoluted_is_complex_false"
            ),
            resolution_note=(
                "The local CellPhoneDB deconvoluted output explicitly "
                "marks this gene as is_complex=False for the same "
                "interaction, so it is preserved as a simple gene even "
                "if complex_input.csv contains a colliding complex name."
            ),
        )

    semantic_record = None

    # Prefer the original CellPhoneDB partner identifier for complex
    # lookup. Fall back to the display label only if needed.
    for candidate in (
        partner_clean,
        display_clean,
    ):
        if candidate is None:
            continue

        semantic_record = (
            semantics.classify(
                candidate
            )
        )

        if semantic_record is not None:
            break

    if semantic_record is None:
        if gene_clean is not None:
            return EnrichedCellPhoneDBEntity(
                display_name=(
                    display_clean
                    or gene_clean
                ),
                partner_id=(
                    partner_clean
                ),
                entity_semantic_type=(
                    "gene"
                ),
                entity_key=(
                    "GENE:"
                    + gene_clean
                ),
                primary_gene=(
                    gene_clean
                ),
                subunits=(
                    gene_clean,
                ),
                subunits_resolved=True,
                proxy_genes=(),
                associated_genes=(),
                cpdb_complex_name=None,
                cpdb_receptor=None,
                cpdb_secreted=None,
                cpdb_transmembrane=None,
                cpdb_integrin=None,
                cpdb_other=None,
                cpdb_other_desc=None,
                resolution_source=(
                    "cellphonedb_standardized_gene"
                ),
                resolution_note=(
                    "Resolved as a simple gene from the "
                    "CellPhoneDB standardized interaction."
                ),
            )

        unresolved_name = (
            partner_clean
            or display_clean
            or "UNRESOLVED"
        )

        return EnrichedCellPhoneDBEntity(
            display_name=(
                display_clean
            ),
            partner_id=(
                partner_clean
            ),
            entity_semantic_type=(
                "unresolved_entity"
            ),
            entity_key=(
                "CPDB_UNRESOLVED:"
                + unresolved_name.upper()
            ),
            primary_gene=None,
            subunits=(),
            subunits_resolved=False,
            proxy_genes=(),
            associated_genes=(),
            cpdb_complex_name=None,
            cpdb_receptor=None,
            cpdb_secreted=None,
            cpdb_transmembrane=None,
            cpdb_integrin=None,
            cpdb_other=None,
            cpdb_other_desc=None,
            resolution_source=(
                "unresolved"
            ),
            resolution_note=(
                "No simple gene symbol or matching local "
                "CellPhoneDB complex metadata was available."
            ),
        )

    complex_name = (
        semantic_record.complex_name
    )

    deconv_record = (
        deconvoluted.resolve(
            complex_name
        )
    )

    deconv_genes = (
        tuple(
            deconv_record.subunits
        )
        if deconv_record is not None
        else ()
    )

    semantic_type = (
        semantic_record.semantic_type
    )

    common = dict(
        display_name=(
            display_clean
            or complex_name
        ),
        partner_id=(
            partner_clean
        ),
        primary_gene=None,
        cpdb_complex_name=(
            complex_name
        ),
        cpdb_receptor=(
            semantic_record.receptor
        ),
        cpdb_secreted=(
            semantic_record.secreted
        ),
        cpdb_transmembrane=(
            semantic_record.transmembrane
        ),
        cpdb_integrin=(
            semantic_record.integrin
        ),
        cpdb_other=(
            semantic_record.other
        ),
        cpdb_other_desc=(
            semantic_record.other_desc
        ),
    )

    if semantic_type == "protein_complex":
        if deconv_genes:
            entity_key = (
                "COMPLEX:"
                + "+".join(
                    sorted(
                        set(
                            deconv_genes
                        )
                    )
                )
            )

            return EnrichedCellPhoneDBEntity(
                **common,
                entity_semantic_type=(
                    "protein_complex"
                ),
                entity_key=(
                    entity_key
                ),
                subunits=(
                    tuple(
                        sorted(
                            set(
                                deconv_genes
                            )
                        )
                    )
                ),
                subunits_resolved=True,
                proxy_genes=(),
                associated_genes=(),
                resolution_source=(
                    "cellphonedb_local_metadata+deconvoluted"
                ),
                resolution_note=(
                    "CellPhoneDB metadata classifies this as a "
                    "protein complex; component genes were resolved "
                    "from the local deconvoluted analysis output."
                ),
            )

        return EnrichedCellPhoneDBEntity(
            **common,
            entity_semantic_type=(
                "protein_complex"
            ),
            entity_key=(
                _unresolved_complex_key(
                    semantic_type,
                    complex_name,
                )
            ),
            subunits=(),
            subunits_resolved=False,
            proxy_genes=(),
            associated_genes=(),
            resolution_source=(
                "cellphonedb_local_metadata"
            ),
            resolution_note=(
                "CellPhoneDB metadata classifies this as a "
                "protein complex, but component genes were not "
                "resolved from the local deconvoluted output."
            ),
        )

    if semantic_type == "nonprotein_ligand_proxy":
        return EnrichedCellPhoneDBEntity(
            **common,
            entity_semantic_type=(
                "nonprotein_ligand_proxy"
            ),
            entity_key=(
                "CPDB_PROXY:"
                + complex_name.upper()
            ),
            subunits=(),
            subunits_resolved=False,
            proxy_genes=(
                tuple(
                    sorted(
                        set(
                            deconv_genes
                        )
                    )
                )
            ),
            associated_genes=(),
            resolution_source=(
                "cellphonedb_local_metadata+deconvoluted"
                if deconv_genes
                else
                "cellphonedb_local_metadata"
            ),
            resolution_note=(
                "CellPhoneDB explicitly marks this entity as "
                "other_desc=biosynthesis_enzyme. Any listed genes "
                "are expression/biosynthesis proxy genes, not "
                "molecular subunits of the non-protein ligand."
            ),
        )

    if semantic_type == "other_cpdb_complex":
        return EnrichedCellPhoneDBEntity(
            **common,
            entity_semantic_type=(
                "other_cpdb_complex"
            ),
            entity_key=(
                "CPDB_OTHER:"
                + complex_name.upper()
            ),
            subunits=(),
            subunits_resolved=False,
            proxy_genes=(),
            associated_genes=(
                tuple(
                    sorted(
                        set(
                            deconv_genes
                        )
                    )
                )
            ),
            resolution_source=(
                "cellphonedb_local_metadata+deconvoluted"
                if deconv_genes
                else
                "cellphonedb_local_metadata"
            ),
            resolution_note=(
                "CellPhoneDB marks this as an 'other' complex "
                "class. Genes from deconvolution are stored only "
                "as associated genes, not automatically as "
                "protein-complex subunits."
            ),
        )

    return EnrichedCellPhoneDBEntity(
        **common,
        entity_semantic_type=(
            "unclassified_cpdb_complex"
        ),
        entity_key=(
            "CPDB_COMPLEX_UNCLASSIFIED:"
            + complex_name.upper()
        ),
        subunits=(),
        subunits_resolved=False,
        proxy_genes=(),
        associated_genes=(
            tuple(
                sorted(
                    set(
                        deconv_genes
                    )
                )
            )
        ),
        resolution_source=(
            "cellphonedb_local_metadata"
        ),
        resolution_note=(
            "CellPhoneDB complex metadata was available but "
            "did not support a more specific semantic class."
        ),
    )


def _directional_source_fields(
    row: pd.Series,
) -> dict[str, Any] | None:
    receptor_a = _bool(
        row.get(
            "receptor_a"
        )
    )

    receptor_b = _bool(
        row.get(
            "receptor_b"
        )
    )

    if (
        receptor_b
        and not receptor_a
    ):
        return {
            "ligand_partner_id":
                row.get(
                    "partner_a"
                ),

            "receptor_partner_id":
                row.get(
                    "partner_b"
                ),

            "ligand_gene":
                row.get(
                    "gene_a"
                ),

            "receptor_gene":
                row.get(
                    "gene_b"
                ),
        }

    if (
        receptor_a
        and not receptor_b
    ):
        return {
            "ligand_partner_id":
                row.get(
                    "partner_b"
                ),

            "receptor_partner_id":
                row.get(
                    "partner_a"
                ),

            "ligand_gene":
                row.get(
                    "gene_b"
                ),

            "receptor_gene":
                row.get(
                    "gene_a"
                ),
        }

    return None


def enrich_standardized_cellphonedb(
    standardized_df: pd.DataFrame,
    *,
    complex_input_path: str | Path,
    deconvoluted_path: str | Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Add complex/proxy-aware columns to standardized CellPhoneDB rows.

    Existing columns are preserved exactly.
    No score, percentile, directionality or consensus logic is changed.
    """
    required_columns = {
        "partner_a",
        "partner_b",
        "gene_a",
        "gene_b",
        "receptor_a",
        "receptor_b",
        "ligand",
        "receptor",
    }

    missing = (
        required_columns
        - set(
            standardized_df.columns
        )
    )

    if missing:
        raise ValueError(
            "Standardized CellPhoneDB table is missing: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )

    semantics = (
        LocalCellPhoneDBEntitySemantics(
            complex_input_path=(
                complex_input_path
            )
        )
    )

    deconvoluted = (
        LocalCellPhoneDBDeconvolutedResolver(
            source_path=(
                deconvoluted_path
            )
        )
    )

    result = (
        standardized_df.copy()
    )

    ligand_rows = []
    receptor_rows = []
    complex_aware_keys = []

    for _, row in result.iterrows():
        source_fields = (
            _directional_source_fields(
                row
            )
        )

        if source_fields is None:
            ligand_rows.append(
                {}
            )

            receptor_rows.append(
                {}
            )

            complex_aware_keys.append(
                None
            )

            continue

        ligand_entity = (
            resolve_cpdb_entity(
                display_name=row.get(
                    "ligand"
                ),
                partner_id=source_fields[
                    "ligand_partner_id"
                ],
                gene=source_fields[
                    "ligand_gene"
                ],
                interaction_id=row.get(
                    "id_cp_interaction"
                ),
                semantics=(
                    semantics
                ),
                deconvoluted=(
                    deconvoluted
                ),
            )
        )

        receptor_entity = (
            resolve_cpdb_entity(
                display_name=row.get(
                    "receptor"
                ),
                partner_id=source_fields[
                    "receptor_partner_id"
                ],
                gene=source_fields[
                    "receptor_gene"
                ],
                interaction_id=row.get(
                    "id_cp_interaction"
                ),
                semantics=(
                    semantics
                ),
                deconvoluted=(
                    deconvoluted
                ),
            )
        )

        ligand_rows.append(
            ligand_entity.to_columns(
                side="ligand"
            )
        )

        receptor_rows.append(
            receptor_entity.to_columns(
                side="receptor"
            )
        )

        sender = _clean(
            row.get(
                "sender_cell_type"
            )
        )

        receiver = _clean(
            row.get(
                "receiver_cell_type"
            )
        )

        if (
            sender
            and receiver
        ):
            complex_aware_keys.append(
                "||".join(
                    [
                        sender,
                        receiver,
                        ligand_entity.entity_key,
                        receptor_entity.entity_key,
                    ]
                )
            )
        else:
            complex_aware_keys.append(
                None
            )

    ligand_df = pd.DataFrame(
        ligand_rows,
        index=result.index,
    )

    receptor_df = pd.DataFrame(
        receptor_rows,
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

    result[
        "complex_aware_event_key_candidate"
    ] = complex_aware_keys

    directional = (
        result[
            "sender_cell_type"
        ].notna()
        &
        result[
            "receiver_cell_type"
        ].notna()
        &
        result[
            "ligand"
        ].notna()
        &
        result[
            "receptor"
        ].notna()
    )

    summary = {
        "rows_total":
            int(
                len(
                    result
                )
            ),

        "rows_directional":
            int(
                directional.sum()
            ),

        "complex_input_path":
            str(
                Path(
                    complex_input_path
                ).resolve()
            ),

        "deconvoluted_path":
            str(
                Path(
                    deconvoluted_path
                ).resolve()
            ),

        "runtime_network_required":
            False,

        "consensus_changed":
            False,

        "existing_scores_changed":
            False,

        "ligand_semantic_counts":
            (
                result.loc[
                    directional,
                    "ligand_entity_semantic_type",
                ]
                .value_counts(
                    dropna=False
                )
                .to_dict()
            ),

        "receptor_semantic_counts":
            (
                result.loc[
                    directional,
                    "receptor_entity_semantic_type",
                ]
                .value_counts(
                    dropna=False
                )
                .to_dict()
            ),

        "protein_complex_rows":
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

        "nonprotein_proxy_rows":
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
    }

    return (
        result,
        summary,
    )


def save_enrichment_summary(
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
