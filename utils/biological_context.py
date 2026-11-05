from functools import lru_cache
from pathlib import Path
import re

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

UNIPROT_FILE = (
    PROJECT_ROOT
    / "data"
    / "biological_context"
    / "uniprot_human_reviewed.tsv"
)


def _clean_text(value):
    if value is None:
        return None

    text = str(value).strip()

    if not text or text.lower() == "nan":
        return None

    # Remove UniProt evidence-code blocks such as:
    # {ECO:0000269|PubMed:12345678}
    text = re.sub(
        r"\s*\{ECO:[^}]+\}",
        "",
        text,
    )

    text = re.sub(
        r"\s+",
        " ",
        text,
    ).strip()

    return text or None


def _strip_prefix(
    value,
    prefix,
):
    text = _clean_text(
        value
    )

    if not text:
        return None

    if text.upper().startswith(
        prefix.upper()
    ):
        text = text[
            len(prefix):
        ].strip()

    return text or None


def _first_sentence(
    value,
):
    text = _clean_text(
        value
    )

    if not text:
        return None

    # Keep the annotation concise for the UI.
    match = re.search(
        r"^(.+?[.!?])(?:\s|$)",
        text,
    )

    if match:
        return match.group(
            1
        ).strip()

    return text


def _short_protein_name(
    value,
):
    text = _clean_text(
        value
    )

    if not text:
        return None

    # Protein names can contain many aliases in parentheses.
    # Keep the first curated name for a concise display.
    if " (" in text:
        text = text.split(
            " (",
            1,
        )[0].strip()

    return text or None


def _normalize_gene_symbol(
    value,
):
    return str(
        value
    ).strip().upper()


def _primary_gene_from_field(
    value,
):
    text = _clean_text(
        value
    )

    if not text:
        return None

    # In the downloaded UniProt TSV the "Gene Names" field
    # lists the primary gene first, followed by synonyms.
    return (
        text.split()[0]
        .strip()
        .upper()
    )


@lru_cache(
    maxsize=1
)
def _load_uniprot_index():
    if not UNIPROT_FILE.exists():
        raise FileNotFoundError(
            "Local UniProt file not found: "
            f"{UNIPROT_FILE}"
        )

    table = pd.read_csv(
        UNIPROT_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    # Your downloaded file currently contains:
    # Entry, Entry Name, Protein names, Gene Names,
    # Function [CC], Pathway, Subcellular location [CC]
    required = {
        "Entry",
        "Protein names",
        "Function [CC]",
        "Pathway",
        "Subcellular location [CC]",
    }

    missing = (
        required
        - set(
            table.columns
        )
    )

    if missing:
        raise ValueError(
            "Local UniProt TSV is missing required "
            f"columns: {sorted(missing)}"
        )

    if "Gene Names (primary)" in table.columns:
        gene_column = (
            "Gene Names (primary)"
        )

    elif "Gene Names" in table.columns:
        gene_column = (
            "Gene Names"
        )

    else:
        raise ValueError(
            "Local UniProt TSV must contain either "
            "'Gene Names' or 'Gene Names (primary)'."
        )

    index = {}

    for _, row in table.iterrows():
        if (
            gene_column
            == "Gene Names (primary)"
        ):
            gene_symbol = (
                _normalize_gene_symbol(
                    row.get(
                        gene_column,
                        "",
                    )
                )
            )

        else:
            gene_symbol = (
                _primary_gene_from_field(
                    row.get(
                        gene_column,
                        "",
                    )
                )
            )

        if not gene_symbol:
            continue

        # Keep the first reviewed human record for a primary
        # gene symbol. This prevents synonyms from silently
        # overriding an existing primary-symbol match.
        if gene_symbol in index:
            continue

        function_text = (
            _strip_prefix(
                row.get(
                    "Function [CC]",
                    "",
                ),
                "FUNCTION:",
            )
        )

        pathway_text = (
            _clean_text(
                row.get(
                    "Pathway",
                    "",
                )
            )
        )

        subcellular_text = (
            _strip_prefix(
                row.get(
                    "Subcellular location [CC]",
                    "",
                ),
                "SUBCELLULAR LOCATION:",
            )
        )

        index[
            gene_symbol
        ] = {
            "gene_symbol":
                gene_symbol,

            "entry":
                _clean_text(
                    row.get(
                        "Entry",
                        "",
                    )
                ),

            "entry_name":
                _clean_text(
                    row.get(
                        "Entry Name",
                        "",
                    )
                ),

            "protein_name":
                _short_protein_name(
                    row.get(
                        "Protein names",
                        "",
                    )
                ),

            "function":
                function_text,

            "function_short":
                _first_sentence(
                    function_text
                ),

            "pathway":
                pathway_text,

            "subcellular_location":
                subcellular_text,

            "source":
                (
                    "UniProtKB/Swiss-Prot "
                    "(local reviewed human snapshot)"
                ),
        }

    return index


def get_gene_annotation(
    gene_symbol,
):
    symbol = (
        _normalize_gene_symbol(
            gene_symbol
        )
    )

    try:
        index = (
            _load_uniprot_index()
        )

    except Exception as error:
        return {
            "gene_symbol":
                symbol,

            "found":
                False,

            "error":
                str(
                    error
                ),

            "source":
                (
                    "UniProtKB/Swiss-Prot "
                    "(local reviewed human snapshot)"
                ),
        }

    annotation = index.get(
        symbol
    )

    if annotation is None:
        return {
            "gene_symbol":
                symbol,

            "found":
                False,

            "error":
                None,

            "source":
                (
                    "UniProtKB/Swiss-Prot "
                    "(local reviewed human snapshot)"
                ),
        }

    result = dict(
        annotation
    )

    result[
        "found"
    ] = True

    result[
        "error"
    ] = None

    return result


def _description_for_gene(
    annotation,
):
    if not annotation.get(
        "found"
    ):
        return (
            "No curated local UniProtKB/Swiss-Prot "
            "annotation was found for this gene."
        )

    protein_name = (
        annotation.get(
            "protein_name"
        )
    )

    function_short = (
        annotation.get(
            "function_short"
        )
    )

    if (
        protein_name
        and function_short
    ):
        return (
            f"{protein_name}. "
            f"{function_short}"
        )

    if protein_name:
        return protein_name

    if function_short:
        return function_short

    return (
        "A reviewed human UniProtKB entry was found, "
        "but no concise function annotation was available."
    )


def _possible_interaction_meaning(
    ligand_symbol,
    receptor_symbol,
    ligand_annotation,
    receptor_annotation,
    sender_cell_type,
    receiver_cell_type,
):
    sender = (
        str(
            sender_cell_type
        ).strip()
        if sender_cell_type
        else "the sender cell"
    )

    receiver = (
        str(
            receiver_cell_type
        ).strip()
        if receiver_cell_type
        else "the receiver cell"
    )

    ligand_function = (
        ligand_annotation.get(
            "function_short"
        )
        if ligand_annotation.get(
            "found"
        )
        else None
    )

    receptor_function = (
        receptor_annotation.get(
            "function_short"
        )
        if receptor_annotation.get(
            "found"
        )
        else None
    )

    if receptor_function:
        text = (
            f"A cautious biological interpretation is that "
            f"{ligand_symbol} from {sender} may signal through "
            f"{receptor_symbol} on {receiver}, potentially "
            f"influencing receptor-associated processes such as: "
            f"{receptor_function}"
        )

        if ligand_function:
            text += (
                " The ligand annotation is also consistent with: "
                f"{ligand_function}"
            )

        return text

    if ligand_function:
        return (
            f"This CCC event proposes that {ligand_symbol} "
            f"from {sender} may engage {receptor_symbol} on "
            f"{receiver}. The local UniProt annotation describes "
            f"{ligand_symbol} as: {ligand_function} "
            "A more specific receptor-mediated interpretation "
            "is not available from the local annotations."
        )

    return (
        f"This CCC event proposes that {ligand_symbol} from "
        f"{sender} may engage {receptor_symbol} on {receiver}. "
        "A more specific curated functional interpretation is "
        "not available in the local UniProt snapshot."
    )


def _parse_biological_context_entity(
    *,
    display_name,
    entity_key,
):
    """
    Resolve molecular identity conservatively.

    Rules
    -----
    GENE:X
        One protein-coding gene candidate X.

    COMPLEX:A+B
        Verified canonical complex identity with required
        subunits A and B.

    Other canonical entity types
        Opaque for local UniProt gene annotation.

    If entity_key is absent, historical behavior is preserved:
    the display name is treated as one gene symbol.

    Display labels are never split on underscores.
    """

    display = (
        str(
            display_name
        ).strip()
        if display_name is not None
        else ""
    )

    key = (
        str(
            entity_key
        ).strip()
        if entity_key is not None
        else ""
    )


    if not key:

        gene = (
            _normalize_gene_symbol(
                display
            )
            if display
            else ""
        )

        return {
            "display_name":
                display,

            "entity_key":
                (
                    "GENE:"
                    + gene
                    if gene
                    else ""
                ),

            "entity_type":
                "gene",

            "required_subunits":
                (
                    [gene]
                    if gene
                    else []
                ),

            "gene_resolvable":
                bool(
                    gene
                ),

            "identity_source":
                "legacy_gene_fallback",
        }


    if key.startswith(
        "GENE:"
    ):

        gene = (
            key[
                len("GENE:"):
            ]
            .strip()
            .upper()
        )

        return {
            "display_name":
                display,

            "entity_key":
                key,

            "entity_type":
                "gene",

            "required_subunits":
                (
                    [gene]
                    if gene
                    else []
                ),

            "gene_resolvable":
                bool(
                    gene
                ),

            "identity_source":
                "canonical_entity_key",
        }


    if key.startswith(
        "COMPLEX:"
    ):

        payload = key[
            len("COMPLEX:"):
        ].strip()

        subunits = []

        for token in payload.split(
            "+"
        ):

            gene = (
                token
                .strip()
                .upper()
            )

            if (
                gene
                and gene not in subunits
            ):
                subunits.append(
                    gene
                )


        return {
            "display_name":
                display,

            "entity_key":
                key,

            "entity_type":
                "complex",

            "required_subunits":
                subunits,

            "gene_resolvable":
                bool(
                    subunits
                ),

            "identity_source":
                "canonical_entity_key",
        }


    entity_type = (
        key.split(
            ":",
            1,
        )[0]
        if ":"
        in key
        else "unknown"
    )


    return {
        "display_name":
            display,

        "entity_key":
            key,

        "entity_type":
            entity_type.lower(),

        "required_subunits":
            [],

        "gene_resolvable":
            False,

        "identity_source":
            "canonical_entity_key_opaque",
    }


def _annotation_for_entity(
    identity,
):
    """
    Return legacy-compatible annotation plus explicit
    per-subunit annotations.

    A complex never receives a fabricated single-protein
    UniProt annotation.
    """

    entity_type = identity.get(
        "entity_type"
    )

    subunits = list(
        identity.get(
            "required_subunits",
            [],
        )
    )


    if (
        entity_type == "gene"
        and len(
            subunits
        )
        == 1
    ):

        annotation = get_gene_annotation(
            subunits[0]
        )

        return {
            "legacy_annotation":
                annotation,

            "subunit_annotations":
                [
                    annotation
                ],

            "annotation_scope":
                "single_gene",
        }


    if entity_type == "complex":

        subunit_annotations = [
            get_gene_annotation(
                gene
            )
            for gene in subunits
        ]

        return {
            "legacy_annotation":
                {
                    "gene_symbol":
                        None,

                    "found":
                        False,

                    "error":
                        None,

                    "source":
                        (
                            "UniProtKB/Swiss-Prot "
                            "(local reviewed human snapshot)"
                        ),

                    "annotation_scope":
                        "complex_subunits",

                    "entity_key":
                        identity.get(
                            "entity_key"
                        ),

                    "required_subunits":
                        subunits,
                },

            "subunit_annotations":
                subunit_annotations,

            "annotation_scope":
                "complex_subunits",
        }


    return {
        "legacy_annotation":
            {
                "gene_symbol":
                    None,

                "found":
                    False,

                "error":
                    None,

                "source":
                    (
                        "UniProtKB/Swiss-Prot "
                        "(local reviewed human snapshot)"
                    ),

                "annotation_scope":
                    "opaque_entity",

                "entity_key":
                    identity.get(
                        "entity_key"
                    ),

                "required_subunits":
                    [],
            },

        "subunit_annotations":
            [],

        "annotation_scope":
            "opaque_entity",
    }


def _description_for_entity(
    *,
    identity,
    annotation_bundle,
):
    entity_type = identity.get(
        "entity_type"
    )

    subunits = list(
        identity.get(
            "required_subunits",
            [],
        )
    )


    if (
        entity_type == "gene"
        and len(
            subunits
        )
        == 1
    ):

        return _description_for_gene(
            annotation_bundle[
                "legacy_annotation"
            ]
        )


    if entity_type == "complex":

        if not subunits:

            return (
                "This receptor/ligand is represented as a "
                "protein complex, but no verified canonical "
                "subunit set is available for local annotation."
            )


        found = [
            annotation
            for annotation
            in annotation_bundle.get(
                "subunit_annotations",
                [],
            )
            if annotation.get(
                "found"
            )
        ]


        return (
            "Protein complex with required canonical subunits "
            + ", ".join(
                subunits
            )
            + ". Curated local UniProtKB/Swiss-Prot context "
            + (
                "is available for "
                + str(
                    len(
                        found
                    )
                )
                + " of "
                + str(
                    len(
                        subunits
                    )
                )
                + " subunits. "
                if subunits
                else ""
            )
            + "Annotations are resolved per subunit; no single "
            "UniProt protein record is assumed to represent the "
            "assembled complex."
        )


    return (
        "This molecular entity cannot be conservatively resolved "
        "to a single protein-coding gene or a verified protein-"
        "complex subunit set from its canonical identity. Local "
        "UniProt gene annotation is therefore not applied to the "
        "display label."
    )


def _entity_function_summaries(
    annotation_bundle,
):
    summaries = []

    for annotation in annotation_bundle.get(
        "subunit_annotations",
        [],
    ):

        if not annotation.get(
            "found"
        ):
            continue

        gene = annotation.get(
            "gene_symbol"
        )

        function_short = annotation.get(
            "function_short"
        )

        if (
            gene
            and function_short
        ):

            summaries.append(
                str(
                    gene
                )
                + ": "
                + str(
                    function_short
                )
            )

    return summaries


def _possible_interaction_meaning_entity(
    *,
    ligand_identity,
    receptor_identity,
    ligand_bundle,
    receptor_bundle,
    sender_cell_type,
    receiver_cell_type,
):
    """
    Preserve historical wording exactly for simple gene-gene
    events. Use cautious entity-aware wording otherwise.
    """

    ligand_type = ligand_identity.get(
        "entity_type"
    )

    receptor_type = receptor_identity.get(
        "entity_type"
    )

    ligand_subunits = list(
        ligand_identity.get(
            "required_subunits",
            [],
        )
    )

    receptor_subunits = list(
        receptor_identity.get(
            "required_subunits",
            [],
        )
    )


    if (
        ligand_type == "gene"
        and receptor_type == "gene"
        and len(
            ligand_subunits
        )
        == 1
        and len(
            receptor_subunits
        )
        == 1
    ):

        return _possible_interaction_meaning(
            ligand_symbol=(
                ligand_subunits[0]
            ),
            receptor_symbol=(
                receptor_subunits[0]
            ),
            ligand_annotation=(
                ligand_bundle[
                    "legacy_annotation"
                ]
            ),
            receptor_annotation=(
                receptor_bundle[
                    "legacy_annotation"
                ]
            ),
            sender_cell_type=(
                sender_cell_type
            ),
            receiver_cell_type=(
                receiver_cell_type
            ),
        )


    sender = (
        str(
            sender_cell_type
        ).strip()
        if sender_cell_type
        else "the sender cell"
    )

    receiver = (
        str(
            receiver_cell_type
        ).strip()
        if receiver_cell_type
        else "the receiver cell"
    )

    ligand_label = (
        ligand_identity.get(
            "display_name"
        )
        or ligand_identity.get(
            "entity_key"
        )
        or "the ligand entity"
    )

    receptor_label = (
        receptor_identity.get(
            "display_name"
        )
        or receptor_identity.get(
            "entity_key"
        )
        or "the receptor entity"
    )


    text = (
        "This CCC event proposes that "
        + str(
            ligand_label
        )
        + " from "
        + sender
        + " may be associated with signaling through "
        + str(
            receptor_label
        )
        + " on "
        + receiver
        + "."
    )


    if receptor_type == "complex":

        text += (
            " The receptor is represented as a protein complex "
            "with required canonical subunits "
            + ", ".join(
                receptor_subunits
            )
            + ". Local UniProt annotations are interpreted "
            "per subunit and do not by themselves establish "
            "binding to, assembly of, or signaling through the "
            "complete receptor complex in this dataset."
        )


    elif receptor_type not in {
        "gene",
        "complex",
    }:

        text += (
            " The receptor entity is opaque for conservative "
            "gene-level UniProt interpretation, so its display "
            "label is not treated as a gene symbol."
        )


    if ligand_type == "complex":

        text += (
            " The ligand is represented as a protein complex "
            "with required canonical subunits "
            + ", ".join(
                ligand_subunits
            )
            + "; local UniProt context is resolved per subunit."
        )


    elif ligand_type not in {
        "gene",
        "complex",
    }:

        text += (
            " The ligand entity is opaque for conservative "
            "gene-level UniProt interpretation."
        )


    receptor_functions = (
        _entity_function_summaries(
            receptor_bundle
        )
    )

    ligand_functions = (
        _entity_function_summaries(
            ligand_bundle
        )
    )


    if receptor_functions:

        text += (
            " Receptor-subunit annotations include: "
            + " | ".join(
                receptor_functions
            )
        )


    if ligand_functions:

        text += (
            " Ligand/subunit annotations include: "
            + " | ".join(
                ligand_functions
            )
        )


    return text


def get_biological_context(
    ligand,
    receptor,
    sender_cell_type=None,
    receiver_cell_type=None,
    ligand_entity_key=None,
    receptor_entity_key=None,
):
    ligand_identity = (
        _parse_biological_context_entity(
            display_name=ligand,
            entity_key=ligand_entity_key,
        )
    )

    receptor_identity = (
        _parse_biological_context_entity(
            display_name=receptor,
            entity_key=receptor_entity_key,
        )
    )


    ligand_bundle = (
        _annotation_for_entity(
            ligand_identity
        )
    )

    receptor_bundle = (
        _annotation_for_entity(
            receptor_identity
        )
    )


    ligand_annotation = (
        ligand_bundle[
            "legacy_annotation"
        ]
    )

    receptor_annotation = (
        receptor_bundle[
            "legacy_annotation"
        ]
    )


    ligand_description = (
        _description_for_entity(
            identity=ligand_identity,
            annotation_bundle=(
                ligand_bundle
            ),
        )
    )

    receptor_description = (
        _description_for_entity(
            identity=receptor_identity,
            annotation_bundle=(
                receptor_bundle
            ),
        )
    )


    possible_role = (
        _possible_interaction_meaning_entity(
            ligand_identity=(
                ligand_identity
            ),
            receptor_identity=(
                receptor_identity
            ),
            ligand_bundle=(
                ligand_bundle
            ),
            receptor_bundle=(
                receptor_bundle
            ),
            sender_cell_type=(
                sender_cell_type
            ),
            receiver_cell_type=(
                receiver_cell_type
            ),
        )
    )


    ligand_type = ligand_identity.get(
        "entity_type"
    )

    receptor_type = receptor_identity.get(
        "entity_type"
    )


    if (
        ligand_type == "gene"
        and receptor_type == "gene"
    ):

        caution = (
            "This section provides curated biological context "
            "for interpretation. It does not prove that the "
            "ligand-receptor interaction physically occurs in "
            "this dataset."
        )

    else:

        caution = (
            "This section provides curated biological context "
            "for interpretation. Gene-level UniProt annotations "
            "for molecular complexes are resolved per verified "
            "subunit and must not be interpreted as evidence that "
            "the assembled complex forms, physically binds the "
            "partner, or signals in this dataset."
        )


    return {
        # ------------------------------------------
        # Historical public fields
        # ------------------------------------------
        "ligand":
            (
                ligand_identity[
                    "required_subunits"
                ][0]
                if (
                    ligand_type == "gene"
                    and ligand_identity[
                        "required_subunits"
                    ]
                )
                else ligand_identity.get(
                    "display_name"
                )
            ),

        "receptor":
            (
                receptor_identity[
                    "required_subunits"
                ][0]
                if (
                    receptor_type == "gene"
                    and receptor_identity[
                        "required_subunits"
                    ]
                )
                else receptor_identity.get(
                    "display_name"
                )
            ),

        "ligand_description":
            ligand_description,

        "receptor_description":
            receptor_description,

        "possible_role":
            possible_role,

        "caution":
            caution,

        "ligand_annotation":
            ligand_annotation,

        "receptor_annotation":
            receptor_annotation,

        "source":
            (
                "UniProtKB/Swiss-Prot "
                "(local reviewed human snapshot)"
            ),

        # ------------------------------------------
        # Entity-aware fields
        # ------------------------------------------
        "ligand_entity_key":
            ligand_identity.get(
                "entity_key"
            ),

        "receptor_entity_key":
            receptor_identity.get(
                "entity_key"
            ),

        "ligand_entity_type":
            ligand_type,

        "receptor_entity_type":
            receptor_type,

        "ligand_identity_source":
            ligand_identity.get(
                "identity_source"
            ),

        "receptor_identity_source":
            receptor_identity.get(
                "identity_source"
            ),

        "ligand_required_subunits":
            list(
                ligand_identity.get(
                    "required_subunits",
                    [],
                )
            ),

        "receptor_required_subunits":
            list(
                receptor_identity.get(
                    "required_subunits",
                    [],
                )
            ),

        "ligand_annotation_scope":
            ligand_bundle.get(
                "annotation_scope"
            ),

        "receptor_annotation_scope":
            receptor_bundle.get(
                "annotation_scope"
            ),

        "ligand_subunit_annotations":
            list(
                ligand_bundle.get(
                    "subunit_annotations",
                    [],
                )
            ),

        "receptor_subunit_annotations":
            list(
                receptor_bundle.get(
                    "subunit_annotations",
                    [],
                )
            ),
    }

