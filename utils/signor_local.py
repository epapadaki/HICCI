from __future__ import annotations

import csv
from pathlib import Path
import sqlite3
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]

SIGNOR_DB = (
    PROJECT_ROOT
    / "data"
    / "biological_context"
    / "signor"
    / "signor_local.sqlite"
)

UNIPROT_HUMAN_REVIEWED = (
    PROJECT_ROOT
    / "data"
    / "biological_context"
    / "uniprot_human_reviewed.tsv"
)


def _clean(
    value: Any,
) -> str:
    if value is None:
        return ""

    return " ".join(
        str(
            value
        ).split()
    ).strip()


def _connect() -> sqlite3.Connection:
    if not SIGNOR_DB.exists():
        raise FileNotFoundError(
            "Local SIGNOR database was not found: "
            f"{SIGNOR_DB}. Import an official SIGNOR release "
            "with scripts/prepare_signor_local.py first."
        )

    con = sqlite3.connect(
        SIGNOR_DB
    )

    con.row_factory = sqlite3.Row

    return con


def load_metadata() -> dict[str, str]:
    con = _connect()

    rows = con.execute(
        """
        SELECT
            key,
            value
        FROM metadata
        """
    ).fetchall()

    con.close()

    return {
        row[
            "key"
        ]:
            row[
                "value"
            ]
        for row
        in rows
    }


def resolve_human_reviewed_uniprot(
    gene_symbol: str,
) -> list[dict[str, str]]:
    """
    Resolve a human gene symbol or UniProt accession against
    the platform's local reviewed-human UniProt snapshot.
    """

    query = _clean(
        gene_symbol
    ).upper()

    if not query:
        return []

    if not UNIPROT_HUMAN_REVIEWED.exists():
        raise FileNotFoundError(
            "Local reviewed human UniProt snapshot was not found: "
            f"{UNIPROT_HUMAN_REVIEWED}"
        )

    exact_accession = []
    primary_gene = []
    alias_gene = []

    with UNIPROT_HUMAN_REVIEWED.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            accession = _clean(
                row.get(
                    "Entry"
                )
            ).upper()

            gene_names = _clean(
                row.get(
                    "Gene Names"
                )
            )

            tokens = [
                token.upper()
                for token
                in gene_names.split()
                if token.strip()
            ]

            record = {
                "accession":
                    accession,

                "primary_gene":
                    (
                        tokens[0]
                        if tokens
                        else ""
                    ),

                "gene_names":
                    gene_names,

                "protein_names":
                    _clean(
                        row.get(
                            "Protein names"
                        )
                    ),
            }

            if accession == query:
                exact_accession.append(
                    record
                )

            elif tokens:
                if tokens[0] == query:
                    primary_gene.append(
                        record
                    )

                elif query in tokens:
                    alias_gene.append(
                        record
                    )

    ordered = (
        exact_accession
        + primary_gene
        + alias_gene
    )

    deduped = []
    seen = set()

    for record in ordered:
        accession = record[
            "accession"
        ]

        if (
            not accession
            or accession in seen
        ):
            continue

        seen.add(
            accession
        )

        deduped.append(
            record
        )

    return deduped


def _row_to_relation(
    row: sqlite3.Row,
) -> dict[str, Any]:
    return {
        "entity_a":
            row[
                "entity_a"
            ],

        "type_a":
            row[
                "type_a"
            ],

        "id_a":
            row[
                "id_a"
            ],

        "database_a":
            row[
                "database_a"
            ],

        "entity_b":
            row[
                "entity_b"
            ],

        "type_b":
            row[
                "type_b"
            ],

        "id_b":
            row[
                "id_b"
            ],

        "database_b":
            row[
                "database_b"
            ],

        "effect":
            row[
                "effect"
            ],

        "mechanism":
            row[
                "mechanism"
            ],

        "tax_id":
            row[
                "tax_id"
            ],

        "cell_context":
            row[
                "cell_data"
            ],

        "tissue_context":
            row[
                "tissue_data"
            ],

        "pmid":
            row[
                "pmid"
            ],

        "direct":
            row[
                "direct"
            ],

        "notes":
            row[
                "notes"
            ],

        "evidence_sentence":
            row[
                "sentence"
            ],

        "signor_id":
            row[
                "signor_id"
            ],

        "score":
            row[
                "score"
            ],

        "source":
            "SIGNOR local snapshot",

        "evidence_scope":
            "receptor_general",

        "interpretation_note":
            (
                "Curated SIGNOR causal relation involving the "
                "receptor. It does not show that the selected "
                "ligand triggers this relation in the user's "
                "dataset."
            ),
    }


def get_outgoing_relations_for_receptor(
    receptor: str,
    *,
    human_only: bool = True,
) -> dict[str, Any]:
    """
    Return causal relations where the receptor is the source
    entity (ENTITYA -> ENTITYB).

    All data are read locally. No network request is performed.
    """

    receptor_clean = _clean(
        receptor
    )

    if not receptor_clean:
        raise ValueError(
            "A receptor identifier is required."
        )

    mappings = resolve_human_reviewed_uniprot(
        receptor_clean
    )

    if not mappings:
        raise ValueError(
            f"Could not resolve receptor '{receptor_clean}' "
            "to a reviewed human UniProt accession. If the "
            "receptor is a complex, do not collapse it to a "
            "single subunit; complex-aware handling should be "
            "implemented separately."
        )

    mapping = mappings[0]

    accession = mapping[
        "accession"
    ]

    primary_gene = mapping[
        "primary_gene"
    ]

    con = _connect()

    if human_only:
        rows = con.execute(
            """
            SELECT *
            FROM relations
            WHERE (
                UPPER(id_a) = UPPER(?)
                OR UPPER(entity_a) = UPPER(?)
                OR UPPER(entity_a) = UPPER(?)
            )
            AND (
                tax_id LIKE '%9606%'
                OR tax_id = ''
            )
            ORDER BY
                CASE
                    WHEN score GLOB '[0-9]*'
                    THEN CAST(score AS REAL)
                    ELSE NULL
                END DESC,
                signor_id ASC
            """,
            (
                accession,
                receptor_clean,
                primary_gene,
            ),
        ).fetchall()

    else:
        rows = con.execute(
            """
            SELECT *
            FROM relations
            WHERE (
                UPPER(id_a) = UPPER(?)
                OR UPPER(entity_a) = UPPER(?)
                OR UPPER(entity_a) = UPPER(?)
            )
            ORDER BY
                signor_id ASC
            """,
            (
                accession,
                receptor_clean,
                primary_gene,
            ),
        ).fetchall()

    con.close()

    relations = [
        _row_to_relation(
            row
        )
        for row
        in rows
    ]

    return {
        "receptor_input":
            receptor_clean,

        "uniprot_mapping":
            mapping,

        "relations":
            relations,

        "metadata":
            load_metadata(),

        "scope":
            "receptor_general",

        "runtime_network_required":
            False,

        "scientific_caution":
            (
                "These are receptor-general curated causal "
                "relations from the local SIGNOR snapshot. "
                "They do not establish that the selected ligand "
                "caused the downstream relation, and they do not "
                "show signaling activity in the user's dataset."
            ),
    }
