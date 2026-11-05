from __future__ import annotations

import csv
from pathlib import Path
import sqlite3
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]

REACTOME_DB = (
    PROJECT_ROOT
    / "data"
    / "biological_context"
    / "reactome"
    / "reactome_local.sqlite"
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
    if not REACTOME_DB.exists():
        raise FileNotFoundError(
            "Local Reactome database was not found: "
            f"{REACTOME_DB}. Import the downloaded Reactome "
            "files with scripts/prepare_reactome_local.py first."
        )

    con = sqlite3.connect(
        REACTOME_DB
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


def get_receptor_pathways(
    receptor: str,
    *,
    species: str = "Homo sapiens",
) -> dict[str, Any]:
    """
    Return local Reactome pathway membership for a receptor.

    Important:
    - Reactome membership is receptor-general prior knowledge.
    - It is not evidence that the selected ligand activates the
      pathway.
    - It is not evidence that the pathway is active in the
      user's scRNA-seq dataset.
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
            "to a reviewed human UniProt accession. If this is "
            "a receptor complex, keep the complex identity and "
            "handle its subunits explicitly rather than "
            "collapsing it to one gene."
        )

    mapping = mappings[0]

    accession = mapping[
        "accession"
    ]

    con = _connect()

    rows = con.execute(
        """
        SELECT DISTINCT
            up.uniprot_id,
            up.reactome_id,
            up.reactome_url,
            up.pathway_name AS mapping_pathway_name,
            up.evidence_code,
            up.species,
            ps.pathway_name AS summary_pathway_name,
            ps.summation,
            pc.pathway_name AS catalog_pathway_name
        FROM uniprot_pathways AS up
        LEFT JOIN pathway_summations AS ps
            ON ps.reactome_id = up.reactome_id
        LEFT JOIN pathway_catalog AS pc
            ON pc.reactome_id = up.reactome_id
        WHERE UPPER(up.uniprot_id) = UPPER(?)
          AND LOWER(up.species) = LOWER(?)
        ORDER BY
            up.pathway_name ASC,
            up.reactome_id ASC
        """,
        (
            accession,
            species,
        ),
    ).fetchall()

    con.close()

    pathways = []

    for row in rows:
        pathway_name = (
            row[
                "summary_pathway_name"
            ]
            or row[
                "catalog_pathway_name"
            ]
            or row[
                "mapping_pathway_name"
            ]
            or ""
        )

        pathways.append(
            {
                "reactome_id":
                    row[
                        "reactome_id"
                    ],

                "pathway_name":
                    pathway_name,

                "description":
                    row[
                        "summation"
                    ]
                    or "",

                "evidence_code":
                    row[
                        "evidence_code"
                    ]
                    or "",

                "species":
                    row[
                        "species"
                    ]
                    or "",

                "reactome_url":
                    row[
                        "reactome_url"
                    ]
                    or (
                        "https://reactome.org/content/detail/"
                        + row[
                            "reactome_id"
                        ]
                    ),

                "source":
                    "Reactome local snapshot",

                "evidence_scope":
                    "receptor_general",

                "effect":
                    (
                        "not specified by pathway membership"
                    ),

                "interpretation_note":
                    (
                        "Reactome maps the receptor protein to "
                        "this curated pathway. Membership does "
                        "not show that the selected ligand "
                        "activates the pathway or that the "
                        "pathway is active in the user's dataset."
                    ),
            }
        )

    return {
        "receptor_input":
            receptor_clean,

        "uniprot_mapping":
            mapping,

        "pathways":
            pathways,

        "metadata":
            load_metadata(),

        "scope":
            "receptor_general",

        "runtime_network_required":
            False,

        "scientific_caution":
            (
                "These Reactome records are receptor-general "
                "curated pathway memberships. They do not "
                "establish ligand-specific activation and do "
                "not demonstrate pathway activity in the "
                "user's dataset."
            ),
    }
