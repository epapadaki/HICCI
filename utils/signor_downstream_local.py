from __future__ import annotations

from collections import deque
from pathlib import Path
import sqlite3
from typing import Any


from utils.signor_local import (
    SIGNOR_DB,
    load_metadata,
    resolve_human_reviewed_uniprot,
)


def _clean(
    value: Any,
) -> str:
    if value is None:
        return ""

    return " ".join(
        str(value).split()
    ).strip()


def _connect() -> sqlite3.Connection:
    if not Path(
        SIGNOR_DB
    ).exists():
        raise FileNotFoundError(
            "Local SIGNOR database was not found: "
            f"{SIGNOR_DB}"
        )

    con = sqlite3.connect(
        SIGNOR_DB
    )

    con.row_factory = sqlite3.Row

    return con


def _fetch_outgoing(
    con: sqlite3.Connection,
    *,
    node_id: str,
    node_name: str,
    human_only: bool = True,
    max_edges: int = 100,
) -> list[sqlite3.Row]:
    node_id_clean = _clean(
        node_id
    )

    node_name_clean = _clean(
        node_name
    )

    if human_only:
        rows = con.execute(
            """
            SELECT *
            FROM relations
            WHERE (
                UPPER(id_a) = UPPER(?)
                OR UPPER(entity_a) = UPPER(?)
            )
            AND (
                tax_id LIKE '%9606%'
                OR tax_id = ''
            )
            ORDER BY
                signor_id ASC
            LIMIT ?
            """,
            (
                node_id_clean,
                node_name_clean,
                int(
                    max_edges
                ),
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
            )
            ORDER BY
                signor_id ASC
            LIMIT ?
            """,
            (
                node_id_clean,
                node_name_clean,
                int(
                    max_edges
                ),
            ),
        ).fetchall()

    return rows


def get_downstream_causal_neighborhood(
    receptor: str,
    *,
    max_depth: int = 2,
    human_only: bool = True,
    max_edges_per_node: int = 100,
) -> dict[str, Any]:
    """
    Traverse the local SIGNOR causal graph downstream from a
    receptor.

    Scientific scope:
    - This returns curated causal paths starting from the
      receptor.
    - It is receptor-general prior knowledge.
    - It does NOT show that the selected ligand triggered these
      paths.
    - It does NOT show pathway activity in the user's dataset.
    """

    receptor_clean = _clean(
        receptor
    )

    if not receptor_clean:
        raise ValueError(
            "A receptor identifier is required."
        )

    if max_depth < 1:
        raise ValueError(
            "max_depth must be at least 1."
        )

    mappings = resolve_human_reviewed_uniprot(
        receptor_clean
    )

    if not mappings:
        raise ValueError(
            f"Could not resolve receptor '{receptor_clean}' "
            "to a reviewed human UniProt accession. "
            "Complex-aware handling should remain separate."
        )

    mapping = mappings[0]

    root_name = (
        mapping[
            "primary_gene"
        ]
        or receptor_clean
    )

    root_id = mapping[
        "accession"
    ]

    con = _connect()

    queue = deque(
        [
            {
                "node_name":
                    root_name,

                "node_id":
                    root_id,

                "depth":
                    0,

                "path_nodes":
                    [
                        root_name
                    ],

                "path_edge_ids":
                    [],
            }
        ]
    )

    results = []
    seen_states = {
        (
            root_id.upper(),
            root_name.upper(),
            0,
        )
    }

    seen_edges = set()

    while queue:
        state = queue.popleft()

        current_depth = int(
            state[
                "depth"
            ]
        )

        if current_depth >= max_depth:
            continue

        outgoing = _fetch_outgoing(
            con,
            node_id=state[
                "node_id"
            ],
            node_name=state[
                "node_name"
            ],
            human_only=human_only,
            max_edges=max_edges_per_node,
        )

        for row in outgoing:
            source_name = _clean(
                row[
                    "entity_a"
                ]
            )

            source_id = _clean(
                row[
                    "id_a"
                ]
            )

            target_name = _clean(
                row[
                    "entity_b"
                ]
            )

            target_id = _clean(
                row[
                    "id_b"
                ]
            )

            signor_id = _clean(
                row[
                    "signor_id"
                ]
            )

            edge_key = (
                signor_id,
                source_id.upper(),
                source_name.upper(),
                target_id.upper(),
                target_name.upper(),
            )

            if edge_key in seen_edges:
                continue

            seen_edges.add(
                edge_key
            )

            edge_depth = (
                current_depth
                + 1
            )

            path_nodes = (
                list(
                    state[
                        "path_nodes"
                    ]
                )
                + [
                    target_name
                    or target_id
                    or "unknown"
                ]
            )

            path_edge_ids = (
                list(
                    state[
                        "path_edge_ids"
                    ]
                )
                + [
                    signor_id
                    or "unlabeled"
                ]
            )

            relation = {
                "depth":
                    edge_depth,

                "source_entity":
                    source_name,

                "source_type":
                    _clean(
                        row[
                            "type_a"
                        ]
                    ),

                "source_id":
                    source_id,

                "source_database":
                    _clean(
                        row[
                            "database_a"
                        ]
                    ),

                "target_entity":
                    target_name,

                "target_type":
                    _clean(
                        row[
                            "type_b"
                        ]
                    ),

                "target_id":
                    target_id,

                "target_database":
                    _clean(
                        row[
                            "database_b"
                        ]
                    ),

                "effect":
                    _clean(
                        row[
                            "effect"
                        ]
                    ),

                "mechanism":
                    _clean(
                        row[
                            "mechanism"
                        ]
                    ),

                "cell_context":
                    _clean(
                        row[
                            "cell_data"
                        ]
                    ),

                "tissue_context":
                    _clean(
                        row[
                            "tissue_data"
                        ]
                    ),

                "pmid":
                    _clean(
                        row[
                            "pmid"
                        ]
                    ),

                "direct":
                    _clean(
                        row[
                            "direct"
                        ]
                    ),

                "signor_id":
                    signor_id,

                "score":
                    _clean(
                        row[
                            "score"
                        ]
                    ),

                "evidence_sentence":
                    _clean(
                        row[
                            "sentence"
                        ]
                    ),

                "path_nodes":
                    path_nodes,

                "path_edge_ids":
                    path_edge_ids,

                "path_text":
                    " -> ".join(
                        path_nodes
                    ),

                "evidence_scope":
                    "receptor_general",

                "source":
                    "SIGNOR local snapshot",

                "interpretation_note":
                    (
                        "This is a curated causal path "
                        "downstream of the receptor in SIGNOR. "
                        "It is not ligand-specific evidence and "
                        "does not demonstrate activity in the "
                        "user's dataset."
                    ),
            }

            results.append(
                relation
            )

            if edge_depth >= max_depth:
                continue

            next_name = (
                target_name
                or target_id
            )

            next_id = target_id

            if not (
                next_name
                or next_id
            ):
                continue

            next_state_key = (
                next_id.upper(),
                next_name.upper(),
                edge_depth,
            )

            if next_state_key in seen_states:
                continue

            seen_states.add(
                next_state_key
            )

            queue.append(
                {
                    "node_name":
                        next_name,

                    "node_id":
                        next_id,

                    "depth":
                        edge_depth,

                    "path_nodes":
                        path_nodes,

                    "path_edge_ids":
                        path_edge_ids,
                }
            )

    con.close()

    results = sorted(
        results,
        key=lambda item: (
            item[
                "depth"
            ],
            item[
                "path_text"
            ],
            item[
                "signor_id"
            ],
        ),
    )

    return {
        "receptor_input":
            receptor_clean,

        "uniprot_mapping":
            mapping,

        "max_depth":
            max_depth,

        "relations":
            results,

        "metadata":
            load_metadata(),

        "runtime_network_required":
            False,

        "scope":
            "receptor_general",

        "scientific_caution":
            (
                "These are curated receptor-general causal "
                "paths from the local SIGNOR snapshot. They do "
                "not establish that the selected ligand "
                "triggered the path and do not demonstrate "
                "pathway activation in the user's dataset."
            ),
    }
