from pathlib import Path
import sqlite3
from datetime import datetime, timezone


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_EVALUATION_DB = (
    PROJECT_ROOT
    / "database"
    / "ccc_evaluation.db"
)


VALID_DECISIONS = {
    "Prioritize",
    "Accept",
    "Reject",
    "Uncertain",
}


def _connect(
    database_file=DEFAULT_EVALUATION_DB,
):
    database_file = Path(
        database_file
    )

    database_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    connection = sqlite3.connect(
        database_file,
        timeout=30,
    )

    connection.row_factory = (
        sqlite3.Row
    )

    connection.execute(
        "PRAGMA foreign_keys = ON"
    )

    return connection


def _utc_now():
    return datetime.now(
        timezone.utc
    ).isoformat(
        timespec="seconds"
    )


def _validate_assessment(
    expert_name,
    decision,
    biological_relevance,
    confidence,
):
    expert_name = str(
        expert_name
    ).strip()

    if not expert_name:
        raise ValueError(
            "Expert identifier is required."
        )

    if decision not in VALID_DECISIONS:
        raise ValueError(
            f"Invalid decision: {decision}"
        )

    biological_relevance = int(
        biological_relevance
    )

    confidence = int(
        confidence
    )

    if biological_relevance not in range(
        1,
        6,
    ):
        raise ValueError(
            "Biological relevance must be 1-5."
        )

    if confidence not in range(
        1,
        6,
    ):
        raise ValueError(
            "Confidence must be 1-5."
        )

    return (
        expert_name,
        biological_relevance,
        confidence,
    )


def init_evaluation_database(
    database_file=DEFAULT_EVALUATION_DB,
):
    connection = _connect(
        database_file
    )

    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS
        blinded_evaluation_responses (
            response_id INTEGER
                PRIMARY KEY AUTOINCREMENT,

            evaluation_set_version TEXT
                NOT NULL,

            evaluation_id TEXT
                NOT NULL,

            expert_name TEXT
                NOT NULL,

            stage INTEGER
                NOT NULL
                CHECK (stage IN (1, 2)),

            dataset_id TEXT
                NOT NULL,

            sender_cell_type TEXT
                NOT NULL,

            receiver_cell_type TEXT
                NOT NULL,

            ligand TEXT
                NOT NULL,

            receptor TEXT
                NOT NULL,

            decision TEXT
                NOT NULL
                CHECK (
                    decision IN (
                        'Prioritize',
                        'Accept',
                        'Reject',
                        'Uncertain'
                    )
                ),

            biological_relevance INTEGER
                NOT NULL
                CHECK (
                    biological_relevance
                    BETWEEN 1 AND 5
                ),

            confidence INTEGER
                NOT NULL
                CHECK (
                    confidence
                    BETWEEN 1 AND 5
                ),

            comment TEXT,

            method_count INTEGER,
            total_selected_methods INTEGER,
            support_fraction REAL,
            method_support TEXT,
            consensus_percentile REAL,
            rank_disagreement REAL,

            created_at TEXT
                NOT NULL,

            UNIQUE (
                evaluation_set_version,
                evaluation_id,
                expert_name,
                stage
            )
        );

        CREATE TABLE IF NOT EXISTS
        blinded_evaluation_method_evidence (
            evidence_id INTEGER
                PRIMARY KEY AUTOINCREMENT,

            response_id INTEGER
                NOT NULL,

            method TEXT
                NOT NULL,

            present INTEGER
                NOT NULL
                CHECK (present IN (0, 1)),

            raw_score REAL,
            score_name TEXT,
            global_percentile REAL,
            cell_pair_percentile REAL,

            FOREIGN KEY (
                response_id
            )
            REFERENCES
                blinded_evaluation_responses(
                    response_id
                )
            ON DELETE CASCADE,

            UNIQUE (
                response_id,
                method
            )
        );

        CREATE INDEX IF NOT EXISTS
        idx_blinded_eval_expert
        ON blinded_evaluation_responses (
            evaluation_set_version,
            expert_name,
            evaluation_id
        );

        CREATE INDEX IF NOT EXISTS
        idx_blinded_eval_stage
        ON blinded_evaluation_responses (
            evaluation_set_version,
            stage
        );
        """
    )

    connection.commit()
    connection.close()


def get_response(
    evaluation_set_version,
    evaluation_id,
    expert_name,
    stage,
    database_file=DEFAULT_EVALUATION_DB,
):
    init_evaluation_database(
        database_file
    )

    connection = _connect(
        database_file
    )

    row = connection.execute(
        """
        SELECT *
        FROM blinded_evaluation_responses
        WHERE
            evaluation_set_version = ?
            AND evaluation_id = ?
            AND expert_name = ?
            AND stage = ?
        LIMIT 1
        """,
        (
            str(
                evaluation_set_version
            ),
            str(
                evaluation_id
            ),
            str(
                expert_name
            ).strip(),
            int(
                stage
            ),
        ),
    ).fetchone()

    connection.close()

    return (
        dict(
            row
        )
        if row is not None
        else None
    )


def load_expert_responses(
    evaluation_set_version,
    expert_name,
    database_file=DEFAULT_EVALUATION_DB,
):
    init_evaluation_database(
        database_file
    )

    connection = _connect(
        database_file
    )

    rows = connection.execute(
        """
        SELECT *
        FROM blinded_evaluation_responses
        WHERE
            evaluation_set_version = ?
            AND expert_name = ?
        ORDER BY
            response_id
        """,
        (
            str(
                evaluation_set_version
            ),
            str(
                expert_name
            ).strip(),
        ),
    ).fetchall()

    connection.close()

    return [
        dict(
            row
        )
        for row in rows
    ]


def load_all_responses(
    evaluation_set_version=None,
    database_file=DEFAULT_EVALUATION_DB,
):
    init_evaluation_database(
        database_file
    )

    connection = _connect(
        database_file
    )

    if evaluation_set_version is None:
        rows = connection.execute(
            """
            SELECT *
            FROM blinded_evaluation_responses
            ORDER BY response_id
            """
        ).fetchall()

    else:
        rows = connection.execute(
            """
            SELECT *
            FROM blinded_evaluation_responses
            WHERE
                evaluation_set_version = ?
            ORDER BY response_id
            """,
            (
                str(
                    evaluation_set_version
                ),
            ),
        ).fetchall()

    connection.close()

    return [
        dict(
            row
        )
        for row in rows
    ]


def load_response_method_evidence(
    response_id,
    database_file=DEFAULT_EVALUATION_DB,
):
    init_evaluation_database(
        database_file
    )

    connection = _connect(
        database_file
    )

    rows = connection.execute(
        """
        SELECT *
        FROM blinded_evaluation_method_evidence
        WHERE response_id = ?
        ORDER BY method
        """,
        (
            int(
                response_id
            ),
        ),
    ).fetchall()

    connection.close()

    return [
        dict(
            row
        )
        for row in rows
    ]


def save_stage1_response(
    *,
    evaluation_set_version,
    evaluation_id,
    expert_name,
    dataset_id,
    sender_cell_type,
    receiver_cell_type,
    ligand,
    receptor,
    decision,
    biological_relevance,
    confidence,
    comment,
    database_file=DEFAULT_EVALUATION_DB,
):
    (
        expert_name,
        biological_relevance,
        confidence,
    ) = _validate_assessment(
        expert_name,
        decision,
        biological_relevance,
        confidence,
    )

    init_evaluation_database(
        database_file
    )

    connection = _connect(
        database_file
    )

    try:
        existing = connection.execute(
            """
            SELECT response_id
            FROM blinded_evaluation_responses
            WHERE
                evaluation_set_version = ?
                AND evaluation_id = ?
                AND expert_name = ?
                AND stage = 1
            """,
            (
                str(
                    evaluation_set_version
                ),
                str(
                    evaluation_id
                ),
                expert_name,
            ),
        ).fetchone()

        if existing is not None:
            raise ValueError(
                "Stage 1 has already been locked "
                "for this expert and event."
            )

        cursor = connection.execute(
            """
            INSERT INTO
            blinded_evaluation_responses (
                evaluation_set_version,
                evaluation_id,
                expert_name,
                stage,
                dataset_id,
                sender_cell_type,
                receiver_cell_type,
                ligand,
                receptor,
                decision,
                biological_relevance,
                confidence,
                comment,
                created_at
            )
            VALUES (
                ?, ?, ?, 1, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
            )
            """,
            (
                str(
                    evaluation_set_version
                ),
                str(
                    evaluation_id
                ),
                expert_name,
                str(
                    dataset_id
                ),
                str(
                    sender_cell_type
                ),
                str(
                    receiver_cell_type
                ),
                str(
                    ligand
                ).upper(),
                str(
                    receptor
                ).upper(),
                str(
                    decision
                ),
                biological_relevance,
                confidence,
                str(
                    comment
                    or ""
                ).strip(),
                _utc_now(),
            ),
        )

        response_id = int(
            cursor.lastrowid
        )

        connection.commit()

        return response_id

    except Exception:
        connection.rollback()
        raise

    finally:
        connection.close()


def save_stage2_response(
    *,
    evaluation_set_version,
    evaluation_id,
    expert_name,
    dataset_id,
    sender_cell_type,
    receiver_cell_type,
    ligand,
    receptor,
    decision,
    biological_relevance,
    confidence,
    comment,
    method_count,
    total_selected_methods,
    support_fraction,
    method_support,
    consensus_percentile,
    rank_disagreement,
    method_evidence,
    database_file=DEFAULT_EVALUATION_DB,
):
    (
        expert_name,
        biological_relevance,
        confidence,
    ) = _validate_assessment(
        expert_name,
        decision,
        biological_relevance,
        confidence,
    )

    init_evaluation_database(
        database_file
    )

    connection = _connect(
        database_file
    )

    try:
        stage1 = connection.execute(
            """
            SELECT response_id
            FROM blinded_evaluation_responses
            WHERE
                evaluation_set_version = ?
                AND evaluation_id = ?
                AND expert_name = ?
                AND stage = 1
            """,
            (
                str(
                    evaluation_set_version
                ),
                str(
                    evaluation_id
                ),
                expert_name,
            ),
        ).fetchone()

        if stage1 is None:
            raise ValueError(
                "Stage 2 cannot be stored before "
                "Stage 1 is locked."
            )

        stage2 = connection.execute(
            """
            SELECT response_id
            FROM blinded_evaluation_responses
            WHERE
                evaluation_set_version = ?
                AND evaluation_id = ?
                AND expert_name = ?
                AND stage = 2
            """,
            (
                str(
                    evaluation_set_version
                ),
                str(
                    evaluation_id
                ),
                expert_name,
            ),
        ).fetchone()

        if stage2 is not None:
            raise ValueError(
                "Stage 2 has already been locked "
                "for this expert and event."
            )

        cursor = connection.execute(
            """
            INSERT INTO
            blinded_evaluation_responses (
                evaluation_set_version,
                evaluation_id,
                expert_name,
                stage,
                dataset_id,
                sender_cell_type,
                receiver_cell_type,
                ligand,
                receptor,
                decision,
                biological_relevance,
                confidence,
                comment,
                method_count,
                total_selected_methods,
                support_fraction,
                method_support,
                consensus_percentile,
                rank_disagreement,
                created_at
            )
            VALUES (
                ?, ?, ?, 2, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                ?, ?, ?, ?, ?, ?, ?
            )
            """,
            (
                str(
                    evaluation_set_version
                ),
                str(
                    evaluation_id
                ),
                expert_name,
                str(
                    dataset_id
                ),
                str(
                    sender_cell_type
                ),
                str(
                    receiver_cell_type
                ),
                str(
                    ligand
                ).upper(),
                str(
                    receptor
                ).upper(),
                str(
                    decision
                ),
                biological_relevance,
                confidence,
                str(
                    comment
                    or ""
                ).strip(),
                (
                    int(
                        method_count
                    )
                    if method_count is not None
                    else None
                ),
                (
                    int(
                        total_selected_methods
                    )
                    if total_selected_methods
                    is not None
                    else None
                ),
                support_fraction,
                str(
                    method_support
                    or ""
                ),
                consensus_percentile,
                rank_disagreement,
                _utc_now(),
            ),
        )

        response_id = int(
            cursor.lastrowid
        )

        for evidence in method_evidence:
            connection.execute(
                """
                INSERT INTO
                blinded_evaluation_method_evidence (
                    response_id,
                    method,
                    present,
                    raw_score,
                    score_name,
                    global_percentile,
                    cell_pair_percentile
                )
                VALUES (
                    ?, ?, ?, ?, ?, ?, ?
                )
                """,
                (
                    response_id,
                    str(
                        evidence[
                            "method"
                        ]
                    ),
                    1
                    if bool(
                        evidence[
                            "present"
                        ]
                    )
                    else 0,
                    evidence.get(
                        "raw_score"
                    ),
                    str(
                        evidence.get(
                            "score_name"
                        )
                        or ""
                    ),
                    evidence.get(
                        "global_percentile"
                    ),
                    evidence.get(
                        "cell_pair_percentile"
                    ),
                ),
            )

        connection.commit()

        return response_id

    except Exception:
        connection.rollback()
        raise

    finally:
        connection.close()
