from pathlib import Path
import sqlite3
from datetime import datetime


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DB_FILE = (
    PROJECT_ROOT
    / "database"
    / "ccc_platform.db"
)


# --------------------------------------------------
# Connection
# --------------------------------------------------
def get_connection():
    """
    Create and return a connection to the SQLite database.
    """

    connection = sqlite3.connect(
        DB_FILE
    )

    connection.row_factory = sqlite3.Row

    # Required for ON DELETE CASCADE in child tables.
    connection.execute(
        "PRAGMA foreign_keys = ON"
    )

    return connection


# --------------------------------------------------
# Migration helper
# --------------------------------------------------
def _get_table_columns(
    connection,
    table_name,
):
    """
    Return the column names currently present
    in a SQLite table.
    """

    rows = connection.execute(
        f"PRAGMA table_info({table_name})"
    ).fetchall()

    return {
        row["name"]
        for row in rows
    }


def _add_column_if_missing(
    connection,
    table_name,
    column_name,
    column_definition,
):
    """
    Add a column only if it does not already exist.

    This allows us to migrate an existing database
    without deleting previous feedback records.
    """

    existing_columns = _get_table_columns(
        connection,
        table_name,
    )

    if column_name not in existing_columns:

        connection.execute(
            f"""
            ALTER TABLE {table_name}
            ADD COLUMN {column_name}
            {column_definition}
            """
        )


# --------------------------------------------------
# Database initialization / migration
# --------------------------------------------------
def init_database():
    """
    Create the feedback database structure and
    safely migrate older versions.

    Existing feedback records are preserved.
    """

    connection = get_connection()

    cursor = connection.cursor()


    # --------------------------------------------------
    # Legacy / parent feedback table
    # --------------------------------------------------
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS expert_feedback (

            feedback_id INTEGER PRIMARY KEY AUTOINCREMENT,

            interaction_key TEXT NOT NULL,

            dataset_id TEXT NOT NULL,

            method TEXT NOT NULL,

            id_cp_interaction TEXT,

            interacting_pair TEXT,

            sender_cell_type TEXT,

            receiver_cell_type TEXT,

            ligand TEXT,

            receptor TEXT,

            computational_score REAL,

            decision TEXT NOT NULL,

            biological_relevance INTEGER NOT NULL,

            confidence INTEGER NOT NULL,

            reason_pathway INTEGER DEFAULT 0,

            reason_literature INTEGER DEFAULT 0,

            reason_disease INTEGER DEFAULT 0,

            reason_cell_context INTEGER DEFAULT 0,

            reason_experimental INTEGER DEFAULT 0,

            comment TEXT,

            expert_name TEXT,

            created_at TEXT NOT NULL
        )
        """
    )


    # --------------------------------------------------
    # New event-level fields
    #
    # These are added safely to an existing DB.
    # --------------------------------------------------
    _add_column_if_missing(
        connection,
        "expert_feedback",
        "feedback_scope",
        "TEXT DEFAULT 'method_prediction'",
    )

    _add_column_if_missing(
        connection,
        "expert_feedback",
        "method_support",
        "TEXT",
    )

    _add_column_if_missing(
        connection,
        "expert_feedback",
        "method_count",
        "INTEGER",
    )

    # Expert judgement that the receiver shows the expected downstream response.
    _add_column_if_missing(
        connection,
        "expert_feedback",
        "reason_response",
        "INTEGER DEFAULT 0",
    )

    _add_column_if_missing(
        connection,
        "expert_feedback",
        "consensus_percentile",
        "REAL",
    )

    _add_column_if_missing(
        connection,
        "expert_feedback",
        "percentile_gap",
        "REAL",
    )

    _add_column_if_missing(
        connection,
        "expert_feedback",
        "canonical_entity_key",
        "TEXT",
    )


    # --------------------------------------------------
    # Child table:
    # computational evidence available to the expert
    # at the moment the feedback was submitted.
    #
    # One feedback can therefore have evidence from
    # any number of CCC methods.
    # --------------------------------------------------
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS feedback_method_evidence (

            evidence_id INTEGER PRIMARY KEY AUTOINCREMENT,

            feedback_id INTEGER NOT NULL,

            method TEXT NOT NULL,

            present INTEGER NOT NULL,

            raw_score REAL,

            score_name TEXT,

            global_percentile REAL,

            cell_pair_percentile REAL,

            created_at TEXT NOT NULL,

            FOREIGN KEY (feedback_id)
                REFERENCES expert_feedback(feedback_id)
                ON DELETE CASCADE
        )
        """
    )


    # --------------------------------------------------
    # Helpful indexes
    # --------------------------------------------------
    cursor.execute(
        """
        CREATE INDEX IF NOT EXISTS
        idx_feedback_interaction_key
        ON expert_feedback(interaction_key)
        """
    )

    cursor.execute(
        """
        CREATE INDEX IF NOT EXISTS
        idx_feedback_canonical_entity_key
        ON expert_feedback(canonical_entity_key)
        """
    )

    cursor.execute(
        """
        CREATE INDEX IF NOT EXISTS
        idx_feedback_dataset
        ON expert_feedback(dataset_id)
        """
    )

    cursor.execute(
        """
        CREATE INDEX IF NOT EXISTS
        idx_evidence_feedback_id
        ON feedback_method_evidence(feedback_id)
        """
    )


    # Explicitly mark any old records.
    cursor.execute(
        """
        UPDATE expert_feedback
        SET feedback_scope = 'method_prediction'
        WHERE feedback_scope IS NULL
        """
    )


    connection.commit()
    connection.close()


# --------------------------------------------------
# Save feedback
# --------------------------------------------------
def save_feedback(
    interaction_key,
    dataset_id,
    method,
    id_cp_interaction,
    interacting_pair,
    sender_cell_type,
    receiver_cell_type,
    ligand,
    receptor,
    computational_score,
    decision,
    biological_relevance,
    confidence,
    reason_pathway,
    reason_literature,
    reason_disease,
    reason_cell_context,
    reason_experimental,
    comment,
    expert_name,

    # New optional fields.
    feedback_scope="canonical_event",
    method_support=None,
    method_count=None,
    consensus_percentile=None,
    percentile_gap=None,
    method_evidence=None,
    canonical_entity_key=None,
):
    """
    Store one expert evaluation.

    The parent record represents the expert judgement.

    method_evidence optionally stores a snapshot of
    each CCC method's computational evidence available
    when the expert submitted the judgement.
    """

    init_database()

    connection = get_connection()

    cursor = connection.cursor()

    created_at = datetime.now().isoformat(
        timespec="seconds"
    )


    cursor.execute(
        """
        INSERT INTO expert_feedback (

            interaction_key,
            canonical_entity_key,
            dataset_id,
            method,
            id_cp_interaction,
            interacting_pair,
            sender_cell_type,
            receiver_cell_type,
            ligand,
            receptor,
            computational_score,

            decision,
            biological_relevance,
            confidence,

            reason_pathway,
            reason_literature,
            reason_disease,
            reason_cell_context,
            reason_experimental,

            comment,
            expert_name,

            created_at,

            feedback_scope,
            method_support,
            method_count,
            consensus_percentile,
            percentile_gap

        )

        VALUES (
            ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?, ?, ?,
            ?, ?,
            ?,
            ?, ?, ?, ?, ?
        )
        """,

        (
            interaction_key,
            canonical_entity_key,
            dataset_id,
            method,
            id_cp_interaction,
            interacting_pair,
            sender_cell_type,
            receiver_cell_type,
            ligand,
            receptor,
            computational_score,

            decision,
            biological_relevance,
            confidence,

            int(reason_pathway),
            int(reason_literature),
            int(reason_disease),
            int(reason_cell_context),
            int(reason_experimental),

            comment,
            expert_name,

            created_at,

            feedback_scope,
            method_support,
            method_count,
            consensus_percentile,
            percentile_gap,
        ),
    )


    feedback_id = cursor.lastrowid


    # --------------------------------------------------
    # Store method-level evidence snapshot
    # --------------------------------------------------
    if method_evidence:

        for evidence in method_evidence:

            cursor.execute(
                """
                INSERT INTO feedback_method_evidence (

                    feedback_id,
                    method,
                    present,
                    raw_score,
                    score_name,
                    global_percentile,
                    cell_pair_percentile,
                    created_at

                )

                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,

                (
                    feedback_id,

                    evidence.get(
                        "method"
                    ),

                    int(
                        bool(
                            evidence.get(
                                "present",
                                False,
                            )
                        )
                    ),

                    evidence.get(
                        "raw_score"
                    ),

                    evidence.get(
                        "score_name"
                    ),

                    evidence.get(
                        "global_percentile"
                    ),

                    evidence.get(
                        "cell_pair_percentile"
                    ),

                    created_at,
                ),
            )


    connection.commit()
    connection.close()

    return feedback_id


# --------------------------------------------------
# Load parent feedback
# --------------------------------------------------
def load_feedback():
    """
    Return all stored expert feedback.
    """

    init_database()

    connection = get_connection()

    rows = connection.execute(
        """
        SELECT *
        FROM expert_feedback
        ORDER BY feedback_id DESC
        """
    ).fetchall()

    connection.close()

    return [
        dict(row)
        for row in rows
    ]


# --------------------------------------------------
# Load method evidence
# --------------------------------------------------
def load_method_evidence(
    feedback_id=None
):
    """
    Return method-level evidence snapshots.

    If feedback_id is supplied, return evidence only
    for that specific expert evaluation.
    """

    init_database()

    connection = get_connection()


    if feedback_id is None:

        rows = connection.execute(
            """
            SELECT *
            FROM feedback_method_evidence
            ORDER BY evidence_id DESC
            """
        ).fetchall()

    else:

        rows = connection.execute(
            """
            SELECT *
            FROM feedback_method_evidence
            WHERE feedback_id = ?
            ORDER BY evidence_id ASC
            """,
            (
                feedback_id,
            ),
        ).fetchall()


    connection.close()

    return [
        dict(row)
        for row in rows
    ]

def wipe_all_feedback():
    """Delete every recorded assessment. Irreversible. Used by demo/ephemeral mode."""
    with get_connection() as connection:
        connection.execute("DELETE FROM feedback_method_evidence")
        connection.execute("DELETE FROM expert_feedback")
        connection.commit()
