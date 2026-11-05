"""
Expert assessments (the Human-in-the-Loop step).

Stored in the platform's existing database (``database/ccc_platform.db``)
through the existing ``database.db`` module, so assessments made here and in
the Streamlit pages are the same records.

An assessment is made for the biological event (sender, receiver, ligand
entity, receptor entity), not for one method's row. A snapshot of every
method's score at that moment is stored with it.
"""

from __future__ import annotations

import csv
import io
from contextlib import closing
from typing import Any

from . import config, datasets, index_store

from database.db import get_connection, init_database, save_feedback

REASONS = {
    "reason_pathway": "Plausible signalling pathway",
    "reason_literature": "Supported by literature",
    "reason_disease": "Relevant to the disease or condition",
    "reason_cell_context": "Plausible for these cell types",
    "reason_experimental": "Worth testing experimentally",
    "reason_response": "Receiver shows the expected response",
}

_initialised = False


def _init() -> None:
    global _initialised
    if not _initialised:
        init_database()
        _initialised = True


def gene_level_key(sender: str, receiver: str, ligand: str, receptor: str) -> str:
    """Key format used by the Streamlit pages (kept for compatibility)."""
    return f"{sender}||{receiver}||{str(ligand).upper()}||{str(receptor).upper()}"


def _event_rows(event_key: str) -> list[dict[str, Any]]:
    with closing(index_store._connect()) as conn:
        rows = conn.execute(
            "SELECT id, method, score, score_name, global_percentile, cell_pair_percentile, "
            "sender, receiver, ligand, receptor, method_count, supporting_methods "
            "FROM results WHERE event_key = ? ORDER BY method",
            (event_key,),
        ).fetchall()
    return [dict(r) for r in rows]


def list_for_result(result_id: int) -> list[dict[str, Any]] | None:
    rec = index_store.get_result(result_id)
    if not rec:
        return None
    _init()
    legacy_key = gene_level_key(rec["sender"], rec["receiver"], rec["ligand"], rec["receptor"])
    with closing(get_connection()) as conn:
        rows = conn.execute(
            "SELECT * FROM expert_feedback WHERE dataset_id = ? AND (canonical_entity_key = ? "
            "OR (canonical_entity_key IS NULL AND interaction_key = ?)) "
            "ORDER BY feedback_id DESC",
            (datasets.active().id, rec["event_key"], legacy_key),
        ).fetchall()
    return [_public(dict(r)) for r in rows]


def _public(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": row["feedback_id"],
        "expert": row.get("expert_name") or "",
        "decision": row.get("decision"),
        "relevance": row.get("biological_relevance"),
        "confidence": row.get("confidence"),
        "reasons": [label for key, label in REASONS.items() if row.get(key)],
        "comment": row.get("comment") or "",
        "created_at": row.get("created_at"),
        "scope": row.get("feedback_scope"),
    }


def save(result_id: int, payload: dict[str, Any]) -> dict[str, Any]:
    rec = index_store.get_result(result_id)
    if not rec:
        raise LookupError("Result not found.")

    decision = str(payload.get("decision") or "")
    if decision not in config.DECISIONS:
        raise ValueError("Choose a decision: " + ", ".join(config.DECISIONS) + ".")
    expert = str(payload.get("expert") or "").strip()
    if not expert:
        raise ValueError("Enter your name or expert identifier.")
    try:
        relevance = int(payload.get("relevance"))
        confidence = int(payload.get("confidence"))
    except (TypeError, ValueError):
        raise ValueError("Relevance and confidence must be whole numbers from 1 to 5.")
    if not (1 <= relevance <= 5 and 1 <= confidence <= 5):
        raise ValueError("Relevance and confidence must be from 1 to 5.")
    reasons = set(payload.get("reasons") or [])

    rows = _event_rows(rec["event_key"])
    present = {r["method"]: r for r in rows}
    evidence = []
    for method in index_store.methods():
        r = present.get(method)
        evidence.append({
            "method": method,
            "present": r is not None,
            "raw_score": r["score"] if r else None,
            "score_name": r["score_name"] if r else None,
            "global_percentile": r["global_percentile"] if r else None,
            "cell_pair_percentile": r["cell_pair_percentile"] if r else None,
        })
    percentiles = [r["global_percentile"] for r in rows if r["global_percentile"] is not None]

    _init()
    feedback_id = save_feedback(
        interaction_key=gene_level_key(rec["sender"], rec["receiver"], rec["ligand"], rec["receptor"]),
        canonical_entity_key=rec["event_key"],
        dataset_id=datasets.active().id,
        method="canonical_event",
        id_cp_interaction=None,
        interacting_pair=f"{rec['ligand']}_{rec['receptor']}",
        sender_cell_type=rec["sender"],
        receiver_cell_type=rec["receiver"],
        ligand=rec["ligand"],
        receptor=rec["receptor"],
        computational_score=None,   # no raw score is comparable across methods
        decision=decision,
        biological_relevance=relevance,
        confidence=confidence,
        reason_pathway=int("reason_pathway" in reasons),
        reason_literature=int("reason_literature" in reasons),
        reason_disease=int("reason_disease" in reasons),
        reason_cell_context=int("reason_cell_context" in reasons),
        reason_experimental=int("reason_experimental" in reasons),
        comment=str(payload.get("comment") or "").strip(),
        expert_name=expert,
        feedback_scope="canonical_event",
        method_support=rec.get("supporting_methods"),
        method_count=rec.get("method_count"),
        consensus_percentile=(sum(percentiles) / len(percentiles)) if percentiles else None,
        percentile_gap=(max(percentiles) - min(percentiles)) if len(percentiles) > 1 else 0.0,
        method_evidence=evidence,
    )
    if "reason_response" in reasons and feedback_id:
        # Stored in its own column (added to the database by init_database).
        from database.db import get_connection
        with closing(get_connection()) as conn:
            conn.execute("UPDATE expert_feedback SET reason_response = 1 WHERE rowid = ?", (feedback_id,))
            conn.commit()
    return {"id": feedback_id}


def set_decision(result_id: int, kept: bool, expert: str | None) -> dict[str, Any]:
    """On/off HITL: mark an interaction as kept (default) or removed.

    Only removals are stored (decision = 'Reject'). A kept interaction is the
    default and holds no record, so re-checking a row simply deletes its record.
    """
    rec = index_store.get_result(result_id)
    if not rec:
        raise LookupError("Result not found.")
    expert = (str(expert or "").strip()) or "expert"
    _init()
    legacy_key = gene_level_key(rec["sender"], rec["receiver"], rec["ligand"], rec["receptor"])
    # one decision per event: clear any previous record (from this on/off UI or the old form)
    with closing(get_connection()) as conn:
        conn.execute(
            "DELETE FROM expert_feedback WHERE dataset_id = ? AND (canonical_entity_key = ? "
            "OR (canonical_entity_key IS NULL AND interaction_key = ?))",
            (datasets.active().id, rec["event_key"], legacy_key))
        conn.commit()
    # Store the explicit choice (Accept = kept, Reject = removed). Both are kept
    # because the default now depends on a rank cutoff, so a kept row can be an
    # explicit override of a below-cutoff interaction the expert wants to add.
    feedback_id = save_feedback(
        interaction_key=legacy_key,
        canonical_entity_key=rec["event_key"],
        dataset_id=datasets.active().id,
        method="canonical_event",
        id_cp_interaction=None,
        interacting_pair=f"{rec['ligand']}_{rec['receptor']}",
        sender_cell_type=rec["sender"],
        receiver_cell_type=rec["receiver"],
        ligand=rec["ligand"],
        receptor=rec["receptor"],
        computational_score=None,
        decision="Accept" if kept else "Reject",
        biological_relevance=3,
        confidence=3,
        reason_pathway=0,
        reason_literature=0,
        reason_disease=0,
        reason_cell_context=0,
        reason_experimental=0,
        comment="",
        expert_name=expert,
        feedback_scope="canonical_event",
        method_support=rec.get("supporting_methods"),
        method_count=rec.get("method_count"),
        consensus_percentile=None,
        percentile_gap=0.0,
        method_evidence=[],
    )
    return {"event_key": rec["event_key"], "kept": bool(kept), "id": feedback_id}


def feedback_with_event_keys() -> list[dict[str, Any]]:
    """
    All assessments with the event they belong to. Older records saved by the
    Streamlit pages have only the gene-level key; they are matched to the
    event with the same populations and ligand/receptor labels.
    """
    if not index_store.index_is_current():
        return []          # no results for this dataset (yet)
    _init()
    with closing(get_connection()) as conn:
        rows = [dict(r) for r in conn.execute(
            "SELECT feedback_id, interaction_key, canonical_entity_key, decision, confidence "
            "FROM expert_feedback WHERE dataset_id = ? ORDER BY feedback_id",
            (datasets.active().id,),
        )]
    out = []
    with closing(index_store._connect()) as idx:
        for r in rows:
            key = r["canonical_entity_key"]
            if not key:
                parts = str(r["interaction_key"] or "").split("||")
                if len(parts) == 4:
                    hit = idx.execute(
                        "SELECT event_key FROM results WHERE sender = ? AND receiver = ? "
                        "AND upper(ligand) = ? AND upper(receptor) = ? LIMIT 1",
                        (parts[0], parts[1], parts[2].upper(), parts[3].upper()),
                    ).fetchone()
                    key = hit[0] if hit else None
            if key:
                out.append({**r, "event_key": key})
    return out


def clear_decisions(decision: str | None = None) -> int:
    """Delete stored expert decisions for the active dataset, so the rank cutoff
    governs pre-selection again. decision=None clears all; "Accept" or "Reject"
    clears only that kind. Returns the number of records removed."""
    _init()
    ds = datasets.active().id
    with closing(get_connection()) as conn:
        if decision in ("Accept", "Reject"):
            cur = conn.execute("DELETE FROM expert_feedback WHERE dataset_id = ? AND decision = ?", (ds, decision))
        else:
            cur = conn.execute("DELETE FROM expert_feedback WHERE dataset_id = ?", (ds,))
        n = cur.rowcount
        conn.commit()
    return int(n if n is not None else 0)


def assessed_keys() -> dict[str, str]:
    """event_key -> latest decision, used to mark rows in the table."""
    return {r["event_key"]: r["decision"] for r in feedback_with_event_keys()}


def export_csv() -> str:
    _init()
    with closing(get_connection()) as conn:
        rows = conn.execute("SELECT * FROM expert_feedback WHERE dataset_id = ? ORDER BY feedback_id",
                            (datasets.active().id,)).fetchall()
        cols = [d[0] for d in conn.execute("SELECT * FROM expert_feedback LIMIT 0").description]
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(cols)
    for r in rows:
        writer.writerow([r[c] for c in cols])
    return buf.getvalue()
