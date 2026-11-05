"""
Fast, read-only index of CCI results.

The existing results are large CSV files (80-170 MB). Reading them on every
page view is what made the Streamlit pages slow. Here they are converted
once into an indexed SQLite table (one row per method result), and every
request reads only the rows it needs: one page of the table, or one record.

The index is rebuilt automatically when the source file changes.

Row model (one row = one method's prediction for one candidate event):
    event_key, sender, receiver,
    ligand/receptor entity keys + display labels,
    method, raw score, score name, within-method percentiles,
    method support for the same event (how many methods found it),
    method-provided signalling category (CellChat only).
"""

from __future__ import annotations

import json
import sqlite3
import threading
import time
from contextlib import closing
from pathlib import Path
from typing import Any

import pandas as pd

from . import config, datasets

SCHEMA_VERSION = 5

_BUILD_LOCK = threading.Lock()
_STATUS: dict[str, Any] = {"state": "idle", "message": "", "rows": 0}


# ---------------------------------------------------------------------------
# Entity display helpers
# ---------------------------------------------------------------------------
def entity_display(entity_key: str | None, fallback: str | None = None) -> dict[str, Any]:
    """
    Turn a canonical entity key into a display record without losing identity.

    GENE:APOE               -> gene, APOE
    COMPLEX:TREM2+TYROBP    -> complex, subunits TREM2 and TYROBP
    CPDB_PROXY:<name>       -> small molecule estimated from enzyme expression
    CPDB_OTHER:<name>       -> CellPhoneDB entity not resolved to verified subunits
    """
    key = (entity_key or "").strip()
    prefix, _, body = key.partition(":")
    prefix = prefix.upper()

    if prefix == "GENE" and body:
        return {"key": key, "kind": "gene", "label": body, "subunits": [body]}

    if prefix == "COMPLEX" and body:
        subunits = [s for s in body.split("+") if s]
        return {
            "key": key,
            "kind": "complex" if len(subunits) > 1 else "gene",
            "label": "+".join(subunits),
            "subunits": subunits,
        }

    if prefix == "CPDB_PROXY" and body:
        return {"key": key, "kind": "proxy", "label": body, "subunits": []}

    if prefix == "CPDB_OTHER" and body:
        return {"key": key, "kind": "other", "label": body, "subunits": []}

    label = (fallback or body or key or "").replace("COMPLEX:", "")
    return {"key": key, "kind": "unknown", "label": label, "subunits": []}


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
def _source_signature(ds: "datasets.Dataset | None" = None) -> dict[str, Any]:
    ds = ds or datasets.active()
    sig: dict[str, Any] = {"schema": SCHEMA_VERSION, "dataset": ds.id}
    for name, path in (("consensus", ds.consensus_file), ("cellchat", ds.cellchat_file)):
        p = Path(path) if path else None
        sig[name] = [str(p), p.stat().st_mtime if p and p.exists() else None]
    return sig


def _read_meta(conn: sqlite3.Connection) -> dict[str, Any] | None:
    try:
        row = conn.execute("SELECT value FROM meta WHERE key='signature'").fetchone()
    except sqlite3.Error:
        return None
    return json.loads(row[0]) if row else None


def methods() -> list[str]:
    """Tools whose results this dataset has, in the platform's fixed order.

    Read from the index each time (a few milliseconds), never cached: after a
    run with different tools the index is rebuilt and the list must follow.
    """
    found = None
    try:
        with closing(_connect()) as conn:
            row = conn.execute("SELECT value FROM meta WHERE key='methods'").fetchone()
            if row:
                found = json.loads(row[0])
    except Exception:
        found = None
    return list(found) if found else []


available_methods = methods      # alias for code where a parameter is called "methods"


def index_is_current() -> bool:
    ds = datasets.active_or_none()
    if ds is None:
        return False
    if not ds.index_db.exists():
        return False
    # closing(): on Windows an open connection keeps the file locked.
    with closing(sqlite3.connect(ds.index_db)) as conn:
        return _read_meta(conn) == _source_signature(ds)


def _cellchat_categories(cellchat_file: Path | None) -> pd.DataFrame:
    if not cellchat_file or not cellchat_file.exists():
        return pd.DataFrame(
            columns=["sender", "receiver", "m_ligand", "m_receptor",
                     "signaling_pathway", "signaling_annotation", "signaling_evidence"]
        )
    cc = pd.read_csv(
        cellchat_file,
        usecols=["sender_cell_type", "receiver_cell_type", "ligand", "receptor",
                 "pathway_name", "annotation", "evidence"],
        dtype=str,
    )
    cc = cc.rename(columns={
        "sender_cell_type": "sender",
        "receiver_cell_type": "receiver",
        "ligand": "m_ligand",
        "receptor": "m_receptor",
        "pathway_name": "signaling_pathway",
        "annotation": "signaling_annotation",
        "evidence": "signaling_evidence",
    })
    return cc.drop_duplicates(["sender", "receiver", "m_ligand", "m_receptor"])


def build_index(force: bool = False) -> dict[str, Any]:
    """Convert the source CSV into the indexed table. Safe to call repeatedly."""
    with _BUILD_LOCK:
        ds = datasets.active_or_none()
        if ds is None:
            _STATUS.update(state="no_dataset", message="No dataset is loaded.", rows=0)
            return dict(_STATUS)
        if not force and index_is_current():
            _STATUS.update(state="ready", message="Index is up to date.", dataset=ds.id)
            return dict(_STATUS)

        if not ds.consensus_file or not ds.consensus_file.exists():
            _STATUS.update(
                state="missing", dataset=ds.id,
                message=(f"No CCI results for “{ds.label}” yet. "
                         + (f"Expected: {ds.consensus_file}" if ds.consensus_file else
                            "Run the CCI methods for this dataset first.")),
            )
            return dict(_STATUS)

        started = time.time()
        _STATUS.update(state="building", message="Reading results file…", rows=0, dataset=ds.id)

        src = pd.read_csv(ds.consensus_file, low_memory=False)

        frames = []
        present_methods = [name for prefix, name in config.METHODS.items() if f"{prefix}_present" in src.columns]
        for prefix, method_name in config.METHODS.items():
            present_col = f"{prefix}_present"
            if present_col not in src.columns:
                continue
            present = src[present_col].astype(str).str.lower().isin(["true", "1"])
            part = src.loc[present]
            frames.append(pd.DataFrame({
                "event_key": part["interaction_key"],
                "sender": part["sender_cell_type"],
                "receiver": part["receiver_cell_type"],
                "ligand_entity_key": part["ligand_entity_key"],
                "receptor_entity_key": part["receptor_entity_key"],
                "method": method_name,
                "m_ligand": part.get(f"{prefix}_ligand_label"),
                "m_receptor": part.get(f"{prefix}_receptor_label"),
                "score": pd.to_numeric(part.get(f"{prefix}_score"), errors="coerce"),
                "score_name": part.get(f"{prefix}_score_name"),
                "global_percentile": pd.to_numeric(part.get(f"{prefix}_global_percentile"), errors="coerce"),
                "cell_pair_percentile": pd.to_numeric(part.get(f"{prefix}_cell_pair_percentile"), errors="coerce"),
                "method_count": pd.to_numeric(part["method_count"], errors="coerce").fillna(1).astype(int),
                "supporting_methods": part.get("supporting_methods"),
            }))

        _STATUS.update(message="Adding method-specific fields…")
        long = pd.concat(frames, ignore_index=True)

        cc = _cellchat_categories(ds.cellchat_file)
        cc["method"] = "CellChat"
        long = long.merge(cc, how="left", on=["method", "sender", "receiver", "m_ligand", "m_receptor"])

        lig = long["ligand_entity_key"].map(lambda k: entity_display(k)["label"])
        rec = long["receptor_entity_key"].map(lambda k: entity_display(k)["label"])
        long["ligand"] = lig
        long["receptor"] = rec
        long["ligand_kind"] = long["ligand_entity_key"].map(lambda k: entity_display(k)["kind"])
        long["receptor_kind"] = long["receptor_entity_key"].map(lambda k: entity_display(k)["kind"])
        long["search_text"] = (
            lig.str.lower() + " " + rec.str.lower() + " "
            + long["m_ligand"].fillna("").str.lower() + " "
            + long["m_receptor"].fillna("").str.lower()
        )

        _STATUS.update(message="Writing index…")
        tmp_path = ds.index_db.with_suffix(".building")
        if tmp_path.exists():
            tmp_path.unlink()

        conn = sqlite3.connect(tmp_path)
        try:
            long.to_sql("results", conn, index=True, index_label="id", chunksize=20000)
            for col in ("method", "sender", "receiver", "ligand", "receptor",
                        "global_percentile", "method_count", "event_key"):
                conn.execute(f"CREATE INDEX idx_{col} ON results({col})")

            # One row per interaction (event), with each tool's own rank side by side.
            _STATUS.update(message="Grouping predictions by interaction…")
            lr = long.reset_index().rename(columns={"index": "id"})
            wide = lr.pivot_table(index="event_key", columns="method", values="global_percentile", aggfunc="max")
            ids = lr.pivot_table(index="event_key", columns="method", values="id", aggfunc="min")
            first = lr.drop_duplicates("event_key").set_index("event_key")
            ev = first[["sender", "receiver", "ligand", "receptor", "ligand_kind", "receptor_kind",
                        "method_count", "supporting_methods", "search_text"]].copy()
            cc = lr[lr["method"] == "CellChat"].drop_duplicates("event_key").set_index("event_key")
            ev["signaling_pathway"] = cc["signaling_pathway"].reindex(ev.index)
            order = present_methods
            for m in order:
                slug = m.lower()
                ev[f"rank_{slug}"] = wide[m].reindex(ev.index) if m in wide else None
                ev[f"id_{slug}"] = ids[m].reindex(ev.index) if m in ids else None
            # Representative prediction for opening the documentation: first tool, in a fixed order.
            ev["open_id"] = ev[[f"id_{m.lower()}" for m in order]].bfill(axis=1).iloc[:, 0].astype("Int64")
            ev.reset_index().to_sql("events", conn, index=True, index_label="row", chunksize=20000)
            for col in ("sender", "receiver", "ligand", "receptor", "method_count", "event_key"):
                conn.execute(f"CREATE INDEX idx_ev_{col} ON events({col})")
            conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
            conn.execute(
                "INSERT INTO meta VALUES ('methods', ?)",
                (json.dumps(present_methods),),
            )
            conn.execute(
                "INSERT INTO meta VALUES ('signature', ?)",
                (json.dumps(_source_signature(ds)),),
            )
            conn.execute(
                "INSERT INTO meta VALUES ('built_at', ?)",
                (time.strftime("%Y-%m-%d %H:%M:%S"),),
            )
            conn.execute("ANALYZE")
            conn.commit()
        finally:
            conn.close()  # must be closed before the file can be renamed on Windows

        tmp_path.replace(ds.index_db)
        _STATUS.update(
            state="ready",
            rows=int(len(long)),
            message=f"Indexed {len(long):,} results in {time.time() - started:.0f} s.",
        )
        return dict(_STATUS)


def status() -> dict[str, Any]:
    ds = datasets.active_or_none()
    if ds is None:
        return {"state": "no_dataset", "message": "No dataset is loaded.", "rows": 0}
    # A "ready" (or "missing") state can go stale when new results are written for the
    # same dataset, for example after running the tools again: check it, and rebuild.
    if (_STATUS.get("dataset") == ds.id and _STATUS["state"] in ("ready", "missing")
            and not _BUILD_LOCK.locked() and not index_is_current()):
        _STATUS.update(state="idle")
    if _STATUS.get("dataset") != ds.id or _STATUS["state"] == "idle":
        if index_is_current():
            _STATUS.update(state="ready", message="Index is up to date.", dataset=ds.id)
        elif not ds.consensus_file or not ds.consensus_file.exists():
            build_index()          # sets the "missing" state and message
        elif _STATUS["state"] != "building":
            threading.Thread(target=build_index, daemon=True).start()
            _STATUS.update(state="building", message="", dataset=ds.id)
    return dict(_STATUS)


# ---------------------------------------------------------------------------
# Queries
# ---------------------------------------------------------------------------
def _connect() -> sqlite3.Connection:
    # as_uri() gives file:///C:/... on Windows; a plain path with backslashes is not a valid URI.
    conn = sqlite3.connect(datasets.active().index_db.as_uri() + "?mode=ro", uri=True, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


SORTS = {
    "percentile": "global_percentile DESC",
    "support": "method_count DESC, global_percentile DESC",
    "score": "score DESC",            # only allowed when one method is selected
    "sender": "sender ASC, global_percentile DESC",
    "receiver": "receiver ASC, global_percentile DESC",
    "ligand": "ligand ASC, global_percentile DESC",
    "receptor": "receptor ASC, global_percentile DESC",
}

LIST_COLUMNS = (
    "id, event_key, sender, receiver, ligand, receptor, ligand_kind, receptor_kind, "
    "method, score, score_name, global_percentile, method_count, supporting_methods, "
    "signaling_pathway, signaling_annotation"
)


def query_results(
    *,
    methods: list[str] | None = None,
    sender: str | None = None,
    receiver: str | None = None,
    q: str | None = None,
    min_percentile: float | None = None,
    min_methods: int | None = None,
    sort: str = "percentile",
    page: int = 1,
    page_size: int = 50,
) -> dict[str, Any]:
    where: list[str] = []
    params: list[Any] = []

    if methods:
        where.append(f"method IN ({','.join('?' * len(methods))})")
        params.extend(methods)
    if sender:
        where.append("sender = ?")
        params.append(sender)
    if receiver:
        where.append("receiver = ?")
        params.append(receiver)
    if q:
        for token in q.lower().split():
            where.append("search_text LIKE ?")
            params.append(f"%{token}%")
    if min_percentile is not None:
        where.append("global_percentile >= ?")
        params.append(float(min_percentile))
    if min_methods:
        where.append("method_count >= ?")
        params.append(int(min_methods))

    # Raw scores are method-specific: sorting by them is only meaningful
    # within a single method.
    single_method = bool(methods) and len(methods) == 1
    if sort == "score" and not single_method:
        sort = "percentile"
    order = SORTS.get(sort, SORTS["percentile"])

    page_size = max(10, min(int(page_size), 200))
    page = max(1, int(page))
    clause = f"WHERE {' AND '.join(where)}" if where else ""

    with closing(_connect()) as conn:
        total = conn.execute(f"SELECT COUNT(*) FROM results {clause}", params).fetchone()[0]
        rows = conn.execute(
            f"SELECT {LIST_COLUMNS} FROM results {clause} ORDER BY {order}, id LIMIT ? OFFSET ?",
            [*params, page_size, (page - 1) * page_size],
        ).fetchall()

    return {
        "total": int(total),
        "page": page,
        "page_size": page_size,
        "pages": max(1, -(-int(total) // page_size)),
        "sort": sort,
        "score_sort_allowed": single_method,
        "rows": [dict(r) for r in rows],
    }


EVENT_SORTS = {
    "tools": "method_count DESC, ligand ASC, receptor ASC",
    "sender": "sender ASC, receiver ASC, ligand ASC",
    "receiver": "receiver ASC, sender ASC, ligand ASC",
    "ligand": "ligand ASC, receptor ASC, sender ASC",
    "receptor": "receptor ASC, ligand ASC, sender ASC",
}


def query_events(
    *,
    methods: list[str] | None = None,
    sender: str | None = None,
    receiver: str | None = None,
    q: str | None = None,
    min_percentile: float | None = None,
    min_methods: int | None = None,
    sort: str = "tools",
    page: int = 1,
    page_size: int = 50,
) -> dict[str, Any]:
    """One row per interaction. Tool filters and rank thresholds apply to each tool on its own."""
    tools = available_methods()
    chosen = [m for m in (methods or []) if m in tools] or tools
    where: list[str] = []
    params: list[Any] = []
    # reported by at least one of the chosen tools (and, with a threshold, in that tool's own top part)
    parts = []
    for m in chosen:
        col = f"rank_{m.lower()}"
        if min_percentile:
            parts.append(f"{col} >= ?"); params.append(float(min_percentile))
        else:
            parts.append(f"{col} IS NOT NULL")
    where.append("(" + " OR ".join(parts) + ")")
    if sender:
        where.append("sender = ?"); params.append(sender)
    if receiver:
        where.append("receiver = ?"); params.append(receiver)
    if q:
        for token in q.lower().split():
            where.append("search_text LIKE ?"); params.append(f"%{token}%")
    if min_methods:
        where.append("method_count >= ?"); params.append(int(min_methods))

    if sort.startswith("rank_") and sort[5:] in [t.lower() for t in tools]:
        order = f"{sort} IS NULL, {sort} DESC, ligand ASC"
    else:
        sort = sort if sort in EVENT_SORTS else "tools"
        order = EVENT_SORTS[sort]
    # Normal paginated views stay capped at 200 rows; the aggregation endpoints
    # (pathway panel, bulk curation) ask for a very large page to scan every
    # interaction at once, so honour large explicit page sizes too.
    page_size = max(10, min(int(page_size), 500000))
    page = max(1, int(page))
    clause = "WHERE " + " AND ".join(where)
    cols = ("row, event_key, sender, receiver, ligand, receptor, ligand_kind, receptor_kind, method_count, "
            "supporting_methods, signaling_pathway, open_id, " +
            ", ".join(f"rank_{t.lower()}, id_{t.lower()}" for t in tools))
    with closing(_connect()) as conn:
        total = conn.execute(f"SELECT COUNT(*) FROM events {clause}", params).fetchone()[0]
        rows = conn.execute(f"SELECT {cols} FROM events {clause} ORDER BY {order}, row LIMIT ? OFFSET ?",
                            [*params, page_size, (page - 1) * page_size]).fetchall()
    out = []
    for r in rows:
        d = dict(r)
        d["tools"] = [{"method": t, "rank": d.pop(f"rank_{t.lower()}"),
                       "id": (int(d[f"id_{t.lower()}"]) if d.get(f"id_{t.lower()}") is not None else None)} for t in tools]
        for t in tools:
            d.pop(f"id_{t.lower()}", None)
        out.append(d)
    return {"total": int(total), "page": page, "page_size": page_size,
            "pages": max(1, -(-int(total) // page_size)), "sort": sort, "methods": tools, "rows": out}


def get_result(result_id: int) -> dict[str, Any] | None:
    with closing(_connect()) as conn:
        row = conn.execute("SELECT * FROM results WHERE id = ?", (int(result_id),)).fetchone()
        if not row:
            return None
        record = dict(row)
        others = conn.execute(
            "SELECT id, method, score, score_name, global_percentile FROM results "
            "WHERE event_key = ? AND id != ? ORDER BY method",
            (record["event_key"], record["id"]),
        ).fetchall()
    record.pop("search_text", None)
    record["other_methods"] = [dict(o) for o in others]
    return record


def facets() -> dict[str, Any]:
    with closing(_connect()) as conn:
        methods = [r[0] for r in conn.execute("SELECT DISTINCT method FROM results ORDER BY method")]
        senders = [r[0] for r in conn.execute("SELECT DISTINCT sender FROM results ORDER BY sender")]
        receivers = [r[0] for r in conn.execute("SELECT DISTINCT receiver FROM results ORDER BY receiver")]
        counts = {r[0]: r[1] for r in conn.execute("SELECT method, COUNT(*) FROM results GROUP BY method")}
        built = conn.execute("SELECT value FROM meta WHERE key='built_at'").fetchone()
    return {
        "methods": methods,
        "method_counts": counts,
        "senders": senders,
        "receivers": receivers,
        "built_at": built[0] if built else None,
    }
