"""
Expert review queue.

Groups candidate interactions by WHICH tools report them, each tool with its
own criteria:

  Expert revisit            earlier assessments conflict, are uncertain or were
                            made with low confidence
  Found by all tools        reported by CellPhoneDB, NATMI and CellChat
  Found by two tools        reported by two of the three
  Found by one tool only    reported by one tool, at the very top of its own list

Ranks of different tools are never combined or compared; each tool's rank
(within its own results) is shown separately. Interactions whose assessments
already agree ("resolved") leave the queue.
"""

from __future__ import annotations

import threading
from contextlib import closing
from typing import Any

import pandas as pd

from . import assessments, config, index_store


DEFAULTS = {
    "single": 99.9,             # a single tool's own rank, for "found by one tool only"
}

# Each tool is taken with its own criteria: events are grouped by WHICH tools
# report them. Ranks of different tools are never combined or compared.
CATEGORIES = [
    ("revisit", "Expert revisit",
     "Earlier assessments of this event conflict, are uncertain or were made with low confidence."),
    ("all", "Found by all tools",
     "Every tool that was run reports this interaction, each with its own criteria."),
    ("two", "Found by several tools",
     "More than one tool, but not all, report this interaction."),
    ("single", "Found by one tool only",
     "Only one tool reports this interaction, at the very top of its own list. Expert review can judge whether it is plausible or specific to that tool or its database."),
]

_cache: dict[str, Any] = {"sig": None, "events": None}
_lock = threading.Lock()


def _events() -> pd.DataFrame:
    """One row per candidate event with each method's rank (cached)."""
    sig = index_store._source_signature()
    with _lock:
        if _cache["sig"] == sig and _cache["events"] is not None:
            return _cache["events"]
        with closing(index_store._connect()) as conn:
            long = pd.read_sql_query(
                "SELECT id, event_key, sender, receiver, ligand, receptor, method, "
                "global_percentile FROM results",
                conn,
            )
        wide = long.pivot_table(index="event_key", columns="method",
                                values="global_percentile", aggfunc="max")
        methods = index_store.methods()
        for m in methods:
            if m not in wide.columns:
                wide[m] = float("nan")
        wide = wide[methods]

        best = long.sort_values("global_percentile", ascending=False).drop_duplicates("event_key")
        best = best.set_index("event_key")[["id", "sender", "receiver", "ligand", "receptor"]]

        ev = best.join(wide)
        ranks = ev[methods]
        ev["method_count"] = ranks.notna().sum(axis=1)
        ev["best_rank"] = ranks.max(axis=1)
        ev = ev.reset_index()
        _cache.update(sig=sig, events=ev)
        return ev


def _expert_status() -> pd.DataFrame:
    fb = pd.DataFrame(assessments.feedback_with_event_keys(),
                      columns=["event_key", "decision", "confidence"])
    if fb.empty:
        return pd.DataFrame(columns=["event_key", "n_reviews", "status"])
    grp = fb.groupby("event_key")
    out = pd.DataFrame({
        "n_reviews": grp.size(),
        "positive": grp["decision"].apply(lambda d: d.isin(["Prioritize", "Accept"]).sum()),
        "reject": grp["decision"].apply(lambda d: (d == "Reject").sum()),
        "uncertain": grp["decision"].apply(lambda d: (d == "Uncertain").sum()),
        "mean_confidence": grp["confidence"].mean(),
        "decisions": grp["decision"].apply(lambda d: ", ".join(d)),
    })

    def status(r):
        if r.positive > 0 and r.reject > 0:
            return "Conflicting assessments"
        if r.uncertain > 0:
            return "Expert uncertainty"
        if r.mean_confidence <= 2:
            return "Low confidence"
        return "Resolved"

    out["status"] = out.apply(status, axis=1)
    return out.reset_index()


def queue(
    *,
    thresholds: dict[str, float] | None = None,
    category: str | None = None,
    sort_by: str | None = None,
    page: int = 1,
    page_size: int = 50,
) -> dict[str, Any]:
    t = {**DEFAULTS, **{k: v for k, v in (thresholds or {}).items() if v is not None}}
    ev = _events()
    n_methods = len(index_store.methods())

    st = _expert_status()
    if st.empty or "event_key" not in st.columns:      # no assessments yet for this dataset
        st = pd.DataFrame(columns=["event_key", "n_reviews", "status", "decisions", "revisit", "revisit_reason"])
    ev = ev.merge(st[["event_key", "n_reviews", "status", "decisions"]], on="event_key", how="left")
    ev["n_reviews"] = ev["n_reviews"].fillna(0).astype(int)
    ev["status"] = ev["status"].fillna("Not reviewed")

    revisit = ev["status"].isin(["Conflicting assessments", "Expert uncertainty", "Low confidence"])
    open_ = ev["status"] == "Not reviewed"
    if n_methods == 1:      # one tool: the queue is the top of its own list, plus expert revisits
        labels = {"single": ("Top of the tool’s results", "The interactions at the very top of the tool’s own list (rank at least the threshold).")}
    else:
        labels = {}
    masks = {
        "revisit": revisit,
        # With a single tool, "found by all tools" would be every interaction: leave it empty.
        "all": open_ & (ev["method_count"] == n_methods) & (n_methods > 1),
        "two": open_ & (ev["method_count"] >= 2) & (ev["method_count"] < n_methods),
        "single": open_ & (ev["method_count"] == 1) & (ev["best_rank"] >= t["single"]),
    }
    # First matching category wins, in priority order.
    ev["category"] = None
    assigned = pd.Series(False, index=ev.index)
    for key, _label, _why in CATEGORIES:
        m = masks[key] & ~assigned
        ev.loc[m, "category"] = key
        assigned |= m

    counts = {key: int((ev["category"] == key).sum()) for key, _, _ in CATEGORIES}
    q = ev[ev["category"].notna()]
    if category in counts:
        q = q[q["category"] == category]

    order = {key: i for i, (key, _, _) in enumerate(CATEGORIES)}
    q = q.assign(_o=q["category"].map(order))
    # Order: category, then (optionally) one tool's own rank, else alphabetically.
    if sort_by in index_store.methods():
        q = q.assign(_r=q[sort_by].fillna(-1)).sort_values(["_o", "_r", "ligand", "receptor"], ascending=[True, False, True, True])
    else:
        q = q.sort_values(["_o", "ligand", "receptor", "sender", "receiver"])

    total = len(q)
    page_size = max(10, min(int(page_size), 200))
    pages = max(1, -(-total // page_size))
    page = max(1, min(int(page), pages))
    part = q.iloc[(page - 1) * page_size: page * page_size]

    methods = index_store.methods()
    rows = []
    for r in part.itertuples(index=False):
        rd = r._asdict()
        rows.append({
            "event_key": rd["event_key"],
            "result_id": int(rd["id"]),
            "sender": rd["sender"], "receiver": rd["receiver"],
            "ligand": rd["ligand"], "receptor": rd["receptor"],
            "ranks": {m: (None if pd.isna(rd[m]) else round(float(rd[m]), 1)) for m in methods},
            "method_count": int(rd["method_count"]),
            "category": rd["category"],
            "status": rd["status"],
            "n_reviews": int(rd["n_reviews"]),
            "decisions": rd["decisions"] if isinstance(rd["decisions"], str) else "",
        })

    reviewed = ev[ev["n_reviews"] > 0]
    return {
        "thresholds": t,
        "categories": [{"key": k, "label": labels.get(k, (l, w))[0], "why": labels.get(k, (l, w))[1], "count": counts[k]}
                       for k, l, w in CATEGORIES if not (n_methods == 1 and k in ("all", "two"))],
        "methods": methods,
        "category": category if category in counts else None,
        "total": total, "page": page, "pages": pages,
        "rows": rows,
        "summary": {
            "events": int(len(ev)),
            "reviewed": int(len(reviewed)),
            "resolved": int((reviewed["status"] == "Resolved").sum()),
            "revisit": counts["revisit"],
        },
        "sort_by": sort_by if sort_by in index_store.methods() else None,
        "definitions": {
            "rank": "Rank of the interaction within that tool's own results (100 = highest). Ranks of different tools are not compared.",
        },
    }
