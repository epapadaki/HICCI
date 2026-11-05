"""
Communication network between cell populations.

An edge sender -> receiver counts the candidate interactions (distinct events)
that pass the filters. Two sources:
  consensus  interactions reported by at least `min_methods` tools, each tool
             with its own criteria (ranks of different tools are not combined)
  <method>   one tool's predictions, optionally limited to the top of its own list
Counts summarise candidates; they are not measured signalling strength.
"""

from __future__ import annotations

from contextlib import closing
from typing import Any

import pandas as pd

from . import config, index_store, review_queue


def _filtered(source: str, min_methods: int, min_rank: float, cell_types: list[str] | None) -> pd.DataFrame:
    methods = index_store.methods()
    if source in methods:
        with closing(index_store._connect()) as conn:
            df = pd.read_sql_query(
                "SELECT id, event_key, sender, receiver, ligand, receptor, method_count, "
                "global_percentile AS rank FROM results WHERE method = ? AND global_percentile >= ?",
                conn, params=(source, float(min_rank)))
    else:
        ev = review_queue._events()
        df = ev[ev["method_count"] >= int(min_methods)].assign(rank=float("nan"))[
            ["id", "event_key", "sender", "receiver", "ligand", "receptor", "method_count", "rank"] + methods]
    if cell_types:
        keep = set(cell_types)
        df = df[df["sender"].isin(keep) & df["receiver"].isin(keep)]
    return df


def network(source: str = "consensus", min_methods: int = 2, min_rank: float = 90.0,
            cell_types: list[str] | None = None) -> dict[str, Any]:
    df = _filtered(source, min_methods, min_rank, cell_types)
    all_types = sorted(set(index_store.facets()["senders"]) | set(index_store.facets()["receivers"]))
    types = sorted(cell_types) if cell_types else all_types

    edges = (df.groupby(["sender", "receiver"])["event_key"].nunique()
               .reset_index(name="count").sort_values("count", ascending=False))
    out_c = edges.groupby("sender")["count"].sum()
    in_c = edges.groupby("receiver")["count"].sum()
    nodes = [{"name": t, "out": int(out_c.get(t, 0)), "in": int(in_c.get(t, 0)),
              "self": int(edges[(edges["sender"] == t) & (edges["receiver"] == t)]["count"].sum())}
             for t in types]
    return {
        "source": source, "min_methods": int(min_methods), "min_rank": float(min_rank),
        "methods": index_store.methods(),
        "all_cell_types": all_types,
        "nodes": nodes,
        "edges": [{"sender": r.sender, "receiver": r.receiver, "count": int(r.count)} for r in edges.itertuples()],
        "n_interactions": int(df["event_key"].nunique()),
        "note": ("Each arrow counts candidate ligand–receptor interactions from the sender to the receiver "
                 "that pass the filters. It is a count of computational candidates, not a measure of how "
                 "strongly the populations communicate."),
    }


def pair(sender: str, receiver: str, source: str = "consensus", min_methods: int = 2,
         min_rank: float = 90.0, limit: int = 100) -> dict[str, Any]:
    df = _filtered(source, min_methods, min_rank, None)
    sub = df[(df["sender"] == sender) & (df["receiver"] == receiver)]
    if source in index_store.methods():
        sub = sub.sort_values("rank", ascending=False)
    else:
        sub = sub.sort_values(["method_count", "ligand", "receptor"], ascending=[False, True, True])
    methods = index_store.methods()
    rows = []
    for _, r in sub.head(limit).iterrows():
        row = {"result_id": int(r["id"]), "ligand": r["ligand"], "receptor": r["receptor"],
               "rank": None if pd.isna(r["rank"]) else round(float(r["rank"]), 1),
               "method_count": int(r["method_count"])}
        if source not in methods:
            row["ranks"] = {m: (None if pd.isna(r[m]) else round(float(r[m]), 1)) for m in methods}
        rows.append(row)
    return {"sender": sender, "receiver": receiver, "total": int(len(sub)), "rows": rows}
