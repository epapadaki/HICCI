"""
Method comparison: which interactions each tool reports.

Built on the event table of the review queue (one row per candidate
interaction, with each tool's rank within its own results). Only WHICH tools
report an interaction is compared; ranks of different tools are not compared,
because each tool uses its own criteria and database.
"""

from __future__ import annotations

import threading
from itertools import combinations
from typing import Any

import pandas as pd

from . import config, index_store, review_queue

_cache: dict[str, Any] = {"sig": None, "summary": None}
_lock = threading.Lock()


def _methods() -> list[str]:
    return index_store.methods()


def _combo_key(row: pd.Series, methods: list[str]) -> str:
    return "+".join(m for m in methods if pd.notna(row[m]))


def summary() -> dict[str, Any]:
    sig = index_store._source_signature()
    with _lock:
        if _cache["sig"] == sig:
            return _cache["summary"]
    ev = review_queue._events()
    methods = _methods()
    present = ev[methods].notna()
    n_methods = len(methods)

    # ---- combinations (UpSet-style)
    keys = present.apply(lambda r: "+".join(m for m, ok in zip(methods, r) if ok), axis=1)
    combos = []
    for k in range(n_methods, 0, -1):
        for combo in combinations(methods, k):
            key = "+".join(combo)
            combos.append({"key": key, "methods": list(combo), "count": int((keys == key).sum())})
    combos.sort(key=lambda c: (-c["count"]))

    # ---- support levels
    support = {int(k): int(v) for k, v in ev["method_count"].value_counts().sort_index().items()}

    # ---- per method coverage
    per_method = []
    for m in methods:
        sub = ev[present[m]]
        others = present.drop(columns=[m]).any(axis=1)[present[m]]
        per_method.append({
            "method": m,
            "events": int(len(sub)),
            "lr_pairs": int(sub[["ligand", "receptor"]].drop_duplicates().shape[0]),
            "cell_pairs": int(sub[["sender", "receiver"]].drop_duplicates().shape[0]),
            "shared_share": float(others.mean()) if len(sub) else 0.0,
            "only": int((~others).sum()),
        })

    # ---- pairwise overlap (which interactions two tools both report)
    pairs = []
    for a, b in combinations(methods, 2):
        both = present[a] & present[b]
        union = present[a] | present[b]
        pairs.append({
            "a": a, "b": b,
            "shared": int(both.sum()),
            "jaccard": float(both.sum() / union.sum()) if union.sum() else None,
            "share_of_a": float(both.sum() / present[a].sum()) if present[a].sum() else None,
            "share_of_b": float(both.sum() / present[b].sum()) if present[b].sum() else None,
        })

    out = {
        "methods": methods,
        "events": int(len(ev)),
        "support": support,
        "combinations": combos,
        "per_method": per_method,
        "pairs": pairs,
        "definitions": {
            "jaccard": "Interactions reported by both tools, divided by those reported by either.",
            "share": "Share of one tool's interactions that the other tool also reports.",
        },
        "caution": (
            "Each tool is shown with its own criteria: CellPhoneDB (non-statistical analysis) reports "
            "expressed pairs, CellChat reports statistically significant ones, NATMI reports pairs with "
            "expression weights. They also use different ligand-receptor databases, so an interaction "
            "reported by one tool only may simply be absent from the others' databases. Overlap describes "
            "agreement between tools, not accuracy."
        ),
    }
    with _lock:
        _cache.update(sig=sig, summary=out)
    return out


def events(combo: str, page: int = 1, page_size: int = 50) -> dict[str, Any]:
    methods = _methods()
    wanted = [m for m in combo.split("+") if m in methods]
    if not wanted:
        raise ValueError("Unknown combination.")
    ev = review_queue._events()
    present = ev[methods].notna()
    mask = pd.Series(True, index=ev.index)
    for m in methods:
        mask &= present[m] if m in wanted else ~present[m]
    sub = ev[mask].sort_values(["ligand", "receptor", "sender", "receiver"])
    total = len(sub)
    pages = max(1, -(-total // page_size))
    page = max(1, min(int(page), pages))
    part = sub.iloc[(page - 1) * page_size: page * page_size]
    rows = [{
        "result_id": int(r["id"]), "sender": r["sender"], "receiver": r["receiver"],
        "ligand": r["ligand"], "receptor": r["receptor"],
        "ranks": {m: (None if pd.isna(r[m]) else round(float(r[m]), 1)) for m in methods},
    } for _, r in part.iterrows()]
    return {"combo": "+".join(wanted), "total": int(total), "page": page, "pages": pages, "rows": rows}
