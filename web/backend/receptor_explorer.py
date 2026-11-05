"""
Receptor-centric view.

1. Top-ranked set: candidate events whose rank (mean of the within-method
   ranks of the methods that found them) is at or above a threshold.
2. Recurrence: for every receptor entity (complexes kept whole), how many
   events it takes part in, how many of them are top-ranked, and whether it is
   over-represented in the top set (fold enrichment, one-sided Fisher exact
   test, Benjamini-Hochberg FDR across receptors).
3. One receptor: its ligands, sender and receiver populations, and events.
4. Optional: over-representation of the ligand partners' functions with
   g:Profiler, against all ligands of the dataset as background.

Recurrence is exploratory: it depends on the ligand-receptor databases, the
number of populations and how broadly genes are expressed.
"""

from __future__ import annotations

import math
import threading
from typing import Any

import pandas as pd

from . import config, index_store, review_queue

_cache: dict[str, Any] = {"key": None, "table": None}
_lock = threading.Lock()


def _events() -> pd.DataFrame:
    ev = review_queue._events().copy()
    parts = ev["event_key"].str.split("||", regex=False)
    ev["ligand_key"] = parts.str[2]
    ev["receptor_key"] = parts.str[3]
    return ev


TOP_SETS = {
    "all": "Reported by all tools",
    "two": "Reported by at least two tools",
}


def _select(ev: pd.DataFrame, top: str) -> tuple[pd.DataFrame, pd.Series, str]:
    """Universe of events and the 'top set' within it. Tools keep their own criteria:
    either WHICH tools report an event, or the top of ONE tool's own list."""
    methods = index_store.methods()
    if top in ("all", "two"):
        need = len(methods) if top == "all" else 2
        return ev, ev["method_count"] >= need, TOP_SETS[top]
    tool, _, thr = str(top).partition(":")
    if tool not in methods:
        raise ValueError("Unknown top-set definition.")
    thr_v = float(thr or 95)
    uni = ev[ev[tool].notna()]
    return uni, uni[tool] >= thr_v, f"Top {100 - thr_v:g}% of {tool}'s own results"


# ---- statistics without scipy -------------------------------------------------
def _log_comb(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_greater(a: int, b: int, c: int, d: int) -> float:
    """One-sided Fisher exact test (over-representation) for [[a, b], [c, d]]."""
    n1, n2, m = a + b, c + d, a + c          # top size, non-top size, receptor total
    N = n1 + n2
    lo, hi = a, min(n1, m)
    if lo > hi:
        return 1.0
    denom = _log_comb(N, m)
    logs = [_log_comb(n1, k) + _log_comb(n2, m - k) - denom for k in range(lo, hi + 1)]
    top = max(logs)
    return float(min(1.0, math.exp(top) * sum(math.exp(x - top) for x in logs)))


def bh(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [0.0] * n
    running = 1.0
    for rank_i in range(n - 1, -1, -1):
        i = order[rank_i]
        running = min(running, pvalues[i] * n / (rank_i + 1))
        q[i] = running
    return q


# ---- 1 + 2 --------------------------------------------------------------------
def recurrence(top: str = "all", min_events: int = 5) -> dict[str, Any]:
    key = (str(index_store._source_signature()), top, min_events)
    with _lock:
        if _cache["key"] == key:
            return _cache["table"]
    uni, top_mask, label = _select(_events(), top)
    uni = uni.assign(top=top_mask)
    total, total_top = len(uni), int(uni["top"].sum())
    g = uni.groupby("receptor_key")
    tab = pd.DataFrame({
        "receptor": g["receptor"].first(),
        "events": g.size(),
        "top": g["top"].sum().astype(int),
        "ligands": g["ligand_key"].nunique(),
        "senders": g["sender"].nunique(),
        "receivers": g["receiver"].nunique(),
        "max_support": g["method_count"].max().astype(int),
    }).reset_index()
    tab = tab[tab["events"] >= min_events].copy()
    bg = tab["events"] / total if total else 0
    tab["top_share"] = tab["top"] / total_top if total_top else 0
    tab["fold"] = (tab["top_share"] / bg).where(bg > 0)
    pvals = [fisher_greater(int(r.top), int(total_top - r.top), int(r.events - r.top),
                            int(total - total_top - (r.events - r.top))) for r in tab.itertuples()]
    tab["p"] = pvals
    tab["fdr"] = bh(pvals) if pvals else []
    tab = tab.sort_values(["fdr", "fold"], ascending=[True, False])
    rows = [{
        "receptor_key": r.receptor_key, "receptor": r.receptor,
        "complex": r.receptor_key.startswith("COMPLEX:") and "+" in r.receptor_key,
        "events": int(r.events), "top": int(r.top), "ligands": int(r.ligands),
        "senders": int(r.senders), "receivers": int(r.receivers), "max_support": int(r.max_support),
        "fold": None if pd.isna(r.fold) else round(float(r.fold), 2),
        "p": float(r.p), "fdr": float(r.fdr),
    } for r in tab.itertuples()]
    out = {
        "top": top, "top_label": label, "min_events": min_events,
        "methods": index_store.methods(),
        "events": int(total), "top_events": total_top,
        "receptors": rows,
        "definitions": {
            "fold": "Share of the top set that involves the receptor, divided by its share of all interactions considered. Above 1 means over-represented.",
            "fdr": "One-sided Fisher exact test for over-representation, corrected for testing all receptors (Benjamini–Hochberg).",
        },
        "caution": ("Recurrence is exploratory. It depends on how many ligands each database lists for a "
                    "receptor, on the number of populations and on how broadly the genes are expressed. It is "
                    "not evidence that a receptor is biologically more important. Each tool keeps its own "
                    "criteria; ranks of different tools are never combined."),
    }
    with _lock:
        _cache.update(key=key, table=out)
    return out


# ---- 3 ------------------------------------------------------------------------
def receptor(receptor_key: str, top: str = "all") -> dict[str, Any] | None:
    uni, top_mask, label = _select(_events(), top)
    uni = uni.assign(top=top_mask)
    sub = uni[uni["receptor_key"] == receptor_key]
    if sub.empty:
        return None
    lig = (sub.groupby("ligand_key").agg(ligand=("ligand", "first"), events=("event_key", "size"),
                                          top=("top", "sum"), tools=("method_count", "max"),
                                          senders=("sender", "nunique"))
              .reset_index().sort_values(["top", "tools", "ligand"], ascending=[False, False, True]))
    def pops(col):
        t = (sub.groupby(col).agg(events=("event_key", "size"), top=("top", "sum"), tools=("method_count", "max"))
                .reset_index().sort_values(["top", "events"], ascending=False))
        return [{"population": r[col], "events": int(r.events), "top": int(r.top), "tools": int(r.tools)}
                for _, r in t.iterrows()]
    methods = index_store.methods()
    top_rows = sub[sub["top"]].sort_values(["method_count", "ligand", "sender"], ascending=[False, True, True]).head(80)
    return {
        "receptor_key": receptor_key, "receptor": sub["receptor"].iloc[0], "top_label": label,
        "subunits": index_store.entity_display(receptor_key)["subunits"],
        "events": int(len(sub)), "top": int(sub["top"].sum()),
        "ligands": [{"ligand_key": r.ligand_key, "ligand": r.ligand, "events": int(r.events), "top": int(r.top),
                     "tools": int(r.tools), "senders": int(r.senders)} for r in lig.itertuples()],
        "senders": pops("sender"), "receivers": pops("receiver"),
        "top_events": [{"result_id": int(r["id"]), "sender": r["sender"], "receiver": r["receiver"],
                        "ligand": r["ligand"], "method_count": int(r["method_count"]),
                        "ranks": {m: (None if pd.isna(r[m]) else round(float(r[m]), 1)) for m in methods}}
                       for _, r in top_rows.iterrows()],
    }


# ---- 4 ------------------------------------------------------------------------
def ligand_enrichment(receptor_key: str, top: str = "all") -> dict[str, Any]:
    from utils.pathway_enrichment import run_gprofiler_enrichment
    ev, top_mask, _label = _select(_events(), top)
    sub = ev[top_mask & (ev["receptor_key"] == receptor_key)]
    def genes(keys):
        out = set()
        for k in keys:
            out.update(index_store.entity_display(k)["subunits"])
        return sorted(out)
    query = genes(sub["ligand_key"].unique())
    background = genes(ev["ligand_key"].unique())
    if len(query) < 3:
        raise ValueError("At least three ligand genes are needed; this receptor has "
                         f"{len(query)} in the selected set.")
    res = run_gprofiler_enrichment(query_genes=query, background_genes=background)
    terms = [{"id": t.get("native"), "name": t.get("name"), "source": t.get("source"),
              "p": t.get("p_value"), "overlap": t.get("intersection_size"), "term_size": t.get("term_size")}
             for t in (res.get("results") or [])[:40]]
    return {"query": query, "background_size": len(background), "terms": terms, "cached": res.get("cached"),
            "note": ("Over-representation among the ligand genes of this receptor, compared with all ligand genes "
                     "of the dataset (g:Profiler, g:SCS/FDR-corrected). It describes the ligands' annotated "
                     "functions, not the receptor's signalling in your cells.")}
