"""
Evidence-based aggregation over the kept set.

For every kept interaction we reuse the same per-interaction evidence the
documentation panel shows, and reduce each source to a plain yes / no / not-
assessable verdict:

  - literature  : at least one supporting publication was found.
  - context     : at least one of those publications is in the same organism
                  AND tissue as this dataset (same biological context).
  - expression  : the ligand is detected in the sender population AND the
                  receptor in the receiver population (>= EXPR_DETECT of cells).
  - signalling  : in the receiver, an expected pathway of the receptor, or a
                  TF causally linked to it, is significantly more active in the
                  cells that express the receptor (FDR < 0.05).
  - mechanism   : the pair has a known interaction type in CellChatDB.

The sources are kept separate and never combined into a single score. Each
verdict is cached per interaction (by the dataset's source signature), so the
first pass is slow but later passes are instant.
"""

from __future__ import annotations

from typing import Any

from . import curated, index_store

EXPR_DETECT = 0.10  # share of cells that must detect a gene to count it as expressed

SOURCES: list[tuple[str, str]] = [
    ("literature", "Literature"),
    ("context", "Same context"),
    ("expression", "Expression"),
    ("signalling", "Signalling response"),
    ("mechanism", "Known mechanism"),
]
SRC_KEYS = [k for k, _ in SOURCES]

_cache: dict[str, Any] = {"sig": None, "by_id": {}}


def _ensure_sig() -> None:
    sig = index_store._source_signature()
    if _cache["sig"] != sig:
        _cache["sig"] = sig
        _cache["by_id"] = {}


# --------------------------------------------------------------- per-source verdicts
def _best_det(side: dict[str, Any] | None) -> float | None:
    if not side:
        return None
    if side.get("joint"):
        return side["joint"].get("detected")
    found = [x["detected"] for x in (side.get("subunits") or []) if x.get("found")]
    return min(found) if found else None


def _v_expression(rid: int) -> bool | None:
    from . import expression
    ex = expression.dataset_evidence(rid)
    if not ex or not ex.get("available"):
        return None
    ld, rd = _best_det(ex.get("ligand")), _best_det(ex.get("receptor"))
    if ld is None or rd is None:
        return None
    return bool(ld >= EXPR_DETECT and rd >= EXPR_DETECT)


def _v_mechanism(rid: int) -> bool | None:
    from . import mechanism
    rec = index_store.get_result(rid)
    if not rec:
        return None
    lig = index_store.entity_display(rec["ligand_entity_key"], rec.get("m_ligand"))
    rcp = index_store.entity_display(rec["receptor_entity_key"], rec.get("m_receptor"))
    cls = mechanism._classify(lig, rcp)
    t = cls.get("type")
    return bool(t and t != "unknown")


def _v_literature(rid: int) -> tuple[bool | None, bool | None]:
    """Returns (literature, context)."""
    from . import literature
    lit = literature.literature(rid)
    if not lit:
        return None, None
    status = lit.get("status")
    if status == "unsupported":
        return False, False
    if status != "ok":
        return None, None
    pubs = lit.get("publications") or []
    has_lit = len(pubs) > 0
    has_ctx = False
    for p in pubs:
        items = ((p.get("relevance") or {}).get("items")) or []
        st = {it.get("aspect"): it.get("status") for it in items}
        if st.get("Organism") == "same" and st.get("Tissue") == "same":
            has_ctx = True
            break
    return has_lit, has_ctx


def _v_signalling(rid: int) -> bool | None:
    from . import response
    rs = response.response(rid)
    if not rs or not rs.get("available"):
        return None
    if not rs.get("testable"):
        return None
    pw = rs.get("pathways") or []
    tl = rs.get("tfs_linked") or []
    sig_pw = any(p.get("expected") and p.get("fdr") is not None and p["fdr"] < 0.05 for p in pw)
    sig_tf = any(t.get("fdr") is not None and t["fdr"] < 0.05 for t in tl)
    return bool(sig_pw or sig_tf)


def _evidence_for(rid: int) -> dict[str, bool | None]:
    out: dict[str, bool | None] = {k: None for k in SRC_KEYS}
    try:
        out["expression"] = _v_expression(rid)
    except Exception:
        pass
    try:
        out["mechanism"] = _v_mechanism(rid)
    except Exception:
        pass
    try:
        out["literature"], out["context"] = _v_literature(rid)
    except Exception:
        pass
    try:
        out["signalling"] = _v_signalling(rid)
    except Exception:
        pass
    return out


# --------------------------------------------------------------- aggregation
def evidence_aggregation(cutoff: float = 95.0, batch: int = 40, heatmap_rows: int = 80) -> dict[str, Any]:
    _ensure_sig()
    ev = curated._kept(cutoff)
    if ev.empty:
        return {"cutoff": float(cutoff), "total_kept": 0, "computed": 0, "remaining": 0,
                "sources": [{"key": k, "label": l, "supported": 0, "assessed": 0} for k, l in SOURCES],
                "combos": [], "heatmap": [], "expr_detect": EXPR_DETECT, "note": _NOTE}

    meta = {int(r.id): {"event_key": r.event_key, "ligand": r.ligand, "receptor": r.receptor,
                        "sender": r.sender, "receiver": r.receiver,
                        "method_count": int(getattr(r, "method_count", 0) or 0)}
            for r in ev.itertuples()}
    ids = list(meta.keys())

    todo = [i for i in ids if i not in _cache["by_id"]][:max(1, int(batch))]
    for i in todo:
        _cache["by_id"][i] = _evidence_for(i)

    computed = [i for i in ids if i in _cache["by_id"]]
    remaining = len(ids) - len(computed)

    # per-source counts
    srcs = []
    for k, label in SOURCES:
        supported = sum(1 for i in computed if _cache["by_id"][i].get(k) is True)
        assessed = sum(1 for i in computed if _cache["by_id"][i].get(k) is not None)
        srcs.append({"key": k, "label": label, "supported": supported, "assessed": assessed})

    # UpSet combinations (over the sources that are TRUE for each interaction)
    combo_count: dict[tuple, int] = {}
    for i in computed:
        d = _cache["by_id"][i]
        key = tuple(k for k in SRC_KEYS if d.get(k) is True)
        combo_count[key] = combo_count.get(key, 0) + 1
    combos = sorted(({"sources": list(k), "count": v} for k, v in combo_count.items()),
                    key=lambda x: x["count"], reverse=True)[:14]

    # heatmap: strongest-supported first, capped
    def nsupp(i: int) -> int:
        d = _cache["by_id"][i]
        return sum(1 for k in SRC_KEYS if d.get(k) is True)
    order = sorted(computed, key=lambda i: (nsupp(i), meta[i]["method_count"]), reverse=True)[:max(1, int(heatmap_rows))]
    heatmap = []
    for i in order:
        d = _cache["by_id"][i]
        heatmap.append({
            "result_id": i, "ligand": meta[i]["ligand"], "receptor": meta[i]["receptor"],
            "sender": meta[i]["sender"], "receiver": meta[i]["receiver"],
            "cells": {k: d.get(k) for k in SRC_KEYS},
        })

    return {
        "cutoff": float(cutoff),
        "total_kept": len(ids),
        "computed": len(computed),
        "remaining": remaining,
        "expr_detect": EXPR_DETECT,
        "sources": srcs,
        "combos": combos,
        "heatmap": heatmap,
        "note": _NOTE,
    }


def evidence_source_detail(source: str, cutoff: float = 95.0, limit: int = 200) -> dict[str, Any]:
    """The kept interactions that one evidence source supports, for drill-down."""
    _ensure_sig()
    ev = curated._kept(cutoff)
    if ev.empty:
        return {"source": source, "total": 0, "rows": []}
    meta = {int(r.id): {"result_id": int(r.id), "ligand": r.ligand, "receptor": r.receptor,
                        "sender": r.sender, "receiver": r.receiver,
                        "method_count": int(getattr(r, "method_count", 0) or 0)}
            for r in ev.itertuples()}
    rows = [meta[i] for i in meta
            if i in _cache["by_id"] and _cache["by_id"][i].get(source) is True]
    rows.sort(key=lambda r: r["method_count"], reverse=True)
    return {"source": source, "total": len(rows), "rows": rows[:limit]}


_NOTE = ("Each source is a separate yes/no verdict for the interaction — they are never merged into one score. "
         "“Assessed” is how many kept interactions could be evaluated for that source; the rest lack the data "
         "(no expression matrix, no literature hit to compare, signalling resources not installed, …). Literature "
         "and context are abstract-level; expression and signalling are consistency in your data, not proof of "
         "causality.")
