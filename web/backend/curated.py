"""
Curated cell–cell communication network.

Unlike the Network page (which summarises the raw candidate interactions), this
is built only from the interactions the expert kept in the Human-in-the-Loop
step: those ranked at or above the cutoff in their own tool, minus any the
expert removed, plus any the expert added explicitly.

The number of kept interactions between two populations is a count of curated
candidates, not a measurement of how strongly they communicate.
"""

from __future__ import annotations

from contextlib import closing
from typing import Any

import pandas as pd

from . import assessments, index_store, review_queue


def _kept(cutoff: float) -> pd.DataFrame:
    """Events kept after review: (best within-tool rank >= cutoff AND not removed)
    OR explicitly added by the expert."""
    ev = review_queue._events().copy()
    if ev.empty:
        return ev
    decisions = assessments.assessed_keys()        # event_key -> "Accept" / "Reject"
    accept = {k for k, d in decisions.items() if d == "Accept"}
    reject = {k for k, d in decisions.items() if d == "Reject"}
    best = ev["best_rank"].fillna(-1)
    keep = ((best >= float(cutoff)) & ~ev["event_key"].isin(reject)) | ev["event_key"].isin(accept)
    return ev[keep]


def _top(items: list[dict], key: str, n: int = 6) -> list[dict]:
    return sorted([x for x in items if x.get(key, 0) > 0], key=lambda x: x[key], reverse=True)[:n]


def _cci_rows(ev: pd.DataFrame, limit: int = 200) -> list[dict[str, Any]]:
    """Format kept-event rows as CCIs for a drill-down list."""
    rows: list[dict[str, Any]] = []
    if ev is None or ev.empty:
        return rows
    ev2 = ev.sort_values(["method_count", "ligand", "receptor"], ascending=[False, True, True])
    for _, r in ev2.head(limit).iterrows():
        rows.append({
            "result_id": int(r["id"]),
            "ligand": r["ligand"], "receptor": r["receptor"],
            "sender": r["sender"], "receiver": r["receiver"],
            "method_count": int(r["method_count"]),
        })
    return rows


def _pairs_of(ev: pd.DataFrame, limit: int = 12) -> list[dict[str, Any]]:
    if ev is None or ev.empty:
        return []
    p = (ev.groupby(["sender", "receiver"])["event_key"].nunique()
           .reset_index(name="count").sort_values("count", ascending=False))
    return [{"sender": r.sender, "receiver": r.receiver, "count": int(r.count)} for r in p.head(limit).itertuples()]


def _net_data(ev: pd.DataFrame, types: list[str]) -> dict[str, Any]:
    """Nodes (out/in/self degree) and edges for a set of events."""
    if ev is None or ev.empty:
        edges = pd.DataFrame(columns=["sender", "receiver", "count"])
    else:
        edges = (ev.groupby(["sender", "receiver"])["event_key"].nunique()
                   .reset_index(name="count").sort_values("count", ascending=False))
    out_c = edges.groupby("sender")["count"].sum() if not edges.empty else {}
    in_c = edges.groupby("receiver")["count"].sum() if not edges.empty else {}

    def self_count(t: str) -> int:
        if edges.empty:
            return 0
        m = edges[(edges["sender"] == t) & (edges["receiver"] == t)]
        return int(m["count"].sum())

    nodes = [{"name": t, "out": int(out_c.get(t, 0)), "in": int(in_c.get(t, 0)), "self": self_count(t)} for t in types]
    return {"nodes": nodes,
            "edges": [{"sender": r.sender, "receiver": r.receiver, "count": int(r.count)}
                      for r in edges.itertuples()] if not edges.empty else []}


def _metrics(ev: pd.DataFrame) -> dict[str, int]:
    if ev is None or ev.empty:
        return {"ccis": 0, "connections": 0, "populations": 0, "molecules": 0}
    return {
        "ccis": int(ev["event_key"].nunique()),
        "connections": int(ev.groupby(["sender", "receiver"]).ngroups),
        "populations": int(len(set(ev["sender"]) | set(ev["receiver"]))),
        "molecules": int(len(set(ev["ligand"].astype(str)) | set(ev["receptor"].astype(str)))),
    }


def hitl_impact(cutoff: float = 95.0) -> dict[str, Any]:
    """What expert review changed: all candidates (before) vs the kept set (after),
    plus an honest breakdown of Accepted / Rejected / Unreviewed."""
    allev = review_queue._events().copy()
    keptev = _kept(cutoff)
    all_types = sorted(set(index_store.facets()["senders"]) | set(index_store.facets()["receivers"]))

    decisions = assessments.assessed_keys()
    accept_keys = {k for k, d in decisions.items() if d == "Accept"}
    reject_keys = {k for k, d in decisions.items() if d == "Reject"}
    n_total = int(allev["event_key"].nunique()) if not allev.empty else 0
    kept_keys = set(keptev["event_key"]) if not keptev.empty else set()
    kept_explicit = len(kept_keys & accept_keys)

    before = _metrics(allev)
    after = _metrics(keptev)
    before_net = _net_data(allev, all_types)
    after_net = _net_data(keptev, all_types)

    # plain-language read-out of what the picture means
    active_after = [nd for nd in after_net["nodes"] if (nd["out"] + nd["in"]) > 0]
    hub = max(active_after, key=lambda nd: nd["out"]) if active_after else None
    if before["ccis"] == 0:
        takeaway = "No candidate interactions to review yet."
    elif after["ccis"] == 0:
        takeaway = (f"Before — {before['ccis']} candidate interactions. "
                    "After review nothing is kept with the current settings.")
    else:
        drop = round(100 * (1 - after["ccis"] / before["ccis"])) if before["ccis"] else 0
        takeaway = (f"Before — {before['ccis']} candidate interactions in a dense, largely all-to-all network. "
                    f"After expert review — {after['ccis']} kept ({drop}% removed as non-specific), a sparser, "
                    "structured network")
        if hub and hub["out"] > 0:
            takeaway += f" centred on {hub['name']} as the dominant signalling source."
        else:
            takeaway += "."

    return {
        "cutoff": float(cutoff),
        "before": before,
        "after": after,
        "before_net": before_net,
        "after_net": after_net,
        "takeaway": takeaway,
        "decisions": {
            "accepted": len(accept_keys),
            "rejected": len(reject_keys),
            "unreviewed": max(0, n_total - len(accept_keys) - len(reject_keys)),
            "total": n_total,
        },
        "kept_breakdown": {
            "explicit_accept": kept_explicit,
            "default_by_cutoff": len(kept_keys) - kept_explicit,
            "total_kept": len(kept_keys),
        },
        "note": ("Before = all candidate interactions. After = the kept set (rank ≥ cutoff in the tool, minus the ones "
                 "you removed, plus the ones you added). Most of the kept set may be kept by the cutoff rule rather than "
                 "individually reviewed — the breakdown shows how many were explicitly accepted, rejected, or left unreviewed."),
    }


def curated_network(cutoff: float = 95.0, cell_types: list[str] | None = None) -> dict[str, Any]:
    ev = _kept(cutoff)
    all_types = sorted(set(index_store.facets()["senders"]) | set(index_store.facets()["receivers"]))
    if cell_types:
        s = set(cell_types)
        ev = ev[ev["sender"].isin(s) & ev["receiver"].isin(s)]
        types = sorted(cell_types)
    else:
        types = all_types

    if ev.empty:
        edges = pd.DataFrame(columns=["sender", "receiver", "count"])
    else:
        edges = (ev.groupby(["sender", "receiver"])["event_key"].nunique()
                   .reset_index(name="count").sort_values("count", ascending=False))
    out_c = edges.groupby("sender")["count"].sum() if not edges.empty else {}
    in_c = edges.groupby("receiver")["count"].sum() if not edges.empty else {}

    def self_count(t: str) -> int:
        if edges.empty:
            return 0
        m = edges[(edges["sender"] == t) & (edges["receiver"] == t)]
        return int(m["count"].sum())

    nodes = [{"name": t, "out": int(out_c.get(t, 0)), "in": int(in_c.get(t, 0)), "self": self_count(t)}
             for t in types]

    summary = {
        "most_total": _top([{"name": n["name"], "value": n["out"] + n["in"]} for n in nodes], "value"),
        "top_senders": _top([{"name": n["name"], "value": n["out"]} for n in nodes], "value"),
        "top_receivers": _top([{"name": n["name"], "value": n["in"]} for n in nodes], "value"),
        "top_autocrine": _top([{"name": n["name"], "value": n["self"]} for n in nodes], "value"),
        "top_pairs": [{"sender": r.sender, "receiver": r.receiver, "count": int(r.count)}
                      for r in edges.head(8).itertuples()] if not edges.empty else [],
    }

    return {
        "cutoff": float(cutoff),
        "methods": index_store.methods(),
        "all_cell_types": all_types,
        "nodes": nodes,
        "edges": [{"sender": r.sender, "receiver": r.receiver, "count": int(r.count)}
                  for r in edges.itertuples()] if not edges.empty else [],
        "n_interactions": int(ev["event_key"].nunique()) if not ev.empty else 0,
        "n_populations": len([n for n in nodes if (n["out"] + n["in"]) > 0]),
        "summary": summary,
        "note": ("Built only from the interactions kept after expert review: those ranked at or above the "
                 "cutoff in their tool, minus the ones you removed, plus the ones you added. A count of kept "
                 "interactions is not a direct measure of signalling strength."),
    }


def lr_rewiring(cutoff: float = 95.0, max_pairs: int = 60) -> dict[str, Any]:
    """Kept ligand->receptor pairs with their CellChat signalling pathway, for the
    Ligand-receptor rewiring views: a ligand->receptor bipartite graph, and a
    ligand->receptor->pathway tripartite graph. Built only from the curated set."""
    ev = _kept(cutoff)
    empty = {"cutoff": float(cutoff), "pairs": [], "n_interactions": 0,
             "n_ligands": 0, "n_receptors": 0, "n_pathways": 0}
    if ev is None or ev.empty:
        return empty
    pw_map = _event_pathways()
    ev = ev.copy()
    ev["pathway"] = ev["event_key"].map(pw_map).fillna("other")
    pairs = []
    for (lig, rcp), grp in ev.groupby(["ligand", "receptor"]):
        pw = grp["pathway"].mode()
        pairs.append({
            "ligand": str(lig), "receptor": str(rcp),
            "pathway": str(pw.iat[0]) if len(pw) else "other",
            "count": int(grp["event_key"].nunique()),
        })
    pairs.sort(key=lambda x: (-x["count"], x["ligand"], x["receptor"]))
    pairs = pairs[:max_pairs]

    # Reactome links: receptor -> Reactome pathway. A receptor is annotated with
    # several Reactome pathways, so we keep the top pathways (by kept-interaction
    # count) and split each receptor's interaction count equally among the ones it
    # belongs to, so the ligand->receptor and receptor->pathway sides stay balanced.
    rec_total: dict[str, int] = {}
    for p in pairs:
        rec_total[p["receptor"]] = rec_total.get(p["receptor"], 0) + p["count"]
    rec_rpws: dict[str, list[str]] = {}
    pw_weight: dict[str, float] = {}
    for rc, tot in rec_total.items():
        rpws = _reactome_pathways(rc)
        rec_rpws[rc] = rpws
        for nm in rpws:
            pw_weight[nm] = pw_weight.get(nm, 0.0) + tot
    top_re = [nm for nm, _ in sorted(pw_weight.items(), key=lambda kv: -kv[1])[:12]]
    top_set = set(top_re)
    reactome_links: list[dict[str, Any]] = []
    for rc, tot in rec_total.items():
        rpws = [nm for nm in rec_rpws.get(rc, []) if nm in top_set]
        if not rpws:
            continue
        share = tot / len(rpws)
        for nm in rpws:
            reactome_links.append({"receptor": rc, "pathway": nm, "weight": round(share, 3)})

    return {
        "cutoff": float(cutoff),
        "pairs": pairs,
        "reactome_links": reactome_links,
        "reactome_pathways": top_re,
        "n_interactions": sum(p["count"] for p in pairs),
        "n_ligands": len({p["ligand"] for p in pairs}),
        "n_receptors": len({p["receptor"] for p in pairs}),
        "n_pathways": len(top_re),
        "note": ("Kept ligand-receptor pairs after expert review. Top: ligand -> receptor. "
                 "Bottom: the same axes extended to the Reactome pathways each receptor "
                 "participates in (ligand -> receptor -> pathway). Reactome annotation is "
                 "method-independent, so these views are the same whichever CCC tool produced "
                 "the predictions."),
    }


def lr_dotplot(cutoff: float = 95.0, max_pairs: int = 40) -> dict[str, Any]:
    """Curated interactions as a dot-plot grid: ligand->receptor pairs (rows) by
    sender->receiver cell pairs (columns), with within-tool rank and CellChat
    signalling pathway for colour."""
    ev = _kept(cutoff)
    empty = {"cutoff": float(cutoff), "cells": [], "lr_order": [], "cp_order": [], "n": 0}
    if ev is None or ev.empty:
        return empty
    pw_map = _event_pathways()
    ev = ev.copy()
    ev["pathway"] = ev["event_key"].map(pw_map).fillna("other")
    ev["lr"] = ev["ligand"].astype(str) + "→" + ev["receptor"].astype(str)
    ev["cp"] = ev["sender"].astype(str) + "→" + ev["receiver"].astype(str)
    counts = ev.groupby("lr")["event_key"].nunique().sort_values(ascending=False)
    top = list(counts.head(max_pairs).index)
    sub = ev[ev["lr"].isin(top)]
    cells = []
    for lr, g in sub.groupby("lr"):
        for cp, g2 in g.groupby("cp"):
            best = g2["best_rank"].max()
            cells.append({"lr": str(lr), "cp": str(cp),
                          "rank": float(best) if pd.notna(best) else 0.0,
                          "pathway": str(g2.iloc[0]["pathway"])})
    cps = sorted(sub["cp"].unique().tolist(),
                 key=lambda s: (s.split("→")[0], s.split("→")[-1]))
    return {
        "cutoff": float(cutoff), "cells": cells, "lr_order": top, "cp_order": cps,
        "n": int(sub["event_key"].nunique()),
        "note": ("Each dot is a kept ligand-receptor interaction from one cell population to another. "
                 "Dot size scales with the within-tool rank; colour marks the receptor, so the view is "
                 "the same whichever CCC tool produced the predictions."),
    }


def molecular_hubs(cutoff: float = 95.0, top: int = 15) -> dict[str, Any]:
    """Ligands and receptors that take part in the most kept interactions, with
    how many distinct cell populations each one links."""
    ev = _kept(cutoff)

    def hubs(col: str) -> list[dict[str, Any]]:
        if ev.empty:
            return []
        rows = []
        for name, grp in ev.groupby(col):
            senders = set(grp["sender"])
            receivers = set(grp["receiver"])
            rows.append({
                "name": str(name),
                "interactions": int(grp["event_key"].nunique()),
                "senders": len(senders),
                "receivers": len(receivers),
                "populations": len(senders | receivers),
            })
        rows.sort(key=lambda r: (r["interactions"], r["populations"]), reverse=True)
        return rows[:top]

    ligands = hubs("ligand")
    receptors = hubs("receptor")
    return {
        "cutoff": float(cutoff),
        "ligands": ligands,
        "receptors": receptors,
        "max_ligand": ligands[0]["interactions"] if ligands else 1,
        "max_receptor": receptors[0]["interactions"] if receptors else 1,
        "note": ("Molecules in the most kept interactions. “Populations” is how many distinct cell populations the "
                 "molecule links as sender or receiver. A high count marks a network hub — a candidate for further "
                 "study — not proof of functional importance."),
    }


def communication_profiles(cutoff: float = 95.0) -> dict[str, Any]:
    """One communication profile per cell population, from the kept set: how much
    it sends vs receives, its main partners, and the ligands / receptors it uses
    most. Derived only from the kept interactions — no biological-process label
    yet (that comes with the pathway aggregation)."""
    ev = _kept(cutoff)
    all_types = sorted(set(index_store.facets()["senders"]) | set(index_store.facets()["receivers"]))

    def top_mol(sub: pd.DataFrame, col: str, n: int = 3) -> list[dict[str, Any]]:
        if sub.empty:
            return []
        g = sub.groupby(col)["event_key"].nunique().sort_values(ascending=False)
        return [{"name": str(k), "count": int(v)} for k, v in g.head(n).items()]

    profiles = []
    if not ev.empty:
        for t in all_types:
            out_ev = ev[ev["sender"] == t]
            in_ev = ev[ev["receiver"] == t]
            outgoing = int(out_ev["event_key"].nunique())
            incoming = int(in_ev["event_key"].nunique())
            if outgoing == 0 and incoming == 0:
                continue
            partner_counts: dict[str, int] = {}
            for q, grp in out_ev.groupby("receiver"):
                if q == t:
                    continue
                partner_counts[q] = partner_counts.get(q, 0) + int(grp["event_key"].nunique())
            for q, grp in in_ev.groupby("sender"):
                if q == t:
                    continue
                partner_counts[q] = partner_counts.get(q, 0) + int(grp["event_key"].nunique())
            partners = sorted(({"name": str(k), "count": v} for k, v in partner_counts.items()),
                              key=lambda x: x["count"], reverse=True)[:5]
            self_ev = ev[(ev["sender"] == t) & (ev["receiver"] == t)]
            total = outgoing + incoming
            profiles.append({
                "name": t,
                "outgoing": outgoing,
                "incoming": incoming,
                "total": total,
                "autocrine": int(self_ev["event_key"].nunique()),
                "balance": round(outgoing / total, 3) if total else 0.5,
                "partners": partners,
                "top_ligands": top_mol(out_ev, "ligand"),
                "top_receptors": top_mol(in_ev, "receptor"),
            })
    profiles.sort(key=lambda p: p["total"], reverse=True)
    return {
        "cutoff": float(cutoff),
        "profiles": profiles,
        "max_total": profiles[0]["total"] if profiles else 1,
        "note": ("One profile per population, from the kept set. Outgoing = interactions where the population is the "
                 "sender (expressing the ligand); incoming = where it is the receiver (expressing the receptor). "
                 "Autocrine (within-population) interactions count on both sides. The most frequent ligand / receptor "
                 "are the molecules the population uses in the most kept interactions."),
    }


# --------------------------------------------------------------------------- #
#  Biological processes & pathways aggregation                                 #
#                                                                              #
#  Two independent, non-merged labelling sources:                             #
#   - CellChat signalling-pathway family (one label per interaction, from      #
#     CellChatDB; only for interactions CellChat reported).                     #
#   - Reactome pathways the receptor is annotated as a participant in          #
#     (gene-level; one receptor can belong to many pathways).                   #
#                                                                              #
#  Both are *annotations*, not evidence that a process is active in the data.  #
# --------------------------------------------------------------------------- #

_PW_SPLIT = ("_", "&", "+", "/", ":", "|")


def _event_pathways() -> dict[str, str]:
    """event_key -> CellChat signalling-pathway family (only where present)."""
    try:
        with closing(index_store._connect()) as conn:
            df = pd.read_sql_query("SELECT event_key, signaling_pathway FROM events", conn)
    except Exception:
        return {}
    out: dict[str, str] = {}
    for ek, pw in zip(df["event_key"], df["signaling_pathway"]):
        s = "" if pw is None else str(pw).strip()
        if s and s.lower() != "nan":
            out[ek] = s
    return out


_reactome_cache: dict[str, list[str]] = {}


def _reactome_subunit(gene: str) -> list[str]:
    g = gene.strip()
    if not g:
        return []
    if g in _reactome_cache:
        return _reactome_cache[g]
    names: list[str] = []
    try:
        from utils.reactome_local import get_receptor_pathways
        res = get_receptor_pathways(g) or {}
        seen: set[str] = set()
        for p in (res.get("pathways") or []):
            nm = (p.get("pathway_name") or "").strip()
            if nm and nm not in seen:
                seen.add(nm)
                names.append(nm)
    except Exception:
        names = []
    _reactome_cache[g] = names
    return names


def _split_receptor(receptor: str) -> list[str]:
    parts = [str(receptor)]
    for sep in _PW_SPLIT:
        parts = [q for p in parts for q in p.split(sep)]
    return [p.strip() for p in parts if p.strip()]


def _reactome_pathways(receptor: str) -> list[str]:
    """Union of Reactome pathway names across the receptor's subunits."""
    names: list[str] = []
    seen: set[str] = set()
    for sub in _split_receptor(receptor):
        for nm in _reactome_subunit(sub):
            if nm not in seen:
                seen.add(nm)
                names.append(nm)
    return names


def _aggregate_processes(ev: pd.DataFrame, label_fn, top: int = 20) -> dict[str, Any]:
    """Group kept events by a process label. `label_fn(row)` returns a list of
    labels for that interaction (possibly empty, or several for Reactome)."""
    proc: dict[str, dict[str, set]] = {}
    sp: dict[tuple, set] = {}
    pr: dict[tuple, set] = {}
    labeled_keys: set[str] = set()

    for r in ev.itertuples():
        labels = label_fn(r)
        if not labels:
            continue
        labeled_keys.add(r.event_key)
        for lb in labels:
            d = proc.setdefault(lb, {"ek": set(), "s": set(), "r": set()})
            d["ek"].add(r.event_key)
            d["s"].add(r.sender)
            d["r"].add(r.receiver)
            sp.setdefault((r.sender, lb), set()).add(r.event_key)
            pr.setdefault((lb, r.receiver), set()).add(r.event_key)

    processes = sorted(
        ({"name": k, "ccis": len(v["ek"]), "senders": len(v["s"]),
          "receivers": len(v["r"]), "populations": len(v["s"] | v["r"])}
         for k, v in proc.items()),
        key=lambda x: (x["ccis"], x["populations"]), reverse=True)

    top_names = [p["name"] for p in processes[:top]]
    top_set = set(top_names)
    # Sankey restricted to the top processes: Sender -> Process -> Receiver
    links_sp = [{"source": s, "target": p, "value": len(ks)}
                for (s, p), ks in sp.items() if p in top_set]
    links_pr = [{"source": p, "target": rr, "value": len(ks)}
                for (p, rr), ks in pr.items() if p in top_set]
    s_val: dict[str, int] = {}
    for l in links_sp:
        s_val[l["source"]] = s_val.get(l["source"], 0) + l["value"]
    r_val: dict[str, int] = {}
    for l in links_pr:
        r_val[l["target"]] = r_val.get(l["target"], 0) + l["value"]
    p_val = {p["name"]: p["ccis"] for p in processes if p["name"] in top_set}

    sankey = {
        "senders": sorted(({"name": k, "value": v} for k, v in s_val.items()),
                          key=lambda x: x["value"], reverse=True),
        "processes": [{"name": n, "value": p_val.get(n, 0)} for n in top_names],
        "receivers": sorted(({"name": k, "value": v} for k, v in r_val.items()),
                           key=lambda x: x["value"], reverse=True),
        "links_sp": sorted(links_sp, key=lambda x: x["value"], reverse=True),
        "links_pr": sorted(links_pr, key=lambda x: x["value"], reverse=True),
    }
    return {
        "processes": processes,
        "max": processes[0]["ccis"] if processes else 1,
        "n_processes": len(processes),
        "labeled": len(labeled_keys),
        "sankey": sankey,
    }


def biological_processes(cutoff: float = 95.0, top: int = 20) -> dict[str, Any]:
    ev = _kept(cutoff)
    total_kept = int(ev["event_key"].nunique()) if not ev.empty else 0

    pw = _event_pathways()
    cc = _aggregate_processes(ev, lambda r: [pw[r.event_key]] if r.event_key in pw else [], top=top)
    rc = _aggregate_processes(ev, lambda r: _reactome_pathways(r.receptor), top=top)

    cc.update({
        "total_kept": total_kept,
        "unlabeled": max(0, total_kept - cc["labeled"]),
        "note": ("CellChat signalling-pathway families, assigned to each interaction by CellChatDB. Available only for "
                 "interactions CellChat reported; the rest are unlabelled. A label is a database annotation, not evidence "
                 "that the pathway is active in your data. Each interaction carries one family, so counts do not overlap."),
    })
    rc.update({
        "total_kept": total_kept,
        "unlabeled": max(0, total_kept - rc["labeled"]),
        "note": ("Reactome pathways the receptor gene is annotated as a participant in. A receptor can belong to many "
                 "pathways, so one interaction is counted under each — totals overlap and are not additive. Participation "
                 "is a gene-level annotation (receptor only), not evidence that the pathway is active in your data."),
    })
    return {"cutoff": float(cutoff), "total_kept": total_kept, "cellchat": cc, "reactome": rc}


def biological_process_detail(source: str, name: str, cutoff: float = 95.0, limit: int = 200) -> dict[str, Any]:
    """The kept interactions carrying one process label, for drill-down."""
    ev = _kept(cutoff)
    if ev.empty:
        return {"source": source, "name": name, "total": 0, "rows": [], "populations": [], "pairs": []}

    if source == "cellchat":
        pw = _event_pathways()
        mask = ev["event_key"].map(lambda k: pw.get(k) == name)
    else:
        mask = ev["receptor"].map(lambda rc: name in _reactome_pathways(rc))
    sub = ev[mask].copy()

    pop = sorted(set(sub["sender"]) | set(sub["receiver"]))
    pairs = (sub.groupby(["sender", "receiver"])["event_key"].nunique()
                .reset_index(name="count").sort_values("count", ascending=False)) if not sub.empty else pd.DataFrame()
    rows = []
    if not sub.empty:
        sub = sub.sort_values(["method_count", "ligand", "receptor"], ascending=[False, True, True])
        for _, r in sub.head(limit).iterrows():
            rows.append({
                "result_id": int(r["id"]),
                "sender": r["sender"], "receiver": r["receiver"],
                "ligand": r["ligand"], "receptor": r["receptor"],
                "method_count": int(r["method_count"]),
            })
    return {
        "source": source, "name": name,
        "total": int(sub["event_key"].nunique()) if not sub.empty else 0,
        "populations": pop,
        "pairs": [{"sender": r.sender, "receiver": r.receiver, "count": int(r.count)}
                  for r in pairs.head(12).itertuples()] if not pairs.empty else [],
        "rows": rows,
    }


def molecular_hub_detail(kind: str, name: str, cutoff: float = 95.0, limit: int = 200) -> dict[str, Any]:
    """The kept interactions that use one ligand or receptor, for drill-down."""
    ev = _kept(cutoff)
    col = "ligand" if kind == "ligand" else "receptor"
    sub = ev[ev[col] == name].copy() if not ev.empty else ev
    pop = sorted(set(sub["sender"]) | set(sub["receiver"])) if (sub is not None and not sub.empty) else []
    return {
        "kind": kind, "name": name,
        "total": int(sub["event_key"].nunique()) if (sub is not None and not sub.empty) else 0,
        "populations": pop,
        "pairs": _pairs_of(sub),
        "rows": _cci_rows(sub, limit),
    }


def communication_profile_detail(population: str, cutoff: float = 95.0, limit: int = 150) -> dict[str, Any]:
    """The kept interactions a population sends and receives, for drill-down."""
    ev = _kept(cutoff)
    out_ev = ev[ev["sender"] == population].copy() if not ev.empty else ev
    in_ev = ev[ev["receiver"] == population].copy() if not ev.empty else ev
    return {
        "population": population,
        "outgoing": {
            "total": int(out_ev["event_key"].nunique()) if (out_ev is not None and not out_ev.empty) else 0,
            "rows": _cci_rows(out_ev, limit),
        },
        "incoming": {
            "total": int(in_ev["event_key"].nunique()) if (in_ev is not None and not in_ev.empty) else 0,
            "rows": _cci_rows(in_ev, limit),
        },
    }


def hitl_impact_detail(bucket: str, cutoff: float = 95.0, limit: int = 200) -> dict[str, Any]:
    """The CCIs behind an expert-decision bucket (accepted / rejected / unreviewed) or the kept set."""
    if bucket == "kept":
        sub = _kept(cutoff)
    else:
        allev = review_queue._events().copy()
        decisions = assessments.assessed_keys()
        accept = {k for k, d in decisions.items() if d == "Accept"}
        reject = {k for k, d in decisions.items() if d == "Reject"}
        if allev.empty:
            sub = allev
        elif bucket == "accepted":
            sub = allev[allev["event_key"].isin(accept)]
        elif bucket == "rejected":
            sub = allev[allev["event_key"].isin(reject)]
        else:  # unreviewed
            sub = allev[~allev["event_key"].isin(accept) & ~allev["event_key"].isin(reject)]
    pop = sorted(set(sub["sender"]) | set(sub["receiver"])) if (sub is not None and not sub.empty) else []
    return {
        "bucket": bucket,
        "total": int(sub["event_key"].nunique()) if (sub is not None and not sub.empty) else 0,
        "populations": pop,
        "pairs": _pairs_of(sub),
        "rows": _cci_rows(sub, limit),
    }


def curated_pair(sender: str, receiver: str, cutoff: float = 95.0, limit: int = 100) -> dict[str, Any]:
    ev = _kept(cutoff)
    if not ev.empty:
        ev = ev[(ev["sender"] == sender) & (ev["receiver"] == receiver)].copy()
    methods = index_store.methods()
    rows = []
    if not ev.empty:
        ev = ev.sort_values(["method_count", "ligand", "receptor"], ascending=[False, True, True])
        for _, r in ev.head(limit).iterrows():
            rows.append({
                "result_id": int(r["id"]),
                "ligand": r["ligand"], "receptor": r["receptor"],
                "method_count": int(r["method_count"]),
                "ranks": {m: (None if (m not in r or pd.isna(r[m])) else round(float(r[m]), 1)) for m in methods},
            })
    return {"sender": sender, "receiver": receiver, "total": len(rows), "rows": rows}
