"""
Biological Insights Summary — a short, automatically written read-out of the
curated network.

It is rule-based, not generative: every sentence is produced from the real
aggregated numbers and is only emitted when the data actually supports it.
Each insight carries a `goto` (and sometimes a `focus`) so the reader can jump
straight to the view and the interactions behind the claim. Pathway families are
reported as annotations ("associated with"), never as proven activation.
"""

from __future__ import annotations

import statistics
from typing import Any

from . import curated


def _join(names: list[str]) -> str:
    names = [str(n) for n in names]
    if not names:
        return ""
    if len(names) == 1:
        return names[0]
    return ", ".join(names[:-1]) + " and " + names[-1]


def biological_insights(cutoff: float = 95.0) -> dict[str, Any]:
    net = curated.curated_network(cutoff)
    n = net.get("n_interactions", 0)
    note = ("Written automatically from the retained results — descriptive counts over unique interactions, "
            "not statistical claims. Pathway families are database annotations (associated with), not evidence "
            "of activation. Click an insight to see the interactions behind it.")
    if not n:
        return {"cutoff": float(cutoff), "empty": True, "insights": [], "note": note}

    insights: list[dict[str, Any]] = []
    nodes = net["nodes"]
    pops = [d for d in nodes if (d["out"] + d["in"]) > 0]

    # 1) scale of the curated network
    insights.append({
        "text": f"The curated network retains {n} unique interactions among "
                f"{net['n_populations']} cell populations across {len(net['edges'])} sender–receiver pairs.",
        "goto": "network",
    })

    # 2) major receivers (clearly above the typical in-degree, and net receivers)
    if pops:
        med_in = statistics.median([d["in"] for d in pops]) or 0
        rec = [d for d in sorted(pops, key=lambda d: d["in"], reverse=True)
               if d["in"] > 0 and d["in"] >= max(1, 1.5 * med_in) and d["in"] >= d["out"]][:2]
        if rec:
            insights.append({
                "text": f"{_join([d['name'] for d in rec])} stand out as major receivers of signalling "
                        f"({rec[0]['in']} incoming interactions for {rec[0]['name']}).",
                "goto": "profiles", "focus": {"population": rec[0]["name"]},
            })

        # 3) prominent senders
        med_out = statistics.median([d["out"] for d in pops]) or 0
        snd = [d for d in sorted(pops, key=lambda d: d["out"], reverse=True)
               if d["out"] > 0 and d["out"] >= max(1, 1.5 * med_out) and d["out"] >= d["in"]][:2]
        if snd:
            insights.append({
                "text": f"{_join([d['name'] for d in snd])} are prominent signal senders "
                        f"({snd[0]['out']} outgoing interactions for {snd[0]['name']}).",
                "goto": "profiles", "focus": {"population": snd[0]["name"]},
            })

    # 4) autocrine
    auto = net["summary"].get("top_autocrine") or []
    if auto and auto[0]["value"] > 0:
        insights.append({
            "text": f"{auto[0]['name']} shows notable autocrine (within-population) signalling "
                    f"({auto[0]['value']} interactions).",
            "goto": "network",
        })

    # 5) molecular hubs
    hubs = curated.molecular_hubs(cutoff)
    if hubs["ligands"] and hubs["receptors"]:
        lg, rc = hubs["ligands"][0], hubs["receptors"][0]
        insights.append({
            "text": f"The busiest molecular hubs are ligand {lg['name']} ({lg['interactions']} interactions) "
                    f"and receptor {rc['name']} ({rc['interactions']}).",
            "goto": "hubs", "focus": {"kind": "ligand", "name": lg["name"]},
        })

    # 6) biological processes (CellChat only — cheap; Reactome is skipped here)
    ev = curated._kept(cutoff)
    total_kept = int(ev["event_key"].nunique()) if not ev.empty else 0
    pw = curated._event_pathways()
    cc = curated._aggregate_processes(ev, lambda r: [pw[r.event_key]] if r.event_key in pw else [], top=5)
    if cc["labeled"] > 0 and cc["processes"]:
        top3 = [p["name"] for p in cc["processes"][:3]]
        frac = round(cc["labeled"] / total_kept * 100) if total_kept else 0
        insights.append({
            "text": f"Among interactions with a CellChat annotation ({frac}% of the kept set), the most frequent "
                    f"pathway families are {_join(top3)}.",
            "goto": "bio", "focus": {"source": "cellchat", "name": top3[0]},
        })

    # 7) HITL impact
    imp = curated.hitl_impact(cutoff)
    b, a = imp["before"]["ccis"], imp["after"]["ccis"]
    if b > 0:
        pct = round(a / b * 100)
        acc, rej = imp["decisions"]["accepted"], imp["decisions"]["rejected"]
        txt = f"Expert review narrowed {b} candidate interactions down to {a} ({pct}% retained)"
        if acc or rej:
            txt += f"; {acc} were explicitly accepted and {rej} explicitly rejected"
        txt += "."
        insights.append({"text": txt, "goto": "impact"})

    return {"cutoff": float(cutoff), "empty": False, "insights": insights, "note": note}
