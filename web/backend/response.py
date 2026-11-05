"""
Does the receiver show the response expected from this interaction?

Within the receiver population, cells in which the receptor is detected (all
subunits, for a complex) are compared with cells in which it is not, for:
  - pathway activity (PROGENy, 14 pathways),
  - transcription factor activity (CollecTRI),
  - expression of the target genes of those transcription factors.
Activities are univariate linear model (ULM) t-values per cell, as in decoupler:
for each cell, expression over all genes is regressed on the pathway's (or TF's)
signed target weights.

Pathways expected for the receptor come from its Reactome pathways (keyword
match to the PROGENy pathway names); transcription factors linked to the
receptor come from SIGNOR causal paths (up to three steps). Both are
receptor-general knowledge.

This shows functional consistency within one snapshot of cells, not causality:
receptor-positive cells may also differ in identity or state, and cells are not
independent observations.
"""

from __future__ import annotations

import math
import time
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import config, documentation, expression, index_store

from utils.signor_downstream_local import get_downstream_causal_neighborhood

RES = config.PROJECT_ROOT / "data" / "biological_context" / "response"
PROGENY = RES / "progeny_human_top500.tsv"
COLLECTRI = RES / "collectri_human.tsv"
MAX_CELLS = 1500          # per group; larger groups are sampled (fixed seed)
MIN_CELLS = 10
MIN_TARGETS = 10        # sources with fewer targets give unstable activities
LOW_POWER = 30          # fewer cells than this in a group: little power
TF_MIN_DETECTED = 0.10  # TFs detected in fewer receiver cells are hidden

# PROGENy pathway -> words that identify it in Reactome pathway names
PATHWAY_WORDS = {
    "JAK-STAT": ["jak", "stat", "interleukin", "interferon", "cytokine signaling"],
    "NFkB": ["nf-kb", "nfkb", "nf-kappab", "toll", "tlr", "interleukin-1", "il-1"],
    "TNFa": ["tnf", "tumor necrosis factor"],
    "MAPK": ["mapk", "erk", "raf", "ras ", "rho gtpase", "p38", "jnk"],
    "PI3K": ["pi3k", "akt", "pip3", "mtor"],
    "EGFR": ["egfr", "erbb", "egf "],
    "TGFb": ["tgf", "smad", "bmp", "activin"],
    "WNT": ["wnt", "beta-catenin", "frizzled"],
    "VEGF": ["vegf"],
    "Hypoxia": ["hypoxia", "hif"],
    "p53": ["p53", "tp53"],
    "Trail": ["trail", "tnfsf10", "death receptor"],
    "Estrogen": ["estrogen", "esr1"],
    "Androgen": ["androgen"],
}


@lru_cache(maxsize=1)
def _net(path_str: str, _stamp: float) -> pd.DataFrame:
    df = pd.read_csv(path_str, sep="\t")
    df["target"] = df["target"].astype(str).str.upper()
    return df


def _load(path: Path) -> pd.DataFrame | None:
    return _net(str(path), path.stat().st_mtime) if path.exists() else None


def _matrix(store: expression.ExpressionStore, rows: np.ndarray) -> np.ndarray:
    """Dense (cells x genes) normalised expression for the given cells."""
    pos = np.full(len(store.cell_type), -1, dtype=np.int64)
    pos[rows] = np.arange(len(rows))
    indptr = np.asarray(store.indptr)
    indices = np.asarray(store.indices)
    sel = pos[indices] >= 0
    gene_of = np.repeat(np.arange(len(indptr) - 1), np.diff(indptr))[sel]
    cell_of = pos[indices[sel]]
    vals = np.asarray(store.data, dtype=np.float32)[sel]
    if store.raw_counts:
        vals = np.log1p(vals / store.library[rows][cell_of] * 1e4).astype(np.float32)
    X = np.zeros((len(rows), len(indptr) - 1), dtype=np.float32)
    X[cell_of, gene_of] = vals
    return X


def _ulm(X: np.ndarray, genes: list[str], net: pd.DataFrame) -> tuple[np.ndarray, list[str], dict[str, int]]:
    """ULM t-values (cells x sources)."""
    gidx = {g.upper(): i for i, g in enumerate(genes)}
    net = net[net["target"].isin(gidx)]
    sizes = net.groupby("source").size()
    sources = sorted(sizes[sizes >= MIN_TARGETS].index)
    if not sources:
        return np.zeros((X.shape[0], 0)), [], {}
    net = net[net["source"].isin(sources)]
    union = sorted(set(net["target"]))
    uidx = {g: i for i, g in enumerate(union)}
    W = np.zeros((len(union), len(sources)), dtype=np.float64)
    sidx = {s: i for i, s in enumerate(sources)}
    for s, t, w in zip(net["source"], net["target"], net["weight"]):
        W[uidx[t], sidx[s]] = w
    n = X.shape[1]
    Xd = X.astype(np.float64)
    mean_x = Xd.mean(axis=1)
    var_x = Xd.var(axis=1, ddof=1)
    mean_w = W.sum(axis=0) / n
    var_w = ((W ** 2).sum(axis=0) - n * mean_w ** 2) / (n - 1)
    XW = Xd[:, [gidx[g] for g in union]] @ W
    cov = (XW - n * np.outer(mean_x, mean_w)) / (n - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = cov / np.sqrt(np.outer(var_x, var_w))
        r = np.clip(np.nan_to_num(r), -0.999999, 0.999999)
        t = r * np.sqrt((n - 2) / (1 - r ** 2))
    return t, sources, {s: int(sizes[s]) for s in sources}


def _mannwhitney(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sided Mann-Whitney U p-value, normal approximation with tie correction."""
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return float("nan")
    ranks = pd.Series(np.concatenate([a, b])).rank().to_numpy()
    u1 = ranks[:n1].sum() - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    _, counts = np.unique(np.concatenate([a, b]), return_counts=True)
    N = n1 + n2
    tie = (counts ** 3 - counts).sum() / (N * (N - 1)) if N > 1 else 0
    sigma = math.sqrt(n1 * n2 / 12 * ((N + 1) - tie))
    if sigma == 0:
        return 1.0
    z = (u1 - mu) / sigma
    return float(math.erfc(abs(z) / math.sqrt(2)))


def _bh(p: list[float]) -> list[float]:
    idx = sorted(range(len(p)), key=lambda i: (p[i] if p[i] == p[i] else 1.0))
    q, run = [1.0] * len(p), 1.0
    for rank in range(len(p) - 1, -1, -1):
        i = idx[rank]
        pv = p[i] if p[i] == p[i] else 1.0
        run = min(run, pv * len(p) / (rank + 1))
        q[i] = run
    return q


def _compare(act: np.ndarray, names: list[str], pos: np.ndarray, sizes: dict[str, int]) -> list[dict[str, Any]]:
    out, pvals = [], []
    for j, s in enumerate(names):
        a, b = act[pos, j], act[~pos, j]
        sd = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / max(len(a) + len(b) - 2, 1)) if len(a) > 1 and len(b) > 1 else 0.0
        diff = float(a.mean() - b.mean())
        p = _mannwhitney(a, b)
        pvals.append(p)
        out.append({"name": s, "mean_positive": float(a.mean()), "mean_negative": float(b.mean()),
                    "difference": diff, "effect_size": float(diff / sd) if sd > 0 else 0.0,
                    "p": p, "targets": sizes.get(s)})
    for o, q in zip(out, _bh(pvals)):
        o["fdr"] = q
    return out


def _expected_pathways(genes: list[str]) -> dict[str, list[str]]:
    hits: dict[str, list[str]] = {}
    for g in genes:
        for p in documentation._reactome_for(g).get("pathways", []):
            name = (p.get("name") or "")
            low = f" {name.lower()} "
            for prog, words in PATHWAY_WORDS.items():
                if any(w in low for w in words):
                    hits.setdefault(prog, [])
                    if name not in hits[prog] and len(hits[prog]) < 3:
                        hits[prog].append(name)
    return hits


def _linked_tfs(genes: list[str], tfs: set[str]) -> dict[str, str]:
    found: dict[str, str] = {}
    for g in genes:
        try:
            res = get_downstream_causal_neighborhood(g, max_depth=3)
        except Exception:
            continue
        for rel in res.get("relations") or []:
            t = str(rel.get("target_entity") or "").upper()
            name = t if t in tfs else ("AP1" if t in {"JUN", "FOS"} and "AP1" in tfs else
                                       "NFKB" if t in {"RELA", "NFKB1"} and "NFKB" in tfs else None)
            if name and name not in found:
                found[name] = str(rel.get("path_text") or "")
    return found


_CACHE: dict[tuple, dict[str, Any]] = {}


def response(result_id: int) -> dict[str, Any] | None:
    key = (result_id, str(index_store._source_signature()),
           PROGENY.stat().st_mtime if PROGENY.exists() else None, COLLECTRI.stat().st_mtime if COLLECTRI.exists() else None)
    if key in _CACHE:
        return _CACHE[key]
    out = _response(result_id)
    if out is not None and out.get("available"):
        if len(_CACHE) > 200:
            _CACHE.clear()
        _CACHE[key] = out
    return out


def _response(result_id: int) -> dict[str, Any] | None:
    rec = index_store.get_result(result_id)
    if not rec:
        return None
    started = time.time()
    prog, coll = _load(PROGENY), _load(COLLECTRI)
    if prog is None or coll is None:
        return {"available": False, "message": ("The pathway and transcription factor resources are not installed. "
                                                "Run “python web/tools/update_response_resources.py” once (internet needed).")}
    store = expression.store_for_active()
    if store is None:
        return {"available": False, "message": "The expression matrix of this dataset has not been prepared."}
    rcp = index_store.entity_display(rec["receptor_entity_key"], rec.get("m_receptor"))
    lig = index_store.entity_display(rec["ligand_entity_key"], rec.get("m_ligand"))
    receiver = rec["receiver"]
    rows = store.type_rows.get(receiver)
    if rows is None or len(rows) == 0:
        return {"available": False, "message": f"No cells of “{receiver}” in the prepared data."}
    subunits = [g for g in rcp["subunits"]] or [rcp["label"]]
    detected = np.ones(len(store.cell_type), dtype=bool)
    missing = []
    for g in subunits:
        st = store.gene_stats(g)
        if st is None:
            missing.append(g); continue
        detected &= st["detected_mask"]
    if missing:
        return {"available": False, "message": f"Not in the expression data: {', '.join(missing)}."}
    pos_rows, neg_rows = rows[detected[rows]], rows[~detected[rows]]
    n_pos, n_neg = len(pos_rows), len(neg_rows)
    base = {"receiver": receiver, "sender": rec["sender"], "ligand": lig["label"], "receptor": rcp["label"],
            "subunits": subunits, "n_positive": n_pos, "n_negative": n_neg}
    if n_pos < MIN_CELLS or n_neg < MIN_CELLS:
        return {"available": True, "testable": False, **base,
                "message": (f"Not enough cells to compare: {n_pos} with and {n_neg} without the receptor "
                            f"(at least {MIN_CELLS} of each are needed).")}
    rng = np.random.default_rng(0)
    sampled = False
    if n_pos > MAX_CELLS:
        pos_rows = np.sort(rng.choice(pos_rows, MAX_CELLS, replace=False)); sampled = True
    if n_neg > MAX_CELLS:
        neg_rows = np.sort(rng.choice(neg_rows, MAX_CELLS, replace=False)); sampled = True
    sub = np.concatenate([pos_rows, neg_rows])
    is_pos = np.zeros(len(sub), dtype=bool); is_pos[:len(pos_rows)] = True
    X = _matrix(store, sub)

    # Data quality: do the two groups differ in how many genes are detected per cell?
    ngenes = (X > 0).sum(axis=1)
    g_pos, g_neg = float(ngenes[is_pos].mean()), float(ngenes[~is_pos].mean())
    quality_warning = None
    if min(g_pos, g_neg) > 0 and max(g_pos, g_neg) / min(g_pos, g_neg) > 1.25:
        quality_warning = (f"Cells with the receptor detect on average {g_pos:.0f} genes, cells without it {g_neg:.0f}. "
                           "Part of any difference may reflect cell quality or size rather than signalling.")
    low_power = min(n_pos, n_neg) < LOW_POWER

    act_p, pnames, psizes = _ulm(X, store.genes, prog)
    pathways = _compare(act_p, pnames, is_pos, psizes)
    expected = _expected_pathways(subunits)
    for p in pathways:
        p["expected"] = p["name"] in expected
        p["expected_because"] = expected.get(p["name"], [])
    pathways.sort(key=lambda p: (not p["expected"], -abs(p["effect_size"])))

    act_t, tnames, tsizes = _ulm(X, store.genes, coll)
    tfs = _compare(act_t, tnames, is_pos, tsizes)
    # Is the TF itself expressed in the receiver? (AP1 / NFKB: any member)
    members = {"AP1": ["JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"], "NFKB": ["NFKB1", "NFKB2", "REL", "RELA", "RELB"]}
    for t in tfs:
        best = 0.0
        for g in members.get(t["name"], [t["name"]]):
            st = store.gene_stats(g)
            if st:
                best = max(best, st["per_type"].get(receiver, {}).get("detected", 0.0))
        t["tf_detected"] = best
    hidden = sum(1 for t in tfs if t["tf_detected"] < TF_MIN_DETECTED)
    tfs = [t for t in tfs if t["tf_detected"] >= TF_MIN_DETECTED]
    linked = _linked_tfs(subunits, set(tnames))
    for t in tfs:
        t["linked"] = t["name"] in linked
        t["path"] = linked.get(t["name"])
    tfs_linked = sorted([t for t in tfs if t["linked"]], key=lambda t: -abs(t["effect_size"]))
    tfs_top = sorted([t for t in tfs if not t["linked"]], key=lambda t: (t["fdr"], -abs(t["effect_size"])))[:12]

    # target genes of the linked TFs (or of the top TFs when none are linked)
    focus = [t["name"] for t in (tfs_linked[:3] or tfs_top[:3])]
    gidx = {g.upper(): i for i, g in enumerate(store.genes)}
    targets = []
    for tf in focus:
        tg = coll[(coll["source"] == tf) & (coll["target"].isin(gidx))]
        rows_t = []
        for g, w in zip(tg["target"], tg["weight"]):
            col = X[:, gidx[g]]
            a, b = col[is_pos], col[~is_pos]
            rows_t.append({"gene": g, "mode": "activates" if w > 0 else "represses",
                           "mean_positive": float(a.mean()), "mean_negative": float(b.mean()),
                           "difference": float(a.mean() - b.mean())})
        rows_t.sort(key=lambda r: -abs(r["difference"]))
        changed = [r for r in rows_t if abs(r["difference"]) > 1e-9]
        consistent = sum(1 for r in changed if (r["difference"] > 0) == (r["mode"] == "activates"))
        # The TF's direction (more or less active with the receptor) sets which way its targets should move.
        tf_dir = next((t["difference"] for t in tfs if t["name"] == tf), 0.0)
        if tf_dir < 0:
            consistent = len(changed) - consistent
        n = len(changed)
        z = (consistent - n / 2) / math.sqrt(n / 4) if n else 0.0
        p_chance = float(0.5 * math.erfc(z / math.sqrt(2))) if n else 1.0
        targets.append({"tf": tf, "n_targets": len(rows_t), "n_changed": n, "consistent": consistent,
                        "share": (consistent / n) if n else None, "p_vs_chance": p_chance,
                        "tf_direction": "more active with the receptor" if tf_dir > 0 else "less active with the receptor",
                        "top": rows_t[:12]})

    return {
        "available": True, "testable": True, **base,
        "sampled": sampled, "cells_used": {"positive": int(len(pos_rows)), "negative": int(len(neg_rows))},
        "low_power": low_power, "quality_warning": quality_warning,
        "genes_per_cell": {"positive": round(g_pos), "negative": round(g_neg)},
        "tfs_hidden": hidden, "tf_min_detected": TF_MIN_DETECTED, "min_targets": MIN_TARGETS,
        "expected_matches": expected,
        "pathways": pathways, "tfs_linked": tfs_linked, "tfs_top": tfs_top, "targets": targets,
        "expected_source": "Reactome pathways of the receptor, matched by name to the PROGENy pathways",
        "linked_source": "SIGNOR causal paths from the receptor (up to three steps)",
        "method": ("Activities are ULM t-values per cell (decoupler's univariate linear model) with PROGENy "
                   "(top 500 genes per pathway) and CollecTRI. Groups are compared with a Mann–Whitney test; "
                   "FDR by Benjamini–Hochberg. Effect size: difference of means divided by the pooled SD."),
        "caution": ("Consistency, not causality. Cells with the receptor may also differ in identity or state, "
                    "and single cells are not independent observations, so p-values are optimistic. "
                    "Use this as supporting evidence for your judgement."),
        "seconds": round(time.time() - started, 1),
    }
