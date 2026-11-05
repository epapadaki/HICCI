"""
Standardise ICELLNET output for the platform.

Input : results/icellnet/icellnet_raw.csv   (pair_name such as "IL2 / IL2RA + IL2RB + IL2RG")
Output: results/icellnet/standardized_interactions.csv

Ligands and receptors with several subunits become complex entities
("COMPLEX:A+B", subunits sorted), the same keys that the complex-aware
CellPhoneDB and CellChat results use, so the tools can be matched.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "results" / "icellnet" / "icellnet_raw.csv"
OUT = PROJECT_ROOT / "results" / "icellnet" / "standardized_interactions.csv"


def split_side(text: str) -> list[str]:
    parts = [p.strip().upper() for p in re.split(r"\s*\+\s*|\s*&\s*", str(text)) if p and p.strip()]
    return [p for p in parts if p not in ("NA", "")]


def entity(genes: list[str]) -> str:
    uniq = sorted(set(genes))
    return f"GENE:{uniq[0]}" if len(uniq) == 1 else "COMPLEX:" + "+".join(uniq)


def main() -> None:
    if not RAW.exists():
        raise SystemExit(f"Input not found: {RAW}")
    df = pd.read_csv(RAW)
    sides = df["pair_name"].astype(str).str.split(r"\s*/\s*", n=1, regex=True)
    df["lig_genes"] = sides.str[0].map(split_side)
    df["rec_genes"] = sides.str[1].map(lambda s: split_side(s) if isinstance(s, str) else [])
    bad = df["lig_genes"].map(len).eq(0) | df["rec_genes"].map(len).eq(0)
    if bad.any():
        print(f"Skipped {int(bad.sum())} rows whose pair name could not be read")
    df = df[~bad].copy()
    df["ligand_entity"] = df["lig_genes"].map(entity)
    df["receptor_entity"] = df["rec_genes"].map(entity)
    df["ligand"] = df["lig_genes"].map(lambda g: "+".join(sorted(set(g))))
    df["receptor"] = df["rec_genes"].map(lambda g: "+".join(sorted(set(g))))
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"])
    df = df[df["score"] > 0]
    df["complex_aware_event_key_candidate"] = (df["sender_cell_type"].astype(str) + "||" + df["receiver_cell_type"].astype(str)
                                               + "||" + df["ligand_entity"] + "||" + df["receptor_entity"])
    df = df.sort_values("score", ascending=False).drop_duplicates("complex_aware_event_key_candidate")
    out = pd.DataFrame({
        "method": "ICELLNET",
        "complex_aware_event_key_candidate": df["complex_aware_event_key_candidate"],
        "interacting_pair": df["pair_name"],
        "interaction_type": "ligand-receptor",
        "sender_cell_type": df["sender_cell_type"],
        "receiver_cell_type": df["receiver_cell_type"],
        "ligand": df["ligand"],
        "receptor": df["receptor"],
        "ligand_entity_key": df["ligand_entity"],
        "receptor_entity_key": df["receptor_entity"],
        "score": df["score"],
        "score_name": "ICELLNET score",
        "score_higher_is_better": True,
    })
    out["method_global_percentile"] = out["score"].rank(method="average", pct=True) * 100
    out["method_cell_pair_percentile"] = (out.groupby(["sender_cell_type", "receiver_cell_type"])["score"]
                                            .rank(method="average", pct=True) * 100)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    n_cx = int((out["ligand_entity_key"].str.startswith("COMPLEX:") | out["receptor_entity_key"].str.startswith("COMPLEX:")).sum())
    print(f"Saved {len(out):,} ICELLNET interactions ({n_cx:,} involving a complex) to {OUT}")


if __name__ == "__main__":
    main()
