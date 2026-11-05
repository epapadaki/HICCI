"""
Standardise SingleCellSignalR output for the platform.

Input : results/singlecellsignalr/singlecellsignalr_raw.csv
Output: results/singlecellsignalr/standardized_interactions.csv

The score is SingleCellSignalR's LRscore (0-1, higher is stronger). Within-tool
percentiles are computed over the interactions SingleCellSignalR reports
(LRscore >= its threshold), as for the other tools.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "results" / "singlecellsignalr" / "singlecellsignalr_raw.csv"
OUT = PROJECT_ROOT / "results" / "singlecellsignalr" / "standardized_interactions.csv"


def main() -> None:
    if not RAW.exists():
        raise SystemExit(f"Input not found: {RAW}")
    df = pd.read_csv(RAW)
    for col in ("sender_cell_type", "receiver_cell_type", "ligand", "receptor"):
        df[col] = df[col].astype(str).str.strip()
    df["ligand"] = df["ligand"].str.upper()
    df["receptor"] = df["receptor"].str.upper()
    df["score"] = pd.to_numeric(df["LRscore"], errors="coerce")
    df = df.dropna(subset=["score"])
    df = (df.sort_values("score", ascending=False)
            .drop_duplicates(["sender_cell_type", "receiver_cell_type", "ligand", "receptor"]))
    out = pd.DataFrame({
        "method": "SingleCellSignalR",
        "interaction_key": (df["sender_cell_type"] + "||" + df["receiver_cell_type"] + "||"
                            + df["ligand"] + "||" + df["receptor"]),
        "interacting_pair": df["ligand"] + "_" + df["receptor"],
        "interaction_type": "ligand-receptor",
        "sender_cell_type": df["sender_cell_type"],
        "receiver_cell_type": df["receiver_cell_type"],
        "ligand": df["ligand"],
        "receptor": df["receptor"],
        "score": df["score"],
        "score_name": "LRscore",
        "score_higher_is_better": True,
        "scsr_signalling_mode": df.get("mode"),
        "scsr_interaction_type": df.get("interaction_type"),
    })
    out["method_global_percentile"] = out["score"].rank(method="average", pct=True) * 100
    out["method_cell_pair_percentile"] = (out.groupby(["sender_cell_type", "receiver_cell_type"])["score"]
                                            .rank(method="average", pct=True) * 100)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"Saved {len(out):,} SingleCellSignalR interactions to {OUT}")
    print(f"Cell pairs: {out[['sender_cell_type', 'receiver_cell_type']].drop_duplicates().shape[0]}, "
          f"L-R pairs: {out[['ligand', 'receptor']].drop_duplicates().shape[0]}")


if __name__ == "__main__":
    main()
