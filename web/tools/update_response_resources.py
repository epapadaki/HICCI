"""
Download, once, the two prior-knowledge resources used to check the receiver's
response to an interaction, and keep them as local snapshots:

  PROGENy   pathway-responsive genes (14 pathways, human, top 500 genes each)
  CollecTRI transcription factor -> target genes (signed), human

    data/biological_context/response/progeny_human_top500.tsv   (source, target, weight, p_value)
    data/biological_context/response/collectri_human.tsv        (source, target, weight)

Both come from OmniPath (Saez-Rodriguez lab), as used by decoupler.

    python web/tools/update_response_resources.py
"""

from __future__ import annotations

import gzip
import io
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT = PROJECT_ROOT / "data" / "biological_context" / "response"
UA = {"User-Agent": "Mozilla/5.0 (compatible; CCC-Platform)"}

PROGENY_URLS = [
    "https://omnipathdb.org/annotations?resources=PROGENy&format=tsv",
    "https://static.omnipathdb.org/resources/annotations_PROGENy_9606.tsv.gz",
    "http://no-tls.static.omnipathdb.org/resources/annotations_PROGENy_9606.tsv.gz",
]
COLLECTRI_URLS = [
    "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&organisms=9606&loops=yes&format=tsv",
    "https://static.omnipathdb.org/resources/interactions_collectri_9606.tsv.gz",
    "http://no-tls.static.omnipathdb.org/resources/interactions_collectri_9606.tsv.gz",
]
AP1 = {"JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"}
NFKB = {"NFKB1", "NFKB2", "REL", "RELA", "RELB"}


def fetch(urls: list[str], need: set[str], min_rows: int) -> pd.DataFrame:
    """First source whose table has the expected columns and size."""
    last = None
    for url in urls:
        short = url.split("?")[0]
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
                raw = r.read()
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            df = pd.read_csv(io.BytesIO(raw), sep="\t", low_memory=False)
            cols = set(df.columns)
            if not need <= cols or len(df) < min_rows:
                preview = raw[:300].decode("utf-8", errors="replace").replace("\n", " | ")
                print(f"  {short}: unexpected answer ({len(df)} rows, columns {sorted(cols)[:8]}). Starts with: {preview}")
                continue
            print(f"  downloaded {len(df):,} rows from {short}")
            return df
        except Exception as exc:          # try the next source
            last = exc
            print(f"  could not use {short}: {exc}")
    raise SystemExit("Download failed from every source" + (f" ({last})." if last else ".") +
                     " Please send these lines to the developer.")


def progeny(top: int = 500) -> pd.DataFrame:
    a = fetch(PROGENY_URLS, need={"label", "value", "record_id"}, min_rows=1000)
    sym = "genesymbol" if "genesymbol" in a.columns else "genesymbols"
    wide = a.pivot_table(index=["record_id", sym], columns="label", values="value", aggfunc="first").reset_index()
    wide = wide.rename(columns={sym: "target", "pathway": "source"})
    wide["weight"] = pd.to_numeric(wide["weight"], errors="coerce")
    wide["p_value"] = pd.to_numeric(wide["p_value"], errors="coerce")
    wide = wide.dropna(subset=["source", "target", "weight"])
    wide = (wide.sort_values(["source", "p_value"]).groupby("source").head(top))
    return wide[["source", "target", "weight", "p_value"]].reset_index(drop=True)


def collectri() -> pd.DataFrame:
    c = fetch(COLLECTRI_URLS, need={"source", "target"}, min_rows=1000)
    src = "source_genesymbol" if "source_genesymbol" in c.columns else "source"
    tgt = "target_genesymbol" if "target_genesymbol" in c.columns else "target"

    def tf_name(s: str) -> str:
        parts = {p for p in str(s).replace("COMPLEX:", "").replace("_", " ").split() if p}
        if parts & AP1 and parts <= AP1:
            return "AP1"
        if parts & NFKB and parts <= NFKB:
            return "NFKB"
        return str(s)

    stim = pd.to_numeric(c.get("is_stimulation", 1), errors="coerce").fillna(0)
    inh = pd.to_numeric(c.get("is_inhibition", 0), errors="coerce").fillna(0)
    out = pd.DataFrame({"source": c[src].map(tf_name), "target": c[tgt].astype(str),
                        "weight": [(-1 if (i > 0 and s == 0) else 1) for s, i in zip(stim, inh)]})
    out = out[~out["target"].str.contains("_|COMPLEX:", regex=True)]
    return out.drop_duplicates(["source", "target"]).reset_index(drop=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("PROGENy …")
    p = progeny()
    p.to_csv(OUT / "progeny_human_top500.tsv", sep="\t", index=False)
    print(f"Saved PROGENy: {p['source'].nunique()} pathways, {len(p):,} gene weights")
    print("CollecTRI …")
    c = collectri()
    c.to_csv(OUT / "collectri_human.tsv", sep="\t", index=False)
    print(f"Saved CollecTRI: {c['source'].nunique():,} transcription factors, {len(c):,} TF–target links")
    (OUT / "downloaded.txt").write_text(time.strftime("%Y-%m-%d %H:%M:%S") + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
