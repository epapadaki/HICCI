"""
Download Reactome's reaction-level mapping (UniProt -> reaction) for human and
store it as a local SQLite snapshot:

    data/biological_context/reactome/reactome_reactions.sqlite

It lets the platform find the curated Reactome reaction in which a ligand and
its receptor subunits take part (for example the binding reaction), and show
its diagram. Run once (internet needed), and again to update the snapshot:

    python web/tools/update_reactome_reactions.py
    python web/tools/update_reactome_reactions.py --file UniProt2ReactomeReactions.txt   (offline)
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
import time
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT = PROJECT_ROOT / "data" / "biological_context" / "reactome" / "reactome_reactions.sqlite"
URL = "https://reactome.org/download/current/UniProt2ReactomeReactions.txt"
SPECIES = "Homo sapiens"


def build(lines, source: str) -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix(".building")
    if tmp.exists():
        tmp.unlink()
    conn = sqlite3.connect(tmp)
    try:
        conn.execute("CREATE TABLE uniprot_reactions (uniprot_id TEXT, reaction_id TEXT, reaction_name TEXT, evidence TEXT)")
        conn.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT)")
        rows, n = [], 0
        for line in lines:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or parts[5] != SPECIES:
                continue
            rows.append((parts[0].split("-")[0], parts[1], parts[3], parts[4]))
            if len(rows) >= 50000:
                conn.executemany("INSERT INTO uniprot_reactions VALUES (?,?,?,?)", rows); n += len(rows); rows = []
        if rows:
            conn.executemany("INSERT INTO uniprot_reactions VALUES (?,?,?,?)", rows); n += len(rows)
        conn.execute("CREATE INDEX idx_up ON uniprot_reactions(uniprot_id)")
        conn.execute("CREATE INDEX idx_rx ON uniprot_reactions(reaction_id)")
        for k, v in (("source", source), ("species", SPECIES), ("rows", str(n)),
                     ("created", time.strftime("%Y-%m-%d %H:%M:%S"))):
            conn.execute("INSERT INTO metadata VALUES (?, ?)", (k, v))
        conn.commit()
    finally:
        conn.close()
    tmp.replace(OUT)
    return n


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", help="Use a downloaded UniProt2ReactomeReactions.txt instead of downloading")
    a = ap.parse_args(argv)
    if a.file:
        with open(a.file, encoding="utf-8", errors="replace") as fh:
            n = build(fh, str(Path(a.file).resolve()))
    else:
        print(f"Downloading {URL} …", flush=True)
        # Reactome refuses Python's default user agent (HTTP 403), so identify as a normal client.
        req = urllib.request.Request(URL, headers={
            "User-Agent": "Mozilla/5.0 (compatible; CCC-Platform reactome snapshot)",
            "Accept": "text/plain,*/*",
        })
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                text = resp.read().decode("utf-8", errors="replace")
        except OSError as exc:
            sys.exit(f"Download failed: {exc}\n"
                     f"Download the file in a browser from {URL} and run again with --file <path>.")
        n = build(text.splitlines(), URL)
    print(f"Saved {n:,} human UniProt–reaction links to {OUT}")


if __name__ == "__main__":
    main()
