"""
Download the ICELLNET ligand-receptor database once and keep it as a local
snapshot, so that ICELLNET runs are reproducible and work offline:

    data/biological_context/icellnet/ICELLNETdb.tsv

    python web/tools/update_icellnet_db.py
"""

from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT = PROJECT_ROOT / "data" / "biological_context" / "icellnet" / "ICELLNETdb.tsv"
URL = "https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"


def main() -> None:
    req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0 (compatible; CCC-Platform)"})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read()
    except OSError as exc:
        sys.exit(f"Download failed: {exc}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(data)
    n = data.decode("utf-8", errors="replace").count("\n") - 1
    print(f"Saved the ICELLNET database ({n:,} interactions) to {OUT}")


if __name__ == "__main__":
    main()
