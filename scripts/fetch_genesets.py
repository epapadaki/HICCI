#!/usr/bin/env python3
"""
Download the gene-set libraries used by the Functional enrichment page.

Run this ONCE on a machine with internet:

    python scripts/fetch_genesets.py

It writes four GMT files into data/genesets/. After that the enrichment page
runs fully offline (no network at query time).

No internet on the server? On any machine with a browser, open each URL printed
below, save the text as the matching filename, and copy the four files into
data/genesets/ by hand.
"""

from __future__ import annotations

import sys
from pathlib import Path
from urllib.request import Request, urlopen

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "data" / "genesets"

ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="

# (output filename, Enrichr library name) — the same 15 libraries and versions
# the GeneCurator tool reports.
LIBRARIES = [
    ("GO_Biological_Process.gmt", "GO_Biological_Process_2023"),
    ("GO_Cellular_Component.gmt", "GO_Cellular_Component_2023"),
    ("GO_Molecular_Function.gmt", "GO_Molecular_Function_2023"),
    ("KEGG.gmt", "KEGG_2021_Human"),
    ("Reactome.gmt", "Reactome_Pathways_2024"),
    ("WikiPathways.gmt", "WikiPathways_2024_Human"),
    ("BioPlanet.gmt", "BioPlanet_2019"),
    ("MSigDB_Hallmark.gmt", "MSigDB_Hallmark_2020"),
    ("DisGeNET.gmt", "DisGeNET"),
    ("Disease_Perturbations_up.gmt", "Disease_Perturbations_from_GEO_up"),
    ("Disease_Perturbations_down.gmt", "Disease_Perturbations_from_GEO_down"),
    ("DSigDB.gmt", "DSigDB"),
    ("Drug_Perturbations_up.gmt", "Drug_Perturbations_from_GEO_up"),
    ("Drug_Perturbations_down.gmt", "Drug_Perturbations_from_GEO_down"),
]


def _download(lib_name: str) -> bytes:
    req = Request(ENRICHR + lib_name, headers={"User-Agent": "CCC-Platform/1.0"})
    with urlopen(req, timeout=120) as resp:
        return resp.read()


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ok, failed = 0, []
    for fname, lib_name in LIBRARIES:
        dest = OUT_DIR / fname
        try:
            data = _download(lib_name)
            if not data or b"\t" not in data:
                raise ValueError("empty or unexpected response")
            dest.write_bytes(data)
            n_terms = data.count(b"\n")
            print(f"  ok   {fname:<28} {n_terms} terms  ({dest})")
            ok += 1
        except Exception as e:  # noqa: BLE001
            failed.append((fname, lib_name, str(e)))
            print(f"  FAIL {fname:<28} {e}")

    print(f"\n{ok}/{len(LIBRARIES)} libraries written to {OUT_DIR}")
    if failed:
        print("\nCould not download these automatically. Open each URL in a browser,")
        print("save the page as the filename shown, and put it in data/genesets/:")
        for fname, lib_name, _ in failed:
            print(f"  {fname}")
            print(f"      {ENRICHR}{lib_name}")
        return 1
    print("Done. Reload the Functional enrichment page.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
