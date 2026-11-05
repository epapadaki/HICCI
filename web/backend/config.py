"""
Paths and settings for the web platform.

The web app lives in ``<project>/web`` and reuses the existing scientific
modules in ``<project>/utils`` and the data in ``<project>/data`` and
``<project>/results``. Nothing in those folders is modified.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

WEB_ROOT = Path(__file__).resolve().parent.parent
PROJECT_ROOT = Path(
    os.environ.get("CCC_PROJECT_ROOT", WEB_ROOT.parent)
).resolve()

# Make the existing ``utils`` package importable.
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

TEMPLATES_DIR = WEB_ROOT / "frontend" / "templates"
STATIC_DIR = WEB_ROOT / "frontend" / "static"

# Writable area for the web app (index database, project settings).
WORKSPACE_DIR = WEB_ROOT / "workspace"
WORKSPACE_DIR.mkdir(parents=True, exist_ok=True)

# Tools the platform knows: column prefix in the consensus file -> display name.
# A dataset uses those that were run for it (see index_store.methods()).
METHODS = {
    "cellphonedb": "CellPhoneDB",
    "natmi": "NATMI",
    "cellchat": "CellChat",
    "singlecellsignalr": "SingleCellSignalR",
    "icellnet": "ICELLNET",
}

# PubMed
PUBMED_MAX_RESULTS = int(os.environ.get("CCC_PUBMED_MAX_RESULTS", "10"))

DECISIONS = ("Prioritize", "Accept", "Reject", "Uncertain")
