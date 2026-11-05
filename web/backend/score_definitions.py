"""
What each method's raw score means.

Rule: if a score is not listed here, the platform says its definition is not
available. It never converts a raw score into a probability or a confidence.
"""

from __future__ import annotations

import csv
from functools import lru_cache
from typing import Any

from pathlib import Path

from . import datasets

_DEFINITIONS: dict[tuple[str, str], dict[str, str]] = {
    ("CellPhoneDB", "mean interaction expression"): {
        "label": "Mean interaction expression",
        "definition": (
            "Mean of the average ligand expression in the sender population and the "
            "average receptor expression in the receiver population, as reported by "
            "CellPhoneDB. For multi-subunit partners CellPhoneDB uses the subunit with "
            "the lowest average expression. It is an expression summary, not a "
            "probability of interaction."
        ),
    },
    ("NATMI", "edge_average_expression_weight"): {
        "label": "Edge expression weight (mean)",
        "definition": (
            "Product of the mean ligand expression in the sender population and the mean "
            "receptor expression in the receiver population, as computed by NATMI. It is "
            "an expression summary, not a probability of interaction."
        ),
    },
    ("CellChat", "communication probability"): {
        "label": "Communication probability",
        "definition": (
            "CellChat's interaction strength, computed from ligand, receptor and "
            "co-factor expression with a mass-action model. Despite the name it is a "
            "relative strength within this dataset, not the probability that the "
            "interaction occurs."
        ),
    },
    ("CellChat", "CellChat communication probability"): None,  # alias, filled below
    ("ICELLNET", "ICELLNET score"): {
        "label": "ICELLNET score",
        "definition": (
            "ICELLNET's product of the ligand expression in the sender and the receptor "
            "expression in the receiver, after each gene is rescaled so that its highest "
            "average expression across populations is 10; it lies between 0 and 100. Genes "
            "expressed in fewer than 10% of a population's cells count as 0. Complexes follow "
            "ICELLNET's rules for multiple subunits. A high score means high relative "
            "expression, not necessarily high absolute expression."
        ),
    },
    ("SingleCellSignalR", "LRscore"): {
        "label": "LRscore",
        "definition": (
            "SingleCellSignalR's regularised product of the mean ligand expression in the "
            "sender and the mean receptor expression in the receiver, scaled by the mean "
            "expression of the whole matrix; it lies between 0 and 1. SingleCellSignalR "
            "reports pairs at or above 0.5. It is an expression-based score, not a "
            "probability of interaction."
        ),
    },
}
_DEFINITIONS[("CellChat", "CellChat communication probability")] = _DEFINITIONS[
    ("CellChat", "communication probability")
]

PERCENTILE_DEFINITION = (
    "Rank of this raw score within the tool's own results in this dataset (100 = highest). "
    "Each tool uses its own criteria, so ranks are not compared between tools. It is a rank, "
    "not a probability or a confidence level."
)


def _cellchat_run_info() -> dict[str, str]:
    path = datasets.active().cellchat_run_info
    if not path or not path.exists():
        return {}
    return _read_run_info(str(path), path.stat().st_mtime)


@lru_cache(maxsize=16)
def _read_run_info(path: str, _mtime: float) -> dict[str, str]:
    with open(path, newline="") as fh:
        return {row["field"]: row["value"] for row in csv.DictReader(fh)}


def run_parameters(method: str) -> list[str]:
    """Run settings recorded by the pipeline, when available."""
    if method != "CellChat":
        return []
    info = _cellchat_run_info()
    out = []
    if info.get("CellChat_version"):
        out.append(f"CellChat {info['CellChat_version']}")
    if info.get("database"):
        out.append(info["database"])
    if info.get("expression_average"):
        out.append(f"expression average: {info['expression_average']}")
    if info.get("nboot"):
        out.append(f"{info['nboot']} permutations")
    if info.get("pvalue_threshold"):
        out.append(f"kept if p < {info['pvalue_threshold']}")
    return out


def describe_score(method: str, score_name: str | None) -> dict[str, Any]:
    name = (score_name or "").strip()
    known = _DEFINITIONS.get((method, name))
    if known:
        return {
            "name": name,
            "label": known["label"],
            "definition": known["definition"],
            "known": True,
            "run_parameters": run_parameters(method),
        }
    return {
        "name": name or None,
        "label": name or "Score",
        "definition": None,
        "known": False,
        "run_parameters": run_parameters(method),
    }
