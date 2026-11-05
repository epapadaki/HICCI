"""
Project-level experiment description (optional).

Stored in ``web/workspace/project_context.json``. The documentation works
without it; when filled in, the Relevance tab compares each publication with
the user's experiment.
"""

from __future__ import annotations

import json
from typing import Any

from . import datasets, score_definitions

FIELDS = (
    "organism",
    "tissue_or_organ",
    "disease_or_condition",
    "experimental_groups",
    "notes",
)


def _defaults() -> dict[str, Any]:
    ctx: dict[str, Any] = {k: "" for k in FIELDS}
    ctx["experimental_groups"] = []
    # Pre-fill from the dataset description, then from the recorded run settings.
    meta = datasets.active().meta
    species = meta.get("organism") or score_definitions._cellchat_run_info().get("species", "")
    prefilled = []
    if species:
        ctx["organism"] = str(species).capitalize()
        prefilled.append("organism")
    if meta.get("tissue"):
        ctx["tissue_or_organ"] = str(meta["tissue"])
        prefilled.append("tissue_or_organ")
    ctx["prefilled"] = prefilled
    return ctx


def load() -> dict[str, Any]:
    if datasets.active_or_none() is None:
        ctx = {k: "" for k in FIELDS}
        ctx["experimental_groups"] = []
        ctx["prefilled"] = []
        return ctx
    path = datasets.active().context_file
    if path.exists():
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            return {**_defaults(), **data, "prefilled": []}
        except (OSError, json.JSONDecodeError):
            pass
    return _defaults()


def save(data: dict[str, Any]) -> dict[str, Any]:
    clean: dict[str, Any] = {}
    for key in FIELDS:
        value = data.get(key, "")
        if key == "experimental_groups":
            if isinstance(value, str):
                value = [v.strip() for v in value.replace(";", ",").split(",")]
            clean[key] = [str(v).strip() for v in value or [] if str(v).strip()]
        else:
            clean[key] = str(value or "").strip()
    datasets.active().context_file.write_text(
        json.dumps(clean, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return {**clean, "prefilled": []}


def is_filled(ctx: dict[str, Any]) -> bool:
    return any(ctx.get(k) for k in FIELDS if k != "notes")
