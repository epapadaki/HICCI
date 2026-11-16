"""
CCC Platform — FastAPI application.

Pages
  GET  /                                  interactions table + documentation panel

API
  GET  /api/status                        index state (building / ready / missing)
  POST /api/index/rebuild                 rebuild the results index
  GET  /api/facets                        methods, sender and receiver populations
  GET  /api/results                       one page of results (filters, sort, paging)
  GET  /api/results/{id}                  header: computational result
  GET  /api/results/{id}/biology          section A: molecules and receptor pathways
  GET  /api/results/{id}/literature       sections B and C: publications, relevance
  GET  /api/results/{id}/assessments      expert assessments of this event
  POST /api/results/{id}/assessments      save an expert assessment
  GET  /api/assessments/keys              assessed events (to mark table rows)
  GET  /api/assessments/export.csv        all assessments as CSV
  GET  /api/results/{id}/dataset-evidence expression in the user's dataset (section D)
  GET  /api/datasets                      datasets and the active one
  PUT  /api/datasets/active               switch dataset
  GET  /api/project-context               experiment description (optional)
  PUT  /api/project-context               save experiment description
"""

from __future__ import annotations

import shutil
import threading
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from . import advisor, assessments, response, config, dataset_admin, mechanism, sources, datasets, method_comparison, network, receptor_explorer, provenance, documentation, expression, index_store, jobs, literature, project_context, review_queue, curated, evidence, insights, functional_enrichment, discussion, settings

@asynccontextmanager
async def lifespan(_app: FastAPI):
    # Build the index in the background so the page opens immediately
    # and can show progress.
    datasets.migrate_legacy_workspace()
    # Start every launch with nothing loaded, so closing the window never leaves
    # a dataset "remembered".
    datasets.clear_active()
    # Demo / ephemeral mode: start each launch completely clean. Irreversible.
    if settings.ephemeral_mode():
        try:
            datasets.wipe_user_datasets()
            from database.db import wipe_all_feedback, init_database
            init_database()
            wipe_all_feedback()
            uploads = config.PROJECT_ROOT / "data" / "uploads"
            keep = {"melanoma_tme.h5ad"}  # the shipped use-case dataset
            if uploads.exists():
                for p in uploads.iterdir():
                    if p.name in keep:
                        continue
                    if p.is_dir():
                        shutil.rmtree(p, ignore_errors=True)
                    else:
                        p.unlink(missing_ok=True)
        except Exception:
            pass
    threading.Thread(target=documentation.warm_up, daemon=True).start()
    yield


app = FastAPI(title="HICCI", docs_url="/api/docs", redoc_url=None, lifespan=lifespan)
app.mount("/static", StaticFiles(directory=str(config.STATIC_DIR)), name="static")
templates = Jinja2Templates(directory=str(config.TEMPLATES_DIR))


def _asset_version() -> str:
    """Changes whenever a CSS/JS file changes, so browsers never use stale copies."""
    newest = max((f.stat().st_mtime for f in config.STATIC_DIR.rglob("*") if f.is_file()), default=0)
    return str(int(newest))


templates.env.globals["v"] = _asset_version()


@app.middleware("http")
async def _refresh_asset_version(request: Request, call_next):
    # Recompute on each request (cheap, local single-user app) so edited CSS/JS
    # is picked up without restarting, and browsers never serve stale copies.
    if request.url.path in ("/", "") or not request.url.path.startswith(("/static", "/api")):
        templates.env.globals["v"] = _asset_version()
    return await call_next(request)


def _require_index() -> None:
    if datasets.active_or_none() is None:
        raise HTTPException(status_code=409, detail="No dataset is loaded.")
    if not index_store.index_is_current():
        state = index_store.status()          # starts a rebuild when the index is out of date
        msg = ("Preparing the results of this dataset…" if state.get("state") == "building"
               else state.get("message") or "Index is not ready.")
        raise HTTPException(status_code=503, detail=msg)


# --------------------------------------------------------------------------- pages
@app.get("/", response_class=HTMLResponse)
def page_home(request: Request):
    return templates.TemplateResponse(request, "home.html", {"active": "home"})


@app.get("/interactions", response_class=HTMLResponse)
def page_interactions(request: Request):
    return templates.TemplateResponse(request, "interactions.html", {"active": "interactions"})


@app.get("/help", response_class=HTMLResponse)
def page_help(request: Request):
    return templates.TemplateResponse(request, "help.html", {"active": "help"})


@app.get("/api/sources")
def api_sources() -> list[dict[str, Any]]:
    return sources.sources()


@app.get("/network", response_class=HTMLResponse)
def page_network(request: Request):
    # The Candidate network page was removed; send any old link to Aggregated results.
    return RedirectResponse(url="/aggregated", status_code=307)


@app.get("/aggregated", response_class=HTMLResponse)
def page_aggregated(request: Request):
    return templates.TemplateResponse(request, "aggregated.html", {"active": "aggregated"})


@app.get("/enrichment", response_class=HTMLResponse)
def page_enrichment(request: Request):
    return templates.TemplateResponse(request, "enrichment.html", {"active": "enrichment"})


@app.get("/discussion", response_class=HTMLResponse)
def page_discussion(request: Request):
    return templates.TemplateResponse(request, "discussion.html", {"active": "discussion"})


@app.get("/methods", response_class=HTMLResponse)
def page_methods(request: Request):
    return templates.TemplateResponse(request, "methods.html", {"active": "methods"})


@app.get("/datasets", response_class=HTMLResponse)
def page_datasets(request: Request):
    return templates.TemplateResponse(request, "datasets.html", {"active": "datasets"})


@app.get("/run", response_class=HTMLResponse)
def page_run(request: Request):
    return templates.TemplateResponse(request, "run.html", {"active": "run"})


# --------------------------------------------------------------------------- index
@app.get("/api/status")
def api_status() -> dict[str, Any]:
    """Index state, plus the tools this dataset has results from."""
    st = dict(index_store.status())
    if st.get("state") == "ready":
        st["methods"] = index_store.methods()
    return st


@app.post("/api/index/rebuild")
def api_rebuild() -> dict[str, Any]:
    threading.Thread(target=index_store.build_index, kwargs={"force": True}, daemon=True).start()
    return {"state": "building"}


@app.get("/api/facets")
def api_facets() -> dict[str, Any]:
    _require_index()
    return index_store.facets()


@app.get("/api/results")
def api_results(
    method: list[str] = Query(default=[]),
    sender: Optional[str] = None,
    receiver: Optional[str] = None,
    q: Optional[str] = None,
    min_percentile: Optional[float] = None,
    min_methods: Optional[int] = None,
    sort: str = "percentile",
    page: int = 1,
    page_size: int = 50,
) -> dict[str, Any]:
    _require_index()
    return index_store.query_results(
        methods=method or None,
        sender=sender or None,
        receiver=receiver or None,
        q=q or None,
        min_percentile=min_percentile,
        min_methods=min_methods,
        sort=sort,
        page=page,
        page_size=page_size,
    )


@app.get("/api/interactions")
def api_interactions(
    method: list[str] = Query(default=[]),
    sender: Optional[str] = None,
    receiver: Optional[str] = None,
    q: Optional[str] = None,
    min_percentile: Optional[float] = None,
    min_methods: Optional[int] = None,
    sort: str = "tools",
    page: int = 1,
    page_size: int = 50,
) -> dict[str, Any]:
    """One row per interaction, with each tool's own rank."""
    _require_index()
    return index_store.query_events(methods=method or None, sender=sender or None, receiver=receiver or None,
                                    q=q or None, min_percentile=min_percentile, min_methods=min_methods,
                                    sort=sort, page=page, page_size=page_size)


# --------------------------------------------------------------------------- documentation
@app.get("/api/results/{result_id}")
def api_overview(result_id: int) -> dict[str, Any]:
    _require_index()
    data = documentation.overview(result_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Result not found.")
    return data


@app.get("/api/results/{result_id}/response")
def api_response(result_id: int) -> dict[str, Any]:
    _require_index()
    data = response.response(result_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Result not found.")
    return data


@app.get("/api/results/{result_id}/mechanism")
def api_mechanism(result_id: int) -> dict[str, Any]:
    _require_index()
    data = mechanism.mechanism(result_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Result not found.")
    return data


@app.get("/api/results/{result_id}/biology")
def api_biology(result_id: int) -> dict[str, Any]:
    _require_index()
    data = documentation.biology(result_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Result not found.")
    return data


@app.get("/api/results/{result_id}/literature")
def api_literature(result_id: int, refresh: bool = False) -> dict[str, Any]:
    _require_index()
    data = literature.literature(result_id, force_refresh=refresh)
    if data is None:
        raise HTTPException(status_code=404, detail="Result not found.")
    return data


# --------------------------------------------------------------------------- project
@app.get("/api/project-context")
def api_get_context() -> dict[str, Any]:
    return project_context.load()


@app.put("/api/project-context")
def api_put_context(payload: dict[str, Any]) -> dict[str, Any]:
    return project_context.save(payload)


# --------------------------------------------------------------------------- expert assessment
@app.get("/api/results/{result_id}/assessments")
def api_list_assessments(result_id: int) -> dict[str, Any]:
    _require_index()
    items = assessments.list_for_result(result_id)
    if items is None:
        raise HTTPException(status_code=404, detail="Result not found.")
    return {"items": items, "reasons": assessments.REASONS, "decisions": list(config.DECISIONS)}


@app.post("/api/results/{result_id}/assessments")
def api_save_assessment(result_id: int, payload: dict[str, Any]) -> dict[str, Any]:
    _require_index()
    try:
        return assessments.save(result_id, payload)
    except LookupError as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))


@app.post("/api/results/{result_id}/keep")
def api_set_keep(result_id: int, payload: dict[str, Any]) -> dict[str, Any]:
    """On/off HITL: keep (checked) or remove (unchecked) this interaction."""
    _require_index()
    try:
        return assessments.set_decision(result_id, bool(payload.get("kept", True)), payload.get("expert"))
    except LookupError as exc:
        raise HTTPException(status_code=404, detail=str(exc))


@app.post("/api/assessments/bulk")
def api_bulk_assess(payload: dict[str, Any]) -> dict[str, Any]:
    """Pathway-level curation: Accept or Reject EVERY interaction that matches the
    current filters (e.g. a whole signalling category), in one action."""
    _require_index()
    kept = bool(payload.get("kept", False))
    expert = payload.get("expert")
    f = payload.get("filters") or {}
    res = index_store.query_events(
        methods=(f.get("methods") or None),
        sender=(f.get("sender") or None),
        receiver=(f.get("receiver") or None),
        q=(f.get("q") or None),
        min_percentile=f.get("min_percentile"),
        min_methods=f.get("min_methods"),
        sort="tools", page=1, page_size=100000,
    )
    rows = res.get("rows", []) if isinstance(res, dict) else []
    # breakdown by signalling category, so the UI can show exactly what will change
    from collections import Counter as _Counter
    by_cat = _Counter(str(r.get("signaling_pathway") or r.get("category") or "?") for r in rows)
    if payload.get("dry_run"):
        return {"total": len(rows), "by_category": by_cat.most_common(), "kept": kept}
    done, errors = 0, 0
    for row in rows:
        rid = row.get("open_id", row.get("id"))
        if rid is None:
            errors += 1
            continue
        try:
            assessments.set_decision(int(rid), kept, expert)
            done += 1
        except Exception:
            errors += 1
    return {"count": done, "errors": errors, "total": len(rows), "by_category": by_cat.most_common(), "kept": kept}


@app.get("/api/pathways")
def api_pathways() -> dict[str, Any]:
    """Signalling pathways (CellChat categories) with counts, best within-tool rank
    and current curation status — for pathway-level HITL curation."""
    _require_index()
    res = index_store.query_events(sort="tools", page=1, page_size=100000)
    rows = res.get("rows", []) if isinstance(res, dict) else []
    assessed = assessments.assessed_keys()
    from collections import defaultdict as _dd
    agg: dict[str, dict] = _dd(lambda: {"count": 0, "best_rank": 0.0, "kept": 0, "rejected": 0})
    for r in rows:
        cat = str(r.get("signaling_pathway") or "?")
        a = agg[cat]
        a["count"] += 1
        br = 0.0
        tools = r.get("tools")
        if isinstance(tools, list):
            rs = [t.get("rank") for t in tools if isinstance(t, dict) and t.get("rank") is not None]
            if rs:
                br = max(rs)
        if not br and r.get("global_percentile") is not None:
            br = float(r.get("global_percentile") or 0)
        if br > a["best_rank"]:
            a["best_rank"] = br
        d = assessed.get(r.get("event_key"))
        if d == "Accept":
            a["kept"] += 1
        elif d == "Reject":
            a["rejected"] += 1
    out = [{"category": c, **v} for c, v in agg.items()]
    out.sort(key=lambda x: -x["best_rank"])
    return {"pathways": out, "total": len(rows)}


@app.post("/api/pathways/decide")
def api_pathway_decide(payload: dict[str, Any]) -> dict[str, Any]:
    """Keep (Accept) or Reject every interaction in one signalling pathway, at once."""
    _require_index()
    cat = str(payload.get("category") or "")
    kept = bool(payload.get("kept", False))
    expert = payload.get("expert")
    res = index_store.query_events(sort="tools", page=1, page_size=100000)
    rows = [r for r in (res.get("rows", []) if isinstance(res, dict) else [])
            if str(r.get("signaling_pathway") or "?") == cat]
    done, errors = 0, 0
    for r in rows:
        rid = r.get("open_id", r.get("id"))
        if rid is None:
            errors += 1
            continue
        try:
            assessments.set_decision(int(rid), kept, expert)
            done += 1
        except Exception:
            errors += 1
    return {"category": cat, "count": done, "errors": errors, "kept": kept}


@app.get("/api/assessments/keys")
def api_assessed_keys() -> dict[str, str]:
    return assessments.assessed_keys()


@app.post("/api/assessments/clear")
def api_clear_assessments(payload: dict[str, Any]) -> dict[str, Any]:
    """Delete stored expert decisions for the active dataset so the rank cutoff
    governs selection again. {"decision": "Accept"|"Reject"} clears only that kind;
    omit it to clear all."""
    _require_index()
    decision = payload.get("decision")
    if decision not in (None, "Accept", "Reject"):
        raise HTTPException(status_code=422, detail="decision must be 'Accept', 'Reject' or omitted.")
    removed = assessments.clear_decisions(decision)
    return {"removed": removed, "decision": decision, "keys": assessments.assessed_keys()}


@app.get("/api/assessments/export.csv")
def api_export_assessments() -> PlainTextResponse:
    return PlainTextResponse(
        assessments.export_csv(),
        media_type="text/csv",
        headers={"Content-Disposition": 'attachment; filename="expert_assessments.csv"'},
    )


# --------------------------------------------------------------------------- datasets
@app.get("/api/datasets")
def api_datasets() -> dict[str, Any]:
    act = datasets.active_or_none()
    return {"active": act.id if act else None, "items": [d.summary() for d in datasets.all_datasets()]}


@app.put("/api/datasets/clear")
def api_clear_active() -> dict[str, Any]:
    """Unload the current dataset and return to the start state."""
    datasets.clear_active()
    return {"active": None}


def _delete_dataset(ds_id: str) -> dict[str, Any]:
    try:
        datasets.delete(ds_id)
    except PermissionError as exc:
        raise HTTPException(status_code=403, detail=str(exc))
    except LookupError as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    act = datasets.active_or_none()
    return {"deleted": ds_id, "active": act.id if act else None}


@app.post("/api/datasets/{ds_id}/delete")
def api_delete_dataset_post(ds_id: str) -> dict[str, Any]:
    """Remove a dataset the user added (not the shipped example)."""
    return _delete_dataset(ds_id)


@app.delete("/api/datasets/{ds_id}")
def api_delete_dataset(ds_id: str) -> dict[str, Any]:
    return _delete_dataset(ds_id)


@app.post("/api/datasets/upload")
async def api_upload(request: Request, name: str = "", kind: str = "h5ad") -> dict[str, Any]:
    """Receive a file uploaded from the browser and save it under data/uploads/.

    The body is streamed straight to disk, so large .h5ad files do not need to
    fit in memory and no extra dependency (python-multipart) is required.
    """
    safe = Path(name or "").name.strip()
    if not safe:
        raise HTTPException(status_code=400, detail="Missing file name.")
    if kind == "h5ad" and not safe.lower().endswith(".h5ad"):
        raise HTTPException(status_code=400, detail="The single-cell file must be a .h5ad file.")
    dest_dir = config.PROJECT_ROOT / "data" / "uploads"
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / safe
    if dest.exists():
        stem, suf, i = dest.stem, dest.suffix, 2
        while (dest_dir / f"{stem}_{i}{suf}").exists():
            i += 1
        dest = dest_dir / f"{stem}_{i}{suf}"
    size = 0
    try:
        with open(dest, "wb") as fh:
            async for chunk in request.stream():
                fh.write(chunk)
                size += len(chunk)
    except Exception as exc:
        dest.unlink(missing_ok=True)
        raise HTTPException(status_code=500, detail=f"Upload failed: {exc}")
    if size == 0:
        dest.unlink(missing_ok=True)
        raise HTTPException(status_code=400, detail="The uploaded file was empty.")
    rel = str(dest.relative_to(config.PROJECT_ROOT)).replace("\\", "/")
    return {"path": rel, "name": dest.name, "size_mb": round(size / 1e6, 1)}


@app.get("/api/datasets/files")
def api_dataset_files() -> dict[str, Any]:
    return dataset_admin.discover_files()


@app.get("/api/datasets/{ds_id}")
def api_dataset_details(ds_id: str) -> dict[str, Any]:
    data = dataset_admin.details(ds_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Dataset not found.")
    return data


@app.get("/api/workflow")
def api_workflow() -> dict[str, Any]:
    """Where the active dataset stands in the four steps of the platform."""
    ds = datasets.active_or_none()
    if ds is None:
        return {"dataset": None, "prepared": False, "has_results": False, "tools": [], "assessed": 0, "running": []}
    has_results = bool(ds.consensus_file and ds.consensus_file.exists())
    n_assessed = len(assessments.feedback_with_event_keys()) if has_results and index_store.index_is_current() else 0
    return {"dataset": ds.label, "prepared": ds.expression_dir is not None, "has_results": has_results,
            "tools": index_store.methods() if has_results and index_store.index_is_current() else [],
            "assessed": n_assessed, "running": [j for j in jobs.running() if j.get("dataset") == ds.id]}


@app.put("/api/datasets/{ds_id}/view")
def api_set_run(ds_id: str, payload: dict[str, Any]) -> dict[str, Any]:
    """Choose which tool's results of this dataset are shown."""
    try:
        datasets.set_run(ds_id, str(payload.get("run") or ""))
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return {"ok": True}


@app.get("/api/datasets/{ds_id}/advice")
def api_advice(ds_id: str) -> dict[str, Any]:
    """Suitability of each tool for the dataset (independent of what is installed)."""
    data = advisor.advise(ds_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Dataset not found.")
    return data


@app.get("/api/datasets/{ds_id}/provenance")
def api_provenance(ds_id: str) -> dict[str, Any]:
    data = provenance.record(ds_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Dataset not found.")
    return data


@app.get("/api/datasets/{ds_id}/provenance.json")
def api_provenance_file(ds_id: str):
    data = provenance.record(ds_id)
    if not data or not data.get("available"):
        raise HTTPException(status_code=404, detail="No analysis record for this dataset.")
    return JSONResponse(data["record"], headers={
        "Content-Disposition": f'attachment; filename="{ds_id}_analysis_record.json"'})


@app.post("/api/datasets/inspect")
def api_inspect_dataset(payload: dict[str, Any]) -> dict[str, Any]:
    try:
        return dataset_admin.inspect(str(payload.get("h5ad") or ""), payload.get("metadata") or None)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))


@app.post("/api/datasets")
def api_add_dataset(payload: dict[str, Any]) -> dict[str, Any]:
    try:
        return dataset_admin.prepare(payload)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))


@app.post("/api/datasets/{ds_id}/run")
def api_run_methods(ds_id: str, payload: Optional[dict[str, Any]] = None) -> dict[str, Any]:
    try:
        return dataset_admin.run_methods(ds_id, (payload or {}).get("tools"))
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))


@app.get("/api/methods/status")
def api_method_status() -> dict[str, Any]:
    """Fast: configured environment + any tool results already cached. No probing."""
    return settings.method_status()


@app.get("/api/methods/check/{tool}")
def api_check_one(tool: str, force: bool = False) -> dict[str, Any]:
    """Check (load) a single tool on demand. Caches the result for the session."""
    if tool not in settings._TOOL_KEYS:
        raise HTTPException(status_code=404, detail=f"Unknown tool: {tool}")
    res = settings.check_one(tool, force=force)
    return {"tool": tool, "rscript": settings._CHECK_CACHE.get("rscript") or settings.rscript() or "", **res}


@app.get("/api/methods/check")
def api_check_methods(force: bool = False) -> dict[str, Any]:
    return settings.check_methods(force=force)


@app.get("/api/jobs")
def api_jobs() -> list[dict[str, Any]]:
    return dataset_admin.running_jobs()


@app.get("/api/jobs/{job_id}")
def api_job(job_id: str) -> dict[str, Any]:
    job = jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found.")
    return job


@app.get("/api/settings")
def api_settings() -> dict[str, Any]:
    cfg = settings.load()
    key = str(cfg.get("ncbi_api_key", "") or "")
    return {
        "analysis_python": settings.analysis_python() or "",
        "rscript": settings.rscript() or "",
        "tools_python": settings.tools_python() or "",
        "ncbi_email": str(cfg.get("ncbi_email", "") or ""),
        "ncbi_api_key_set": bool(key),
        "ephemeral_mode": settings.ephemeral_mode(),
    }


@app.put("/api/settings/ephemeral")
def api_save_ephemeral(payload: dict[str, Any]) -> dict[str, Any]:
    settings.save({"ephemeral_mode": bool(payload.get("ephemeral_mode"))})
    return {"ephemeral_mode": settings.ephemeral_mode()}


@app.put("/api/settings/ncbi")
def api_save_ncbi(payload: dict[str, Any]) -> dict[str, Any]:
    values: dict[str, Any] = {"ncbi_email": str(payload.get("ncbi_email") or "").strip()}
    if "ncbi_api_key" in payload:  # only overwrite when provided
        values["ncbi_api_key"] = str(payload.get("ncbi_api_key") or "").strip()
    settings.save(values)
    cfg = settings.load()
    return {"ncbi_email": cfg.get("ncbi_email", ""), "ncbi_api_key_set": bool(cfg.get("ncbi_api_key"))}


@app.put("/api/settings/rscript")
def api_save_rscript(payload: dict[str, Any]) -> dict[str, Any]:
    path = settings.clean_path(payload.get("rscript") or "")
    if not settings.rlauncher.exists(path):
        raise HTTPException(status_code=422, detail=f"File not found: {path}")
    settings.save({"rscript": path})
    return settings.check_methods()


@app.put("/api/settings")
def api_save_settings(payload: dict[str, Any]) -> dict[str, Any]:
    path = settings.clean_path(payload.get("analysis_python") or "")
    result = settings.test_analysis_python(path)
    if result["ok"]:
        settings.save({"analysis_python": path})
    return result


@app.put("/api/settings/tools-python")
def api_save_tools_python(payload: dict[str, Any]) -> dict[str, Any]:
    path = settings.clean_path(payload.get("tools_python") or "")
    if not path:                       # clearing it: CellPhoneDB falls back to the analysis Python
        settings.save({"tools_python": ""})
        settings.clear_check_cache()
        return {"ok": True, "cleared": True, "methods": settings.check_methods()}
    result = settings.test_tools_python(path)
    if result["ok"]:
        settings.save({"tools_python": path})
        settings.clear_check_cache()
        result["methods"] = settings.check_methods()
    return result


@app.post("/api/settings/test")
def api_test_settings() -> dict[str, Any]:
    return settings.test_analysis_python()


@app.put("/api/datasets/active")
def api_set_dataset(payload: dict[str, Any]) -> dict[str, Any]:
    try:
        ds = datasets.set_active(str(payload.get("id") or ""))
    except LookupError as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    index_store.status()          # starts building the index if needed
    return ds.summary()


@app.get("/api/results/{result_id}/dataset-evidence")
def api_dataset_evidence(result_id: int) -> dict[str, Any]:
    _require_index()
    data = expression.dataset_evidence(result_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Result not found.")
    return data


# --------------------------------------------------------------------------- method comparison
@app.get("/api/method-comparison")
def api_method_comparison() -> dict[str, Any]:
    _require_index()
    return method_comparison.summary()


@app.get("/api/method-comparison/events")
def api_method_events(combo: str, page: int = 1) -> dict[str, Any]:
    _require_index()
    try:
        return method_comparison.events(combo, page=page)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))


# --------------------------------------------------------------------------- network
@app.get("/api/network")
def api_network(source: str = "consensus", min_methods: int = 2, min_rank: float = 90.0,
                cell_type: list[str] = Query(default=[])) -> dict[str, Any]:
    _require_index()
    return network.network(source, min_methods, min_rank, cell_type or None)


@app.get("/api/network/pair")
def api_network_pair(sender: str, receiver: str, source: str = "consensus",
                     min_methods: int = 2, min_rank: float = 90.0) -> dict[str, Any]:
    _require_index()
    return network.pair(sender, receiver, source, min_methods, min_rank)


# --------------------------------------------------------------------------- curated (aggregated) network
@app.get("/api/curated-network")
def api_curated_network(cutoff: float = 95.0, cell_type: list[str] = Query(default=[])) -> dict[str, Any]:
    _require_index()
    return curated.curated_network(cutoff=cutoff, cell_types=cell_type or None)


@app.get("/api/curated-network/pair")
def api_curated_pair(sender: str, receiver: str, cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return curated.curated_pair(sender, receiver, cutoff=cutoff)


@app.get("/api/hitl-impact")
def api_hitl_impact(cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return curated.hitl_impact(cutoff=cutoff)


@app.get("/api/lr-rewiring")
def api_lr_rewiring(cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return curated.lr_rewiring(cutoff=cutoff)


@app.get("/api/lr-dotplot")
def api_lr_dotplot(cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return curated.lr_dotplot(cutoff=cutoff)


@app.get("/api/functional-enrichment")
def api_functional_enrichment(cutoff: float = 95.0, descriptive: bool = False) -> dict[str, Any]:
    _require_index()
    return functional_enrichment.functional_enrichment(cutoff=cutoff, descriptive=descriptive)


@app.get("/api/discussion/status")
def api_discussion_status() -> dict[str, Any]:
    return discussion.status()


@app.get("/api/discussion")
def api_discussion(cutoff: float = 95.0, refresh: bool = False) -> dict[str, Any]:
    _require_index()
    return discussion.discussion(cutoff=cutoff, force_refresh=refresh)


@app.get("/api/molecular-hubs")
def api_molecular_hubs(cutoff: float = 95.0, top: int = 15) -> dict[str, Any]:
    _require_index()
    return curated.molecular_hubs(cutoff=cutoff, top=top)


@app.get("/api/communication-profiles")
def api_communication_profiles(cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return curated.communication_profiles(cutoff=cutoff)


@app.get("/api/biological-processes")
def api_biological_processes(cutoff: float = 95.0, top: int = 20) -> dict[str, Any]:
    _require_index()
    return curated.biological_processes(cutoff=cutoff, top=top)


@app.get("/api/biological-processes/detail")
def api_biological_process_detail(source: str, name: str, cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return curated.biological_process_detail(source, name, cutoff=cutoff)


@app.get("/api/molecular-hubs/detail")
def api_molecular_hub_detail(kind: str, name: str, cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return curated.molecular_hub_detail(kind, name, cutoff=cutoff)


@app.get("/api/communication-profiles/detail")
def api_communication_profile_detail(population: str, cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return curated.communication_profile_detail(population, cutoff=cutoff)


@app.get("/api/hitl-impact/detail")
def api_hitl_impact_detail(bucket: str, cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return curated.hitl_impact_detail(bucket, cutoff=cutoff)


@app.get("/api/evidence-aggregation")
def api_evidence_aggregation(cutoff: float = 95.0, batch: int = 40) -> dict[str, Any]:
    _require_index()
    return evidence.evidence_aggregation(cutoff=cutoff, batch=batch)


@app.get("/api/evidence-aggregation/detail")
def api_evidence_source_detail(source: str, cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return evidence.evidence_source_detail(source, cutoff=cutoff)


@app.get("/api/biological-insights")
def api_biological_insights(cutoff: float = 95.0) -> dict[str, Any]:
    _require_index()
    return insights.biological_insights(cutoff=cutoff)


# --------------------------------------------------------------------------- receptor explorer
@app.get("/api/receptors")
def api_receptors(top: str = "all", min_events: int = 5) -> dict[str, Any]:
    _require_index()
    try:
        return receptor_explorer.recurrence(top, min_events)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))


@app.get("/api/receptors/detail")
def api_receptor_detail(key: str, top: str = "all") -> dict[str, Any]:
    _require_index()
    data = receptor_explorer.receptor(key, top)
    if data is None:
        raise HTTPException(status_code=404, detail="Receptor not found.")
    return data


@app.get("/api/receptors/ligand-enrichment")
def api_receptor_enrichment(key: str, top: str = "all") -> dict[str, Any]:
    _require_index()
    try:
        return receptor_explorer.ligand_enrichment(key, top)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except Exception as exc:          # network problems reaching g:Profiler, etc.
        raise HTTPException(status_code=502, detail=f"g:Profiler could not be reached: {exc}")


