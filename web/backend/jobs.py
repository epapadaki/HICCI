"""
Background jobs: long steps (preparing a dataset, later running CCI methods)
run as separate processes of the analysis environment. Their output is
written to web/workspace/jobs/<id>.log and can be followed from the browser.
"""

from __future__ import annotations

import json
import os
import subprocess
import threading
import time
import uuid
from typing import Any, Callable

from . import config

JOBS_DIR = config.WORKSPACE_DIR / "jobs"
JOBS_DIR.mkdir(parents=True, exist_ok=True)
RESULT_MARK = "@@RESULT@@ "
STEP_MARK = "@@STEP@@ "

_jobs: dict[str, dict[str, Any]] = {}
_lock = threading.Lock()


def start(kind: str, title: str, cmd: list[str], *, cwd: str | None = None, dataset: str | None = None,
          on_success: Callable[[dict[str, Any] | None], None] | None = None) -> dict[str, Any]:
    job_id = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
    job = {"id": job_id, "kind": kind, "title": title, "state": "running",
           "started": time.strftime("%Y-%m-%d %H:%M:%S"), "finished": None,
           "result": None, "error": None, "step": None, "dataset": None}
    job["dataset"] = dataset
    with _lock:
        _jobs[job_id] = job
    log_path = JOBS_DIR / f"{job_id}.log"

    def run():
        result = None
        try:
            with open(log_path, "w", encoding="utf-8") as log:
                env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUTF8="1", PYTHONUNBUFFERED="1")
                proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        text=True, encoding="utf-8", errors="replace", bufsize=1)
                for line in proc.stdout:
                    if line.startswith(STEP_MARK):
                        job["step"] = line[len(STEP_MARK):].strip()
                        continue
                    if line.startswith(RESULT_MARK):
                        try:
                            result = json.loads(line[len(RESULT_MARK):])
                        except json.JSONDecodeError:
                            pass
                        continue
                    log.write(line)
                    log.flush()
                code = proc.wait()
            if code == 0:
                if on_success:
                    on_success(result)
                job.update(state="done", result=result)
            else:
                job.update(state="failed", error=_last_lines(log_path, 3) or f"Exit code {code}")
        except Exception as exc:          # could not start the process, etc.
            job.update(state="failed", error=str(exc))
        job["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")

    threading.Thread(target=run, daemon=True).start()
    return dict(job)


def _last_lines(path, n: int) -> str:
    try:
        lines = [l for l in path.read_text(encoding="utf-8", errors="replace").splitlines() if l.strip()]
    except OSError:
        return ""
    return "\n".join(lines[-n:])


def get(job_id: str) -> dict[str, Any] | None:
    with _lock:
        job = _jobs.get(job_id)
    if not job:
        return None
    out = dict(job)
    try:
        out["log"] = (JOBS_DIR / f"{job_id}.log").read_text(encoding="utf-8", errors="replace")[-6000:]
    except OSError:
        out["log"] = ""
    return out


def running(kind: str | None = None) -> list[dict[str, Any]]:
    with _lock:
        return [dict(j) for j in _jobs.values() if j["state"] == "running" and (kind is None or j["kind"] == kind)]
