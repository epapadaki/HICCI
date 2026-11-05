"""
How to start Rscript: a normal executable (R on Windows/Linux/macOS), or an
R that lives inside WSL (Ubuntu on Windows), written as "wsl:<linux path>".

For a conda R inside WSL (.../envs/<name>/bin/Rscript) the environment is
activated first, as it would be in a terminal.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

WSL_PREFIX = "wsl:"


def is_wsl(rscript: str) -> bool:
    return str(rscript or "").lower().startswith(WSL_PREFIX)


def to_wsl_path(path: str | Path) -> str:
    """C:\\Users\\x\\y  ->  /mnt/c/Users/x/y"""
    p = str(path)
    if not re.match(r"^[A-Za-z]:[\\\\/]", p):
        p = str(Path(p).resolve())
    p = p.replace("\\", "/")
    m = re.match(r"^([A-Za-z]):/(.*)$", p)
    return f"/mnt/{m.group(1).lower()}/{m.group(2)}" if m else p


def command(rscript: str, args: list[str], cwd: str | Path | None = None) -> list[str]:
    """Command list that runs Rscript with args (args are passed through unchanged)."""
    if not is_wsl(rscript):
        return [rscript, *args]
    linux = rscript[len(WSL_PREFIX):].strip()
    m = re.match(r"^(.*)/envs/([^/]+)/bin/Rscript$", linux)
    parts = []
    if cwd is not None:
        parts.append("cd " + shlex.quote(to_wsl_path(cwd)))
    if m:
        base, env = m.group(1), m.group(2)
        parts.append(f"source {shlex.quote(base + '/etc/profile.d/conda.sh')}")
        parts.append(f"conda activate {shlex.quote(env)}")
        parts.append("Rscript " + " ".join(shlex.quote(a) for a in args))
    else:
        parts.append(shlex.quote(linux) + " " + " ".join(shlex.quote(a) for a in args))
    return ["wsl", "bash", "-lc", " && ".join(parts)]


def exists(rscript: str) -> bool:
    return True if is_wsl(rscript) else Path(rscript).exists()
