"""
How to start a Python interpreter for the Python-based CCI tools (CellPhoneDB,
NATMI): a normal executable (the Windows/Linux/macOS Python), or a Python that
lives inside WSL (Ubuntu on Windows), written as "wsl:<linux path>".

This mirrors ``rlauncher`` but for Python. The important difference from R is
that the Python tool steps are invoked with *absolute* path arguments (the
script to run, ``--emFile``, ``--out`` …). When the interpreter is inside WSL,
every absolute Windows path argument is rewritten to its ``/mnt/<drive>/…``
form so the Linux Python can find it; plain flags and values are passed
through unchanged.

A WSL Python given as a conda env (``.../envs/<name>/bin/python``) has that
environment activated first, exactly as in a terminal. A plain venv (any other
``.../bin/python``) is run directly.

Why this exists: on a locked-down Windows (Application Control / WDAC) the
freshly pip-installed native DLLs that SciPy/scikit-learn ship are blocked, so
CellPhoneDB cannot import there. Running it inside WSL side-steps the Windows
policy, the same way the R tools already run inside WSL.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

# Reuse the shared helpers so the "wsl:" convention stays identical to R.
from rlauncher import WSL_PREFIX, is_wsl, to_wsl_path  # noqa: F401

_WIN_ABS = re.compile(r"^[A-Za-z]:[\\/]")


def _conv(arg: str) -> str:
    """Rewrite an absolute Windows path argument to its /mnt form; leave
    flags, option values and relative tokens untouched."""
    return to_wsl_path(arg) if _WIN_ABS.match(str(arg)) else str(arg)


def command(python: str, args: list[str], cwd: str | Path | None = None) -> list[str]:
    """Command list that runs ``python args``.

    For a normal interpreter this is just ``[python, *args]`` — identical to the
    previous behaviour, so nothing changes when no WSL Python is configured.
    For a ``wsl:`` Python the command runs inside WSL, with absolute Windows
    path arguments converted to /mnt paths and a conda env activated if needed.
    """
    if not is_wsl(python):
        return [python, *args]
    linux = python[len(WSL_PREFIX):].strip()
    conv = [_conv(a) for a in args]
    m = re.match(r"^(.*)/envs/([^/]+)/bin/(python[0-9.]*)$", linux)
    parts: list[str] = []
    if cwd is not None:
        parts.append("cd " + shlex.quote(to_wsl_path(cwd)))
    if m:
        base, env = m.group(1), m.group(2)
        parts.append("source " + shlex.quote(base + "/etc/profile.d/conda.sh"))
        parts.append("conda activate " + shlex.quote(env))
        parts.append("python " + " ".join(shlex.quote(a) for a in conv))
    else:
        parts.append(shlex.quote(linux) + " " + " ".join(shlex.quote(a) for a in conv))
    return ["wsl", "bash", "-lc", " && ".join(parts)]


def exists(python: str) -> bool:
    return True if is_wsl(python) else Path(python).exists()
