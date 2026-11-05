"""
Minimal Gemini client shared by the platform's AI-assisted components.

As in GeneCurator, the API key is never written in the code: it is read from
the environment variable GEMINI_API_KEY, normally set in a .env file in the
project root (ignored by git). Without a key every AI-assisted component is
simply off and the platform uses its rule-based behaviour.

Plain HTTPS (urllib) — no extra packages are needed.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError
from urllib.request import Request, urlopen

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = "gemini-2.5-flash"
FALLBACK_MODELS = ["gemini-2.5-flash", "gemini-2.0-flash", "gemini-flash-latest"]
API = "https://generativelanguage.googleapis.com/v1beta/models"


def load_dotenv() -> None:
    """Load PROJECT_ROOT/.env into the environment (python-dotenv if installed,
    otherwise a minimal KEY=VALUE reader). Existing variables are kept."""
    env_file = PROJECT_ROOT / ".env"
    try:
        from dotenv import load_dotenv as _ld
        _ld(env_file)
        return
    except ImportError:
        pass
    try:
        for line in env_file.read_text(encoding="utf-8-sig").splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, _, v = line.partition("=")
            k, v = k.strip(), v.strip().strip('"').strip("'")
            if k and k not in os.environ:
                os.environ[k] = v
    except OSError:
        pass


load_dotenv()


def api_key() -> str:
    return os.environ.get("GEMINI_API_KEY", "").strip()


def available() -> bool:
    return bool(api_key())


def model_name() -> str:
    return os.environ.get("GEMINI_MODEL", "").strip() or DEFAULT_MODEL


def _post(model: str, payload: dict[str, Any], timeout: int = 120) -> dict[str, Any]:
    req = Request(f"{API}/{model}:generateContent",
                  data=json.dumps(payload).encode("utf-8"),
                  headers={"Content-Type": "application/json", "x-goog-api-key": api_key()},
                  method="POST")
    with urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def generate_json(system: str, user: str, *, temperature: float = 0.2,
                  max_tokens: int = 8192) -> tuple[dict[str, Any], str]:
    """One prompt in, one JSON object out. Returns (parsed_json, model_used).
    Retries briefly on overload/rate limit and skips models the key cannot use.
    Raises RuntimeError if no valid JSON could be obtained."""
    if not available():
        raise RuntimeError("GEMINI_API_KEY not set.")
    payload = {
        "systemInstruction": {"parts": [{"text": system}]},
        "contents": [{"role": "user", "parts": [{"text": user}]}],
        "generationConfig": {"temperature": temperature, "maxOutputTokens": max_tokens,
                             "responseMimeType": "application/json"},
    }
    want = model_name()
    candidates = [want] + [m for m in FALLBACK_MODELS if m != want]
    last: Exception | None = None
    for model in candidates:
        for attempt in range(3):
            try:
                out = _post(model, payload)
            except HTTPError as e:
                last = e
                if e.code in (429, 503):
                    time.sleep(1.5 * (attempt + 1))
                    continue
                if e.code == 404:
                    break           # model not available on this key: try the next one
                raise RuntimeError(f"Gemini request failed (HTTP {e.code}).") from e
            cands = out.get("candidates") or []
            text = "".join(p.get("text", "") for p in
                           ((cands[0].get("content") or {}).get("parts") or [])) if cands else ""
            text = text.strip()
            if text.startswith("```"):
                text = text.strip("`")
                text = text[text.find("{"):] if "{" in text else text
            try:
                return json.loads(text), model
            except json.JSONDecodeError as e:
                last = e
                break
    raise RuntimeError(f"Gemini returned no usable output ({last}).")
