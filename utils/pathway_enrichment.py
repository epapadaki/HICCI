from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[1]

CACHE_DB = (
    PROJECT_ROOT
    / "data"
    / "enrichment"
    / "gprofiler_cache.sqlite"
)

GPROFILER_URL = (
    "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
)

DEFAULT_SOURCES = [
    "GO:BP",
    "REAC",
    "KEGG",
    "WP",
]


def _utc_now() -> str:
    return datetime.now(
        timezone.utc
    ).isoformat(
        timespec="seconds"
    )


def _clean_genes(
    genes,
) -> list[str]:
    cleaned = {
        str(
            gene
        )
        .strip()
        .upper()

        for gene
        in genes

        if str(
            gene
        ).strip()
    }

    return sorted(
        cleaned
    )


def _cache_key(
    *,
    organism: str,
    query_genes: list[str],
    background_genes: list[str],
    sources: list[str],
    significance_method: str,
    user_threshold: float,
) -> str:
    payload = {
        "organism":
            organism,

        "query_genes":
            query_genes,

        "background_genes":
            background_genes,

        "sources":
            sorted(
                sources
            ),

        "significance_method":
            significance_method,

        "user_threshold":
            user_threshold,
    }

    encoded = json.dumps(
        payload,
        sort_keys=True,
        ensure_ascii=False,
    ).encode(
        "utf-8"
    )

    return hashlib.sha256(
        encoded
    ).hexdigest()


def _connect() -> sqlite3.Connection:
    CACHE_DB.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    con = sqlite3.connect(
        CACHE_DB,
        timeout=30,
    )

    con.row_factory = sqlite3.Row

    return con


def init_cache() -> None:
    con = _connect()

    con.execute(
        """
        CREATE TABLE IF NOT EXISTS gprofiler_results (
            cache_key TEXT PRIMARY KEY,
            organism TEXT NOT NULL,
            query_genes_json TEXT NOT NULL,
            background_genes_json TEXT NOT NULL,
            sources_json TEXT NOT NULL,
            result_json TEXT NOT NULL,
            meta_json TEXT,
            created_at TEXT NOT NULL
        )
        """
    )

    con.commit()
    con.close()


def _load_cached(
    cache_key: str,
) -> dict[str, Any] | None:
    init_cache()

    con = _connect()

    row = con.execute(
        """
        SELECT
            result_json,
            meta_json,
            created_at
        FROM gprofiler_results
        WHERE cache_key = ?
        LIMIT 1
        """,
        (
            cache_key,
        ),
    ).fetchone()

    con.close()

    if row is None:
        return None

    return {
        "results":
            json.loads(
                row[
                    "result_json"
                ]
            ),

        "meta":
            (
                json.loads(
                    row[
                        "meta_json"
                    ]
                )
                if row[
                    "meta_json"
                ]
                else {}
            ),

        "cached":
            True,

        "cached_at":
            row[
                "created_at"
            ],
    }


def _save_cached(
    *,
    cache_key: str,
    organism: str,
    query_genes: list[str],
    background_genes: list[str],
    sources: list[str],
    results: list[dict[str, Any]],
    meta: dict[str, Any],
) -> None:
    init_cache()

    con = _connect()

    con.execute(
        """
        INSERT INTO gprofiler_results (
            cache_key,
            organism,
            query_genes_json,
            background_genes_json,
            sources_json,
            result_json,
            meta_json,
            created_at
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(cache_key)
        DO UPDATE SET
            result_json = excluded.result_json,
            meta_json = excluded.meta_json,
            created_at = excluded.created_at
        """,
        (
            cache_key,
            organism,
            json.dumps(
                query_genes,
                ensure_ascii=False,
            ),
            json.dumps(
                background_genes,
                ensure_ascii=False,
            ),
            json.dumps(
                sources,
                ensure_ascii=False,
            ),
            json.dumps(
                results,
                ensure_ascii=False,
            ),
            json.dumps(
                meta,
                ensure_ascii=False,
            ),
            _utc_now(),
        ),
    )

    con.commit()
    con.close()


def _post_gprofiler(
    payload: dict[str, Any],
) -> dict[str, Any]:
    body = json.dumps(
        payload,
        ensure_ascii=False,
    ).encode(
        "utf-8"
    )

    request = Request(
        GPROFILER_URL,
        data=body,
        method="POST",
        headers={
            "Content-Type":
                "application/json",

            "Accept":
                "application/json",

            "User-Agent":
                "CCC-Platform/1.0",
        },
    )

    try:
        with urlopen(
            request,
            timeout=60,
        ) as response:

            raw = response.read().decode(
                "utf-8"
            )

    except HTTPError as error:
        detail = error.read().decode(
            "utf-8",
            errors="replace",
        )

        raise RuntimeError(
            "g:Profiler enrichment request failed "
            f"with HTTP {error.code}: {detail}"
        ) from error

    except URLError as error:
        raise RuntimeError(
            "Could not connect to g:Profiler: "
            f"{error}"
        ) from error

    try:
        return json.loads(
            raw
        )

    except json.JSONDecodeError as error:
        raise RuntimeError(
            "g:Profiler returned invalid JSON."
        ) from error


def run_gprofiler_enrichment(
    *,
    query_genes,
    background_genes,
    organism: str = "hsapiens",
    sources: list[str] | None = None,
    significance_method: str = "fdr",
    user_threshold: float = 0.05,
    force_refresh: bool = False,
) -> dict[str, Any]:
    """
    Run g:Profiler g:GOSt over-representation analysis
    using the supplied custom gene background.

    The default use in the CCC platform is:
      query = unique ligand partners of one receptor
              in the selected top-ranked CCC event set
      background = all unique ligands in the active
                   ranked canonical CCC event universe
    """

    organism = str(
        organism
    ).strip()

    if not organism:
        raise ValueError(
            "A g:Profiler organism identifier is required."
        )

    query = _clean_genes(
        query_genes
    )

    background = _clean_genes(
        background_genes
    )

    if not query:
        raise ValueError(
            "The enrichment query contains no genes."
        )

    if len(
        query
    ) < 3:
        raise ValueError(
            "At least three unique query genes are required "
            "for this exploratory enrichment analysis."
        )

    if not background:
        raise ValueError(
            "The enrichment background contains no genes."
        )

    background_set = set(
        background
    )

    query = [
        gene
        for gene
        in query
        if gene in background_set
    ]

    if len(
        query
    ) < 3:
        raise ValueError(
            "Fewer than three query genes remain after "
            "intersecting the query with the CCC background."
        )

    resolved_sources = (
        list(
            sources
        )
        if sources
        else DEFAULT_SOURCES.copy()
    )

    cache_key = _cache_key(
        organism=organism,
        query_genes=query,
        background_genes=background,
        sources=resolved_sources,
        significance_method=significance_method,
        user_threshold=user_threshold,
    )

    if not force_refresh:
        cached = _load_cached(
            cache_key
        )

        if cached is not None:
            return {
                "organism":
                    organism,

                "query_genes":
                    query,

                "background_genes":
                    background,

                "sources":
                    resolved_sources,

                "results":
                    cached[
                        "results"
                    ],

                "meta":
                    cached[
                        "meta"
                    ],

                "cached":
                    True,

                "cached_at":
                    cached[
                        "cached_at"
                    ],
            }

    payload = {
        "organism":
            organism,

        "query":
            query,

        "sources":
            resolved_sources,

        "user_threshold":
            float(
                user_threshold
            ),

        "significance_threshold_method":
            significance_method,

        "ordered":
            False,

        "no_evidences":
            True,

        "domain_scope":
            "custom",

        "background":
            background,

        "all_results":
            False,
    }

    response = _post_gprofiler(
        payload
    )

    results = response.get(
        "result",
        [],
    )

    meta = response.get(
        "meta",
        {},
    )

    if not isinstance(
        results,
        list,
    ):
        raise RuntimeError(
            "Unexpected g:Profiler response format."
        )

    _save_cached(
        cache_key=cache_key,
        organism=organism,
        query_genes=query,
        background_genes=background,
        sources=resolved_sources,
        results=results,
        meta=meta,
    )

    return {
        "organism":
            organism,

        "query_genes":
            query,

        "background_genes":
            background,

        "sources":
            resolved_sources,

        "results":
            results,

        "meta":
            meta,

        "cached":
            False,

        "cached_at":
            None,
    }
