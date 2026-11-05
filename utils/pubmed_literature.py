from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import sqlite3
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parents[1]

CACHE_DB = (
    PROJECT_ROOT
    / "data"
    / "literature"
    / "pubmed_cache.sqlite"
)

EUTILS_BASE = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
)

TOOL_NAME = "ccc_platform"

DEFAULT_MAX_RESULTS = 8


# --------------------------------------------------
# Small helpers
# --------------------------------------------------
def _utc_now() -> str:
    return datetime.now(
        timezone.utc
    ).isoformat(
        timespec="seconds"
    )


def _clean_symbol(
    value: str,
) -> str:
    return str(
        value
    ).strip().upper()


def _clean_text(
    value: Any,
) -> str | None:
    if value is None:
        return None

    text = re.sub(
        r"\s+",
        " ",
        str(
            value
        ),
    ).strip()

    return text or None


def _xml_text(
    element: ET.Element | None,
) -> str | None:
    if element is None:
        return None

    text = "".join(
        element.itertext()
    )

    return _clean_text(
        text
    )


def _json_dumps(
    value: Any,
) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
    )


def _json_loads(
    value: str | None,
    default: Any,
) -> Any:
    if not value:
        return default

    try:
        return json.loads(
            value
        )

    except Exception:
        return default


# --------------------------------------------------
# SQLite cache
# --------------------------------------------------
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

    con.executescript(
        """
        CREATE TABLE IF NOT EXISTS pubmed_articles (
            pmid TEXT PRIMARY KEY,
            title TEXT,
            abstract TEXT,
            journal TEXT,
            year TEXT,
            authors_json TEXT,
            mesh_json TEXT,
            publication_types_json TEXT,
            doi TEXT,
            fetched_at TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS pair_search_cache (
            ligand TEXT NOT NULL,
            receptor TEXT NOT NULL,
            query TEXT NOT NULL,
            max_results INTEGER NOT NULL,
            pmids_json TEXT NOT NULL,
            fetched_at TEXT NOT NULL,
            PRIMARY KEY (
                ligand,
                receptor,
                max_results
            )
        );
        """
    )

    con.commit()
    con.close()


def clear_pair_cache(
    ligand: str,
    receptor: str,
) -> None:
    init_cache()

    ligand = _clean_symbol(
        ligand
    )

    receptor = _clean_symbol(
        receptor
    )

    con = _connect()

    con.execute(
        """
        DELETE FROM pair_search_cache
        WHERE ligand = ?
          AND receptor = ?
        """,
        (
            ligand,
            receptor,
        ),
    )

    con.commit()
    con.close()


# --------------------------------------------------
# Optional UniProt-derived aliases
# --------------------------------------------------
def _protein_name_alias(
    gene_symbol: str,
) -> str | None:
    try:
        from utils.biological_context import (
            get_gene_annotation,
        )

        annotation = (
            get_gene_annotation(
                gene_symbol
            )
        )

        if not annotation.get(
            "found"
        ):
            return None

        protein_name = (
            annotation.get(
                "protein_name"
            )
        )

        if not protein_name:
            return None

        # Keep the alias concise enough for PubMed search.
        protein_name = str(
            protein_name
        ).strip()

        if (
            len(
                protein_name
            )
            > 120
        ):
            return None

        return protein_name

    except Exception:
        return None


def _escape_pubmed_phrase(
    value: str,
) -> str:
    return str(
        value
    ).replace(
        '"',
        "",
    ).strip()


def build_pair_query(
    ligand: str,
    receptor: str,
) -> str:
    ligand = _clean_symbol(
        ligand
    )

    receptor = _clean_symbol(
        receptor
    )

    ligand_alias = (
        _protein_name_alias(
            ligand
        )
    )

    receptor_alias = (
        _protein_name_alias(
            receptor
        )
    )

    ligand_terms = [
        f'"{_escape_pubmed_phrase(ligand)}"[Title/Abstract]'
    ]

    receptor_terms = [
        f'"{_escape_pubmed_phrase(receptor)}"[Title/Abstract]'
    ]

    if (
        ligand_alias
        and ligand_alias.upper()
        != ligand
    ):
        ligand_terms.append(
            f'"{_escape_pubmed_phrase(ligand_alias)}"[Title/Abstract]'
        )

    if (
        receptor_alias
        and receptor_alias.upper()
        != receptor
    ):
        receptor_terms.append(
            f'"{_escape_pubmed_phrase(receptor_alias)}"[Title/Abstract]'
        )

    ligand_block = (
        "("
        + " OR ".join(
            ligand_terms
        )
        + ")"
    )

    receptor_block = (
        "("
        + " OR ".join(
            receptor_terms
        )
        + ")"
    )

    return (
        ligand_block
        + " AND "
        + receptor_block
    )



def _normalize_gene_symbols(
    values: list[str] | tuple[str, ...] | None,
) -> list[str]:
    cleaned = []
    seen = set()

    for value in values or []:
        symbol = _clean_symbol(
            value
        )

        if not symbol:
            continue

        if symbol in seen:
            continue

        seen.add(
            symbol
        )
        cleaned.append(
            symbol
        )

    return cleaned


def _normalize_complex_lookup_name(
    value: Any,
) -> str | None:
    """
    Normalize a method-specific complex label for local lookup.

    Important: this helper is only used to match against curated
    local complex snapshots. It is not used to infer subunits by
    splitting arbitrary display names.
    """
    value = _clean_text(
        value
    )

    if value is None:
        return None

    lowered = value.lower()

    for prefix in (
        "complex:",
        "heteromer:",
        "multimer:",
    ):
        if lowered.startswith(
            prefix
        ):
            value = value[
                len(prefix):
            ].strip()
            break

    value = value.replace(
        " ",
        "_",
    )

    return value.upper() or None


_COMPLEX_SNAPSHOT_CACHE: dict[str, list[str]] | None = None


def _load_local_complex_subunit_lookup() -> dict[str, list[str]]:
    """
    Load local complex-name -> gene-subunit mappings.

    This uses platform-local snapshots only. It allows PubMed search
    to use documented complex subunits instead of artificial display
    labels such as COMPLEX:PLEXINA4_COMPLEX1.
    """
    global _COMPLEX_SNAPSHOT_CACHE

    if _COMPLEX_SNAPSHOT_CACHE is not None:
        return _COMPLEX_SNAPSHOT_CACHE

    lookup: dict[str, list[str]] = {}

    snapshot_paths = [
        PROJECT_ROOT
        / "data"
        / "biological_context"
        / "cellphonedb_complexes"
        / "cellphonedb_database_1.json",
    ]

    for snapshot_path in snapshot_paths:
        if not snapshot_path.exists():
            continue

        try:
            payload = json.loads(
                snapshot_path.read_text(
                    encoding="utf-8"
                )
            )
        except Exception:
            continue

        for record in payload.get(
            "complexes",
            [],
        ):
            if not isinstance(
                record,
                dict,
            ):
                continue

            subunits = _normalize_gene_symbols(
                record.get(
                    "subunits"
                )
                or []
            )

            if not subunits:
                continue

            candidate_names = [
                record.get(
                    "complex_name"
                ),
                record.get(
                    "canonical_entity_key"
                ),
            ]

            for candidate in candidate_names:
                key = _normalize_complex_lookup_name(
                    candidate
                )

                if key:
                    lookup[
                        key
                    ] = subunits

    _COMPLEX_SNAPSHOT_CACHE = lookup
    return lookup


def _subunits_from_local_complex_snapshot(
    *values: Any,
) -> list[str]:
    lookup = _load_local_complex_subunit_lookup()

    for value in values:
        key = _normalize_complex_lookup_name(
            value
        )

        if key and key in lookup:
            return list(
                lookup[
                    key
                ]
            )

    return []




def _looks_like_complex_label(
    value: Any,
) -> bool:
    text = _clean_text(
        value
    )

    if not text:
        return False

    upper = text.upper()

    return (
        upper.startswith("COMPLEX:")
        or "_COMPLEX" in upper
        or upper.endswith("_COMPLEX")
        or "RECEPTOR_TYPE" in upper
    )


def _resolve_required_subunits_or_complex_labels(
    required_subunits: list[str] | tuple[str, ...] | None,
    *fallback_values: Any,
) -> list[str]:
    """
    Resolve a required-subunit argument into real gene symbols.

    Older complex-aware UI code may pass a complex display label such as
    COMPLEX:IL20_RECEPTOR_TYPE_I in the required_subunits slot. That value is
    not a gene subunit and must not be searched directly in PubMed. When this
    happens, resolve the label through the local complex snapshot first.
    """
    raw_values = [
        value
        for value
        in (required_subunits or [])
        if _clean_text(
            value
        )
    ]

    if raw_values and all(
        _looks_like_complex_label(
            value
        )
        for value
        in raw_values
    ):
        resolved = _subunits_from_local_complex_snapshot(
            *raw_values,
            *fallback_values,
        )

        if resolved:
            return resolved

        for value in raw_values:
            resolved = _subunits_from_entity_key(
                str(
                    value
                )
            )

            if resolved:
                return resolved

        return []

    return _normalize_gene_symbols(
        raw_values
    )


def _subunits_from_entity_key(
    entity_key: str | None,
) -> list[str]:
    if not entity_key:
        return []

    value = str(
        entity_key
    ).strip()

    if not value.upper().startswith(
        "COMPLEX:"
    ):
        return []

    # First try platform-local curated complex snapshots. This is
    # required for CellPhoneDB display labels such as
    # COMPLEX:PLEXINA4_COMPLEX1, where the real subunits are stored
    # in the local complex database rather than encoded in the label.
    snapshot_subunits = _subunits_from_local_complex_snapshot(
        value
    )

    if snapshot_subunits:
        return snapshot_subunits

    payload = value.split(
        ":",
        1,
    )[1]

    # Only plus-separated canonical complex keys are interpreted as
    # explicit subunits. We deliberately do NOT split underscores here:
    # labels such as PLEXINA4_COMPLEX1 are method-specific complex
    # identifiers, not gene lists.
    if "+" not in payload:
        return []

    return _normalize_gene_symbols(
        payload.split(
            "+"
        )
    )


def _pubmed_terms_for_gene(
    gene_symbol: str,
) -> list[str]:
    gene_symbol = _clean_symbol(
        gene_symbol
    )

    if not gene_symbol:
        return []

    terms = [
        gene_symbol,
    ]

    alias = _protein_name_alias(
        gene_symbol
    )

    if (
        alias
        and alias.upper()
        != gene_symbol
    ):
        terms.append(
            alias
        )

    return _unique_preserve_order(
        terms
    )


def _pubmed_or_block(
    terms: list[str],
) -> str:
    cleaned_terms = [
        _escape_pubmed_phrase(
            term
        )
        for term
        in terms
        if _escape_pubmed_phrase(
            term
        )
    ]

    if not cleaned_terms:
        return ""

    return (
        "("
        + " OR ".join(
            f'"{term}"[Title/Abstract]'
            for term
            in cleaned_terms
        )
        + ")"
    )


def _entity_query_block(
    *,
    display_name: str,
    entity_key: str | None = None,
    required_subunits: list[str] | tuple[str, ...] | None = None,
) -> tuple[str, dict[str, Any]]:
    """
    Build a PubMed query block for a molecular entity.

    For verified complexes, do not search the artificial display
    label as if it were a gene. Instead require the documented
    subunits to be present in the title/abstract query.
    """
    display_symbol = _clean_symbol(
        display_name
    )

    subunits = _resolve_required_subunits_or_complex_labels(
        required_subunits,
        display_name,
        entity_key,
    )

    if not subunits:
        subunits = _subunits_from_entity_key(
            entity_key
        )

    if not subunits:
        subunits = _subunits_from_local_complex_snapshot(
            display_name,
            entity_key,
        )

    is_complex = bool(
        subunits
    ) or (
        str(
            entity_key
            or ""
        )
        .strip()
        .upper()
        .startswith(
            "COMPLEX:"
        )
    )

    if is_complex and subunits:
        subunit_blocks = []

        for subunit in subunits:
            gene_terms = _pubmed_terms_for_gene(
                subunit
            )
            block = _pubmed_or_block(
                gene_terms
            )

            if block:
                subunit_blocks.append(
                    block
                )

        if subunit_blocks:
            return (
                "("
                + " AND ".join(
                    subunit_blocks
                )
                + ")",
                {
                    "entity_type":
                        "complex",
                    "display_name":
                        display_name,
                    "entity_key":
                        entity_key,
                    "required_subunits":
                        subunits,
                    "query_rule":
                        "all_required_subunits",
                },
            )

    terms = _pubmed_terms_for_gene(
        display_symbol
    )

    return (
        _pubmed_or_block(
            terms
        ),
        {
            "entity_type":
                "gene",
            "display_name":
                display_name,
            "entity_key":
                entity_key,
            "required_subunits":
                [],
            "query_rule":
                "single_gene_or_alias",
        },
    )


def build_entity_aware_pair_query(
    ligand: str,
    receptor: str,
    *,
    ligand_entity_key: str | None = None,
    receptor_entity_key: str | None = None,
    ligand_required_subunits: list[str] | tuple[str, ...] | None = None,
    receptor_required_subunits: list[str] | tuple[str, ...] | None = None,
) -> tuple[str, dict[str, Any]]:
    ligand_block, ligand_metadata = _entity_query_block(
        display_name=ligand,
        entity_key=ligand_entity_key,
        required_subunits=ligand_required_subunits,
    )

    receptor_block, receptor_metadata = _entity_query_block(
        display_name=receptor,
        entity_key=receptor_entity_key,
        required_subunits=receptor_required_subunits,
    )

    if not ligand_block or not receptor_block:
        raise ValueError(
            "Ligand and receptor query terms are required."
        )

    return (
        ligand_block
        + " AND "
        + receptor_block,
        {
            "ligand":
                ligand_metadata,
            "receptor":
                receptor_metadata,
        },
    )


def _entity_cache_key(
    display_name: str,
    entity_key: str | None,
    required_subunits: list[str] | tuple[str, ...] | None,
) -> str:
    subunits = _resolve_required_subunits_or_complex_labels(
        required_subunits,
        display_name,
        entity_key,
    )

    if not subunits:
        subunits = _subunits_from_entity_key(
            entity_key
        )

    if not subunits:
        subunits = _subunits_from_local_complex_snapshot(
            display_name,
            entity_key,
        )

    if subunits:
        return (
            "COMPLEX:"
            + "+".join(
                subunits
            )
        )

    if entity_key:
        key = str(
            entity_key
        ).strip()

        if key:
            return key.upper()

    return _clean_symbol(
        display_name
    )


def _extract_entity_pair_snippet(
    article: dict[str, Any],
    *,
    ligand_terms: list[str],
    receptor_term_groups: list[list[str]],
) -> str | None:
    abstract = _clean_text(
        article.get(
            "abstract"
        )
    )

    if not abstract:
        return None

    ligand_patterns = [
        _term_pattern(
            term
        )
        for term
        in ligand_terms
        if term
    ]

    receptor_patterns_by_group = [
        [
            _term_pattern(
                term
            )
            for term
            in group
            if term
        ]
        for group
        in receptor_term_groups
        if group
    ]

    if not ligand_patterns or not receptor_patterns_by_group:
        return None

    sentences = re.split(
        r"(?<=[.!?])\s+",
        abstract,
    )

    for sentence in sentences:
        ligand_hit = any(
            pattern.search(
                sentence
            )
            for pattern
            in ligand_patterns
        )

        receptor_hit = all(
            any(
                pattern.search(
                    sentence
                )
                for pattern
                in group
            )
            for group
            in receptor_patterns_by_group
        )

        if ligand_hit and receptor_hit:
            return _clean_text(
                sentence
            )

    return None


# --------------------------------------------------
# NCBI request layer
# --------------------------------------------------
def _ncbi_settings() -> dict[str, str]:
    """Email and API key for NCBI, from the web app settings or the environment.

    Both are optional. NCBI E-utilities work without them; an email is
    recommended and an API key raises the rate limit and reduces refusals.
    """
    email = str(os.environ.get("NCBI_EMAIL", "") or "").strip()
    api_key = str(os.environ.get("NCBI_API_KEY", "") or "").strip()
    try:
        from web.backend import settings as _settings  # late import; avoids cycles
        cfg = _settings.load()
        email = email or str(cfg.get("ncbi_email", "") or "").strip()
        api_key = api_key or str(cfg.get("ncbi_api_key", "") or "").strip()
    except Exception:
        pass
    return {"email": email, "api_key": api_key}


def _get_email(
    email: str | None,
) -> str:
    # Email is no longer required; kept for backward compatibility.
    if email:
        return str(email).strip()
    return _ncbi_settings()["email"]


def _http_get(
    endpoint: str,
    params: dict[str, Any],
) -> bytes:
    params = {k: v for k, v in params.items() if v not in (None, "")}
    url = (
        EUTILS_BASE
        + endpoint
        + "?"
        + urlencode(
            params
        )
    )

    request = Request(
        url,
        headers={
            "User-Agent":
                "ccc-platform/1.0 (https://github.com/; mailto:"
                + str(params.get("email") or "anonymous") + ")",
        },
    )

    with urlopen(
        request,
        timeout=15,
    ) as response:
        return response.read()


def _esearch(
    query: str,
    max_results: int,
    email: str,
) -> list[str]:
    payload = _http_get(
        "esearch.fcgi",
        {
            "db":
                "pubmed",

            "term":
                query,

            "retmode":
                "json",

            "retmax":
                int(
                    max_results
                ),

            "sort":
                "relevance",

            "tool":
                TOOL_NAME,

            "email":
                email,

            "api_key":
                _ncbi_settings()["api_key"],
        },
    )

    parsed = json.loads(
        payload.decode(
            "utf-8"
        )
    )

    result = parsed.get("esearchresult", {}) or {}
    return [str(pmid) for pmid in result.get("idlist", [])]


def _extract_year(
    article: ET.Element,
) -> str | None:
    paths = [
        ".//JournalIssue/PubDate/Year",
        ".//ArticleDate/Year",
        ".//PubDate/MedlineDate",
    ]

    for path in paths:
        value = _xml_text(
            article.find(
                path
            )
        )

        if value:
            match = re.search(
                r"\b(19|20)\d{2}\b",
                value,
            )

            if match:
                return match.group(
                    0
                )

            return value

    return None


def _extract_doi(
    article: ET.Element,
) -> str | None:
    for item in article.findall(
        ".//ArticleId"
    ):
        if (
            item.attrib.get(
                "IdType"
            )
            == "doi"
        ):
            value = _xml_text(
                item
            )

            if value:
                return value

    return None


def _extract_abstract(
    article: ET.Element,
) -> str | None:
    sections = []

    for element in article.findall(
        ".//Abstract/AbstractText"
    ):
        text = _xml_text(
            element
        )

        if not text:
            continue

        label = _clean_text(
            element.attrib.get(
                "Label"
            )
        )

        if label:
            sections.append(
                f"{label}: {text}"
            )

        else:
            sections.append(
                text
            )

    if not sections:
        return None

    return " ".join(
        sections
    )


def _extract_authors(
    article: ET.Element,
) -> list[str]:
    authors = []

    for author in article.findall(
        ".//AuthorList/Author"
    ):
        collective = _xml_text(
            author.find(
                "CollectiveName"
            )
        )

        if collective:
            authors.append(
                collective
            )
            continue

        last = _xml_text(
            author.find(
                "LastName"
            )
        )

        initials = _xml_text(
            author.find(
                "Initials"
            )
        )

        if last and initials:
            authors.append(
                f"{last} {initials}"
            )

        elif last:
            authors.append(
                last
            )

    return authors


def _extract_mesh(
    article: ET.Element,
) -> list[str]:
    terms = []

    for heading in article.findall(
        ".//MeshHeadingList/MeshHeading"
    ):
        descriptor = _xml_text(
            heading.find(
                "DescriptorName"
            )
        )

        if descriptor:
            terms.append(
                descriptor
            )

    return sorted(
        set(
            terms
        )
    )


def _extract_publication_types(
    article: ET.Element,
) -> list[str]:
    values = []

    for element in article.findall(
        ".//PublicationTypeList/PublicationType"
    ):
        value = _xml_text(
            element
        )

        if value:
            values.append(
                value
            )

    return sorted(
        set(
            values
        )
    )


def _parse_pubmed_xml(
    xml_bytes: bytes,
) -> list[dict[str, Any]]:
    root = ET.fromstring(
        xml_bytes
    )

    results = []

    for record in root.findall(
        ".//PubmedArticle"
    ):
        pmid = _xml_text(
            record.find(
                ".//MedlineCitation/PMID"
            )
        )

        if not pmid:
            continue

        title = _xml_text(
            record.find(
                ".//Article/ArticleTitle"
            )
        )

        journal = _xml_text(
            record.find(
                ".//Article/Journal/Title"
            )
        )

        results.append(
            {
                "pmid":
                    pmid,

                "title":
                    title,

                "abstract":
                    _extract_abstract(
                        record
                    ),

                "journal":
                    journal,

                "year":
                    _extract_year(
                        record
                    ),

                "authors":
                    _extract_authors(
                        record
                    ),

                "mesh_terms":
                    _extract_mesh(
                        record
                    ),

                "publication_types":
                    _extract_publication_types(
                        record
                    ),

                "doi":
                    _extract_doi(
                        record
                    ),

                "pubmed_url":
                    (
                        "https://pubmed.ncbi.nlm.nih.gov/"
                        + pmid
                        + "/"
                    ),
            }
        )

    return results


def _efetch(
    pmids: list[str],
    email: str,
) -> list[dict[str, Any]]:
    if not pmids:
        return []

    payload = _http_get(
        "efetch.fcgi",
        {
            "db":
                "pubmed",

            "id":
                ",".join(
                    pmids
                ),

            "retmode":
                "xml",

            "tool":
                TOOL_NAME,

            "email":
                email,

            "api_key":
                _ncbi_settings()["api_key"],
        },
    )

    return _parse_pubmed_xml(
        payload
    )


# --------------------------------------------------
# Cache persistence
# --------------------------------------------------
def _save_articles(
    articles: list[dict[str, Any]],
) -> None:
    if not articles:
        return

    init_cache()

    con = _connect()

    for article in articles:
        con.execute(
            """
            INSERT INTO pubmed_articles (
                pmid,
                title,
                abstract,
                journal,
                year,
                authors_json,
                mesh_json,
                publication_types_json,
                doi,
                fetched_at
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(pmid) DO UPDATE SET
                title = excluded.title,
                abstract = excluded.abstract,
                journal = excluded.journal,
                year = excluded.year,
                authors_json = excluded.authors_json,
                mesh_json = excluded.mesh_json,
                publication_types_json =
                    excluded.publication_types_json,
                doi = excluded.doi,
                fetched_at = excluded.fetched_at
            """,
            (
                article[
                    "pmid"
                ],
                article.get(
                    "title"
                ),
                article.get(
                    "abstract"
                ),
                article.get(
                    "journal"
                ),
                article.get(
                    "year"
                ),
                _json_dumps(
                    article.get(
                        "authors",
                        [],
                    )
                ),
                _json_dumps(
                    article.get(
                        "mesh_terms",
                        [],
                    )
                ),
                _json_dumps(
                    article.get(
                        "publication_types",
                        [],
                    )
                ),
                article.get(
                    "doi"
                ),
                _utc_now(),
            ),
        )

    con.commit()
    con.close()


def _load_articles(
    pmids: list[str],
) -> list[dict[str, Any]]:
    if not pmids:
        return []

    init_cache()

    con = _connect()

    placeholders = ",".join(
        "?"
        for _ in pmids
    )

    rows = con.execute(
        f"""
        SELECT *
        FROM pubmed_articles
        WHERE pmid IN ({placeholders})
        """,
        tuple(
            pmids
        ),
    ).fetchall()

    con.close()

    by_pmid = {}

    for row in rows:
        by_pmid[
            str(
                row[
                    "pmid"
                ]
            )
        ] = {
            "pmid":
                str(
                    row[
                        "pmid"
                    ]
                ),

            "title":
                row[
                    "title"
                ],

            "abstract":
                row[
                    "abstract"
                ],

            "journal":
                row[
                    "journal"
                ],

            "year":
                row[
                    "year"
                ],

            "authors":
                _json_loads(
                    row[
                        "authors_json"
                    ],
                    [],
                ),

            "mesh_terms":
                _json_loads(
                    row[
                        "mesh_json"
                    ],
                    [],
                ),

            "publication_types":
                _json_loads(
                    row[
                        "publication_types_json"
                    ],
                    [],
                ),

            "doi":
                row[
                    "doi"
                ],

            "pubmed_url":
                (
                    "https://pubmed.ncbi.nlm.nih.gov/"
                    + str(
                        row[
                            "pmid"
                        ]
                    )
                    + "/"
                ),
        }

    return [
        by_pmid[
            pmid
        ]
        for pmid in pmids
        if pmid in by_pmid
    ]


def _save_pair_search(
    ligand: str,
    receptor: str,
    query: str,
    max_results: int,
    pmids: list[str],
) -> None:
    init_cache()

    con = _connect()

    con.execute(
        """
        INSERT INTO pair_search_cache (
            ligand,
            receptor,
            query,
            max_results,
            pmids_json,
            fetched_at
        )
        VALUES (?, ?, ?, ?, ?, ?)
        ON CONFLICT(
            ligand,
            receptor,
            max_results
        )
        DO UPDATE SET
            query = excluded.query,
            pmids_json = excluded.pmids_json,
            fetched_at = excluded.fetched_at
        """,
        (
            ligand,
            receptor,
            query,
            int(
                max_results
            ),
            _json_dumps(
                pmids
            ),
            _utc_now(),
        ),
    )

    con.commit()
    con.close()


def _load_pair_search(
    ligand: str,
    receptor: str,
    max_results: int,
) -> dict[str, Any] | None:
    init_cache()

    con = _connect()

    row = con.execute(
        """
        SELECT *
        FROM pair_search_cache
        WHERE ligand = ?
          AND receptor = ?
          AND max_results = ?
        LIMIT 1
        """,
        (
            ligand,
            receptor,
            int(
                max_results
            ),
        ),
    ).fetchone()

    con.close()

    if row is None:
        return None

    return {
        "query":
            row[
                "query"
            ],

        "pmids":
            _json_loads(
                row[
                    "pmids_json"
                ],
                [],
            ),

        "fetched_at":
            row[
                "fetched_at"
            ],
    }


# --------------------------------------------------
# Automatic context extraction from PubMed metadata
# --------------------------------------------------
_SPECIES_MESH_MAP = {
    "Humans": "Human",
    "Mice": "Mouse",
    "Rats": "Rat",
    "Rabbits": "Rabbit",
    "Dogs": "Dog",
    "Cats": "Cat",
    "Swine": "Pig",
    "Cattle": "Cattle",
    "Sheep": "Sheep",
    "Goats": "Goat",
    "Zebrafish": "Zebrafish",
    "Drosophila": "Drosophila",
    "Primates": "Primate",
}

_CELL_TERM_PATTERN = re.compile(
    r"(cell|cells|cell line|fibroblast|macrophage|"
    r"monocyte|neutrophil|lymphocyte|t[\s-]?cell|"
    r"b[\s-]?cell|endothelial|epithelial|astrocyte|"
    r"microglia|neuron|hepatocyte|keratinocyte|"
    r"adipocyte|osteoblast|osteoclast|chondrocyte|"
    r"myocyte|myoblast|trophoblast|dendritic|"
    r"mast cell|platelet|erythrocyte|beta cell|β-cell)",
    flags=re.IGNORECASE,
)

_TISSUE_KEYWORDS = [
    "brain", "breast", "lung", "liver", "kidney", "heart",
    "skin", "bone", "bone marrow", "blood", "placenta",
    "endomet", "uterus", "uterine", "ovary", "ovarian",
    "prostate", "pancreas", "pancreatic", "colon",
    "colorectal", "intestin", "stomach", "gastric",
    "thyroid", "adipose", "muscle", "retina", "glioma",
    "glioblastoma", "sarcoma", "carcinoma", "neoplasm",
    "tumor", "tumour",
]

_EXPERIMENT_PATTERNS = [
    ("In vitro", re.compile(r"\bin vitro\b", re.IGNORECASE)),
    ("In vivo", re.compile(r"\bin vivo\b", re.IGNORECASE)),
    ("Ex vivo", re.compile(r"\bex vivo\b", re.IGNORECASE)),
    (
        "Primary cells / tissue",
        re.compile(
            r"\bprimary (?:cells?|tissue|samples?)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "Patient-derived material",
        re.compile(
            r"\bpatient[- ]derived\b|\bpatients?\b",
            re.IGNORECASE,
        ),
    ),
    ("Cell line", re.compile(r"\bcell line\b", re.IGNORECASE)),
    ("Organoid model", re.compile(r"\borganoid", re.IGNORECASE)),
    ("Xenograft model", re.compile(r"\bxenograft", re.IGNORECASE)),
    (
        "Knockout / genetic deletion",
        re.compile(
            r"\bknockout\b|\bknock-out\b|\bdeletion\b",
            re.IGNORECASE,
        ),
    ),
    (
        "RNA sequencing",
        re.compile(
            r"\bRNA[- ]?seq\b|\bRNA sequencing\b",
            re.IGNORECASE,
        ),
    ),
    (
        "Single-cell analysis",
        re.compile(
            r"\bsingle[- ]cell\b|\bscRNA[- ]?seq\b",
            re.IGNORECASE,
        ),
    ),
    (
        "Flow cytometry",
        re.compile(r"\bflow cytometr", re.IGNORECASE),
    ),
    (
        "Immunohistochemistry",
        re.compile(r"\bimmunohistochem", re.IGNORECASE),
    ),
]


def _unique_preserve_order(values: list[str]) -> list[str]:
    seen = set()
    output = []

    for value in values:
        cleaned = _clean_text(value)

        if not cleaned:
            continue

        key = cleaned.lower()

        if key in seen:
            continue

        seen.add(key)
        output.append(cleaned)

    return output


def _extract_organisms(article: dict[str, Any]) -> list[str]:
    organisms = []

    for term in article.get("mesh_terms", []):
        if term in _SPECIES_MESH_MAP:
            organisms.append(_SPECIES_MESH_MAP[term])

    if organisms:
        return _unique_preserve_order(organisms)

    text = " ".join(
        [
            str(article.get("title") or ""),
            str(article.get("abstract") or ""),
        ]
    )

    fallbacks = [
        ("Human", r"\bhuman(?:s)?\b|\bpatients?\b"),
        ("Mouse", r"\bmice\b|\bmouse\b|\bmurine\b"),
        ("Rat", r"\brats?\b"),
        ("Zebrafish", r"\bzebrafish\b"),
    ]

    for label, pattern in fallbacks:
        if re.search(pattern, text, flags=re.IGNORECASE):
            organisms.append(label)

    return _unique_preserve_order(organisms)


def _extract_cell_types(
    article: dict[str, Any],
) -> list[str]:
    """
    Extract likely biological cell types from PubMed metadata.

    We deliberately avoid a generic "contains the word cell"
    rule because MeSH terms such as "Cell Proliferation" or
    "Cell Nucleus" are processes/structures, not cell types.
    """
    candidates = []

    mesh_patterns = [
        r"fibroblast",
        r"macrophage",
        r"monocyte",
        r"neutrophil",
        r"lymphocyte",
        r"t-lymphocyte",
        r"b-lymphocyte",
        r"endothelial",
        r"epithelial",
        r"astrocyte",
        r"microglia",
        r"neuron",
        r"hepatocyte",
        r"keratinocyte",
        r"adipocyte",
        r"osteoblast",
        r"osteoclast",
        r"chondrocyte",
        r"myocyte",
        r"myoblast",
        r"trophoblast",
        r"dendritic cell",
        r"mast cell",
        r"platelet",
        r"erythrocyte",
        r"endocrine cell",
        r"insulin-secreting cell",
        r"tumor cells",
        r"cancer-associated fibroblast",
    ]

    for term in article.get(
        "mesh_terms",
        [],
    ):
        lowered = term.lower()

        if any(
            re.search(
                pattern,
                lowered,
                flags=re.IGNORECASE,
            )
            for pattern
            in mesh_patterns
        ):
            candidates.append(
                term
            )

    text = " ".join(
        [
            str(
                article.get(
                    "title"
                )
                or ""
            ),
            str(
                article.get(
                    "abstract"
                )
                or ""
            ),
        ]
    )

    phrase_patterns = [
        (
            "Cancer-associated fibroblasts",
            r"\bcancer[- ]associated fibroblasts?\b|\bCAFs\b",
        ),
        (
            "Fibroblasts",
            r"\bfibroblasts?\b",
        ),
        (
            "CD8+ T cells",
            r"\bCD8\+?\s*T cells?\b",
        ),
        (
            "CD4+ T cells",
            r"\bCD4\+?\s*T cells?\b",
        ),
        (
            "Macrophages",
            r"\bmacrophages?\b",
        ),
        (
            "Endothelial cells",
            r"\bendothelial cells?\b",
        ),
        (
            "Epithelial cells",
            r"\bepithelial cells?\b",
        ),
        (
            "β-cells",
            r"\b(?:beta|β)[- ]cells?\b",
        ),
    ]

    for label, pattern in phrase_patterns:
        if re.search(
            pattern,
            text,
            flags=re.IGNORECASE,
        ):
            candidates.append(
                label
            )

    return _unique_preserve_order(
        candidates
    )[:12]
def _extract_tissue_context(
    article: dict[str, Any],
) -> list[str]:
    """
    Keep tissue, organ, cancer/disease-site and microenvironment
    terms, while dropping generic molecular/process terms.
    """
    candidates = []

    excluded_terms = {
        "Cell Line, Tumor",
        "Neoplasm Proteins",
        "Drug Resistance, Neoplasm",
        "Tumor Escape",
        "Tumor Burden",
        "Blood Glucose",
    }

    for term in article.get(
        "mesh_terms",
        [],
    ):
        if term in excluded_terms:
            continue

        lowered = term.lower()

        if any(
            keyword in lowered
            for keyword
            in _TISSUE_KEYWORDS
        ):
            candidates.append(
                term
            )

    return _unique_preserve_order(
        candidates
    )[:12]
def _extract_experimental_context(
    article: dict[str, Any],
) -> list[str]:
    values = []

    for publication_type in article.get("publication_types", []):
        if publication_type in {
            "Clinical Trial",
            "Clinical Trial, Phase I",
            "Clinical Trial, Phase II",
            "Clinical Trial, Phase III",
            "Clinical Trial, Phase IV",
            "Randomized Controlled Trial",
            "Observational Study",
            "Comparative Study",
            "Case Reports",
        }:
            values.append(publication_type)

    for term in article.get("mesh_terms", []):
        if term in {
            "Cell Line",
            "Cell Line, Tumor",
            "Cells, Cultured",
            "Animal Experimentation",
            "Models, Animal",
            "Disease Models, Animal",
            "Organoids",
            "Patient-Derived Xenograft Models",
            "Gene Knockdown Techniques",
            "Gene Expression Profiling",
            "RNA-Seq",
            "Single-Cell Analysis",
        }:
            values.append(term)

    text = " ".join(
        [
            str(article.get("title") or ""),
            str(article.get("abstract") or ""),
        ]
    )

    for label, pattern in _EXPERIMENT_PATTERNS:
        if pattern.search(text):
            values.append(label)

    return _unique_preserve_order(values)[:12]


def _term_pattern(
    term: str,
) -> re.Pattern:
    escaped = re.escape(
        str(
            term
        ).strip()
    )

    # Gene/protein symbols such as IGF2 must not match
    # longer symbols such as IGF2R or IGF2BP1.
    return re.compile(
        rf"(?<![A-Za-z0-9]){escaped}(?![A-Za-z0-9])",
        flags=re.IGNORECASE,
    )


def _extract_pair_snippet(
    article: dict[str, Any],
    ligand: str,
    receptor: str,
) -> str | None:
    abstract = _clean_text(
        article.get(
            "abstract"
        )
    )

    if not abstract:
        return None

    ligand_terms = [
        ligand,
    ]

    receptor_terms = [
        receptor,
    ]

    ligand_alias = (
        _protein_name_alias(
            ligand
        )
    )

    receptor_alias = (
        _protein_name_alias(
            receptor
        )
    )

    if ligand_alias:
        ligand_terms.append(
            ligand_alias
        )

    if receptor_alias:
        receptor_terms.append(
            receptor_alias
        )

    ligand_patterns = [
        _term_pattern(
            term
        )
        for term in ligand_terms
        if term
    ]

    receptor_patterns = [
        _term_pattern(
            term
        )
        for term in receptor_terms
        if term
    ]

    sentences = re.split(
        r"(?<=[.!?])\s+",
        abstract,
    )

    for sentence in sentences:
        ligand_hit = any(
            pattern.search(
                sentence
            )
            for pattern
            in ligand_patterns
        )

        receptor_hit = any(
            pattern.search(
                sentence
            )
            for pattern
            in receptor_patterns
        )

        if (
            ligand_hit
            and receptor_hit
        ):
            return _clean_text(
                sentence
            )

    return None
def _add_context_summary(
    article: dict[str, Any],
    ligand: str,
    receptor: str,
) -> dict[str, Any]:
    enriched = dict(article)

    enriched["organisms"] = _extract_organisms(article)
    enriched["tissue_context"] = _extract_tissue_context(article)
    enriched["cell_types"] = _extract_cell_types(article)
    enriched["experimental_context"] = (
        _extract_experimental_context(article)
    )
    enriched["pair_evidence_snippet"] = _extract_pair_snippet(
        article,
        ligand=ligand,
        receptor=receptor,
    )
    enriched["context_extraction_note"] = (
        "Automatically extracted from PubMed MeSH terms, "
        "publication metadata and abstract text. "
        "Verify against the original publication."
    )

    return enriched

# --------------------------------------------------
# Public interface
# --------------------------------------------------
def _relaxed_pair_query(query_metadata: dict[str, Any]) -> str | None:
    """A looser query: require the ligand and at least ONE receptor subunit
    (and at least one ligand subunit if the ligand is a complex), instead of
    every subunit in the same abstract. Returns None if nothing to relax."""
    def side_block(meta: dict[str, Any]) -> str:
        subs = list(meta.get("required_subunits") or [])
        terms = [s for s in subs if s]
        if not terms:
            return ""
        return "(" + " OR ".join(
            f'"{_escape_pubmed_phrase(t)}"[Title/Abstract]' for t in terms
        ) + ")"
    lig = side_block(query_metadata.get("ligand") or {})
    rcp = side_block(query_metadata.get("receptor") or {})
    if not lig or not rcp:
        return None
    return lig + " AND " + rcp


def search_pubmed_pair_entity_aware(
    ligand: str,
    receptor: str,
    *,
    ligand_entity_key: str | None = None,
    receptor_entity_key: str | None = None,
    ligand_required_subunits: list[str] | tuple[str, ...] | None = None,
    receptor_required_subunits: list[str] | tuple[str, ...] | None = None,
    max_results: int = DEFAULT_MAX_RESULTS,
    email: str | None = None,
    force_refresh: bool = False,
) -> dict[str, Any]:
    """
    Search PubMed using molecular-entity-aware query construction.

    Verified complexes are not collapsed to a single display label.
    Instead, the query requires the documented complex subunits.
    """
    ligand_display = _clean_symbol(
        ligand
    )

    receptor_display = _clean_symbol(
        receptor
    )

    if not ligand_display or not receptor_display:
        raise ValueError(
            "Ligand and receptor are required."
        )

    max_results = max(
        1,
        min(
            int(
                max_results
            ),
            25,
        ),
    )

    query, query_metadata = build_entity_aware_pair_query(
        ligand_display,
        receptor_display,
        ligand_entity_key=ligand_entity_key,
        receptor_entity_key=receptor_entity_key,
        ligand_required_subunits=ligand_required_subunits,
        receptor_required_subunits=receptor_required_subunits,
    )

    ligand_cache_key = _entity_cache_key(
        ligand_display,
        ligand_entity_key,
        ligand_required_subunits,
    )

    receptor_cache_key = _entity_cache_key(
        receptor_display,
        receptor_entity_key,
        receptor_required_subunits,
    )

    init_cache()

    cached = None

    if not force_refresh:
        cached = _load_pair_search(
            ligand_cache_key,
            receptor_cache_key,
            max_results,
        )

    if cached is not None:
        articles = _load_articles(
            cached[
                "pmids"
            ]
        )
        fetched_at = cached[
            "fetched_at"
        ]
        cache_status = True

    else:
        resolved_email = _get_email(
            email
        )

        used_query = query
        query_mode = "strict"
        pmids = _esearch(
            query,
            max_results,
            resolved_email,
        )

        is_complex_pair = (
            query_metadata["ligand"].get("entity_type") == "complex"
            or query_metadata["receptor"].get("entity_type") == "complex"
        )
        if not pmids and is_complex_pair:
            relaxed = _relaxed_pair_query(query_metadata)
            if relaxed and relaxed != query:
                relaxed_pmids = _esearch(relaxed, max_results, resolved_email)
                if relaxed_pmids:
                    pmids = relaxed_pmids
                    used_query = relaxed
                    query_mode = "relaxed"

        articles = _efetch(
            pmids,
            resolved_email,
        )

        _save_articles(
            articles
        )

        _save_pair_search(
            ligand_cache_key,
            receptor_cache_key,
            used_query,
            max_results,
            pmids,
        )

        query = used_query
        query_metadata = {**query_metadata, "query_mode": query_mode}
        fetched_at = _utc_now()
        cache_status = False

    ligand_terms = []

    if query_metadata["ligand"].get(
        "entity_type"
    ) == "complex":
        for subunit in query_metadata["ligand"].get(
            "required_subunits",
            [],
        ):
            ligand_terms.extend(
                _pubmed_terms_for_gene(
                    subunit
                )
            )
    else:
        ligand_terms = _pubmed_terms_for_gene(
            ligand_display
        )

    receptor_term_groups = []

    if query_metadata["receptor"].get(
        "entity_type"
    ) == "complex":
        for subunit in query_metadata["receptor"].get(
            "required_subunits",
            [],
        ):
            receptor_term_groups.append(
                _pubmed_terms_for_gene(
                    subunit
                )
            )
    else:
        receptor_term_groups.append(
            _pubmed_terms_for_gene(
                receptor_display
            )
        )

    enriched_articles = []

    for article in articles:
        enriched = _add_context_summary(
            article,
            ligand=ligand_display,
            receptor=receptor_display,
        )

        complex_snippet = _extract_entity_pair_snippet(
            article,
            ligand_terms=ligand_terms,
            receptor_term_groups=receptor_term_groups,
        )

        if complex_snippet:
            enriched[
                "pair_evidence_snippet"
            ] = complex_snippet

        enriched[
            "pubmed_query_scope"
        ] = (
            "complex-aware"
            if (
                query_metadata["ligand"].get("entity_type") == "complex"
                or query_metadata["receptor"].get("entity_type") == "complex"
            )
            else "gene"
        )

        enriched_articles.append(
            enriched
        )

    return {
        "ligand":
            ligand_display,

        "receptor":
            receptor_display,

        "query":
            query,

        "query_metadata":
            query_metadata,

        "ligand_cache_key":
            ligand_cache_key,

        "receptor_cache_key":
            receptor_cache_key,

        "articles":
            enriched_articles,

        "result_count":
            len(
                enriched_articles
            ),

        "cached":
            cache_status,

        "fetched_at":
            fetched_at,
    }


def search_pubmed_pair(
    ligand: str,
    receptor: str,
    *,
    max_results: int = DEFAULT_MAX_RESULTS,
    email: str | None = None,
    force_refresh: bool = False,
) -> dict[str, Any]:
    """
    Search PubMed for publications mentioning both members
    of a ligand-receptor pair.

    Runtime behaviour:
    - cached pair -> local SQLite only, no network request
    - uncached pair -> one ESearch + one batched EFetch request
    """

    ligand = _clean_symbol(
        ligand
    )

    receptor = _clean_symbol(
        receptor
    )

    if not ligand or not receptor:
        raise ValueError(
            "Ligand and receptor are required."
        )

    max_results = max(
        1,
        min(
            int(
                max_results
            ),
            25,
        ),
    )

    init_cache()

    if not force_refresh:
        cached_search = (
            _load_pair_search(
                ligand,
                receptor,
                max_results,
            )
        )

        if cached_search is not None:
            articles = (
                _load_articles(
                    cached_search[
                        "pmids"
                    ]
                )
            )

            articles = [
                _add_context_summary(
                    article,
                    ligand=ligand,
                    receptor=receptor,
                )
                for article
                in articles
            ]

            return {
                "ligand":
                    ligand,

                "receptor":
                    receptor,

                "query":
                    cached_search[
                        "query"
                    ],

                "articles":
                    articles,

                "result_count":
                    len(
                        articles
                    ),

                "cached":
                    True,

                "fetched_at":
                    cached_search[
                        "fetched_at"
                    ],
            }

    resolved_email = (
        _get_email(
            email
        )
    )

    query = build_pair_query(
        ligand,
        receptor,
    )

    pmids = _esearch(
        query,
        max_results,
        resolved_email,
    )

    articles = _efetch(
        pmids,
        resolved_email,
    )

    _save_articles(
        articles
    )

    articles = [
        _add_context_summary(
            article,
            ligand=ligand,
            receptor=receptor,
        )
        for article
        in articles
    ]

    _save_pair_search(
        ligand,
        receptor,
        query,
        max_results,
        pmids,
    )

    return {
        "ligand":
            ligand,

        "receptor":
            receptor,

        "query":
            query,

        "articles":
            articles,

        "result_count":
            len(
                articles
            ),

        "cached":
            False,

        "fetched_at":
            _utc_now(),
    }
