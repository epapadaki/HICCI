from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sqlite3
from typing import Any, Literal

from pydantic import BaseModel, Field


PROJECT_ROOT = Path(__file__).resolve().parents[1]

CACHE_DB = (
    PROJECT_ROOT
    / "data"
    / "literature"
    / "publication_evidence_cache.sqlite"
)

DEFAULT_MODEL = os.environ.get(
    "GEMINI_MODEL",
    "gemini-2.5-flash",
)

PROMPT_VERSION = "publication-evidence-classification-v3"


# --------------------------------------------------
# Controlled vocabulary
# --------------------------------------------------
EvidenceType = Literal[
    "direct_binding",
    "perturbation_evidence",
    "functional_response",
    "expression_or_coexpression",
    "computational_prediction",
    "pair_co_mention_only",
    "insufficient_information",
]

PairSpecificStatus = Literal[
    "supported",
    "not_supported",
    "unclear",
]

DirectBindingStatus = Literal[
    "supported_in_abstract",
    "not_supported_in_abstract",
    "unclear",
]

StudyMode = Literal[
    "experimental",
    "computational",
    "mixed",
    "unclear",
]

EvidenceScope = Literal[
    "specific_ligand_receptor_pair",
    "ligand_or_receptor_general",
    "co_mention_only",
    "unclear",
]

ExtractionCertainty = Literal[
    "high",
    "moderate",
    "low",
]

PathwayEffectDirection = Literal[
    "activation",
    "inhibition",
    "modulation",
    "association",
    "unclear",
]

PathwayEvidenceScope = Literal[
    "specific_ligand_receptor_pair",
    "receptor_general",
    "ligand_general",
    "unclear",
]


class PairSpecificPathwayEvidence(BaseModel):
    pathway_name: str = Field(
        description=(
            "Name of the intracellular pathway or downstream "
            "response described in the supplied abstract."
        )
    )

    pathway_description: str = Field(
        description=(
            "Short description of the pathway or response, based "
            "only on the supplied abstract/metadata."
        )
    )

    relationship_to_pair: str = Field(
        description=(
            "How the supplied abstract links the pathway or response "
            "to the candidate ligand-receptor pair. State explicitly "
            "whether the link is pair-specific or indirect."
        )
    )

    effect_direction: PathwayEffectDirection

    evidence_scope: PathwayEvidenceScope

    biological_context: str = Field(
        description=(
            "Cell type, tissue, organism, disease/model or experimental "
            "context for this pathway statement, when available."
        )
    )

    source_sentence: str = Field(
        description=(
            "Short supporting abstract sentence or paraphrased clause. "
            "Do not invent full-text evidence."
        )
    )

    limitations: list[str] = Field(
        default_factory=list,
        description=(
            "Limitations of this pathway extraction, especially whether "
            "it is abstract-only, indirect, receptor-general, or unclear."
        ),
    )


class PublicationEvidenceClassification(BaseModel):
    pmid: str

    primary_evidence_type: EvidenceType

    evidence_types: list[
        EvidenceType
    ] = Field(
        description=(
            "All evidence types explicitly supported by the "
            "available PubMed metadata/abstract. Do not add "
            "categories based on outside knowledge."
        )
    )

    pair_specific_evidence: PairSpecificStatus

    direct_binding_evidence: DirectBindingStatus

    study_mode: StudyMode

    evidence_scope: EvidenceScope

    what_was_studied: str = Field(
        description=(
            "Short description of what the study actually "
            "tested or analysed, based only on the supplied "
            "metadata and abstract."
        )
    )

    finding_summary: str = Field(
        description=(
            "Short, cautious paraphrase of the pair-relevant "
            "finding supported by the available abstract."
        )
    )

    experimental_context: str = Field(
        description=(
            "Organism, tissue, cell system, disease/condition "
            "or model only when available from the supplied "
            "publication information."
        )
    )

    classification_reason: str = Field(
        description=(
            "Why this evidence category was chosen, explicitly "
            "distinguishing pair-specific evidence from simple "
            "co-mention."
        )
    )

    limitations: list[str]

    extraction_certainty: ExtractionCertainty

    pair_specific_pathways: list[PairSpecificPathwayEvidence] = Field(
        default_factory=list,
        description=(
            "Pathway or downstream-response statements that the supplied "
            "abstract explicitly links to the candidate ligand-receptor "
            "pair. Leave empty when the abstract only gives receptor-"
            "general pathways, ligand-general pathways, simple co-mention, "
            "or insufficient information."
        ),
    )


class ClassificationResponse(BaseModel):
    classifications: list[
        PublicationEvidenceClassification
    ]


# --------------------------------------------------
# Helpers
# --------------------------------------------------
def _utc_now() -> str:
    return datetime.now(
        timezone.utc
    ).isoformat(
        timespec="seconds"
    )


def _clean_text(
    value: Any,
) -> str:
    if value is None:
        return ""

    return " ".join(
        str(
            value
        ).split()
    ).strip()


def _get_api_key() -> str:
    key = (
        os.environ.get(
            "GEMINI_API_KEY",
            "",
        )
        .strip()
    )

    if not key:
        raise RuntimeError(
            "GEMINI_API_KEY is not configured in the "
            "server environment."
        )

    return key


def _json_dumps(
    value: Any,
) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
    )


def _article_payload(
    article: dict[str, Any],
) -> dict[str, Any]:
    abstract = _clean_text(
        article.get(
            "abstract"
        )
    )

    # Keep each request bounded while preserving a substantial
    # amount of abstract text.
    if len(
        abstract
    ) > 7000:
        abstract = (
            abstract[:7000]
            + "..."
        )

    return {
        "pmid":
            str(
                article.get(
                    "pmid"
                )
                or ""
            ),

        "title":
            _clean_text(
                article.get(
                    "title"
                )
            ),

        "year":
            _clean_text(
                article.get(
                    "year"
                )
            ),

        "journal":
            _clean_text(
                article.get(
                    "journal"
                )
            ),

        "organisms":
            article.get(
                "organisms"
            )
            or [],

        "tissue_context":
            article.get(
                "tissue_context"
            )
            or [],

        "cell_types":
            article.get(
                "cell_types"
            )
            or [],

        "experimental_context":
            article.get(
                "experimental_context"
            )
            or [],

        "pair_specific_abstract_context":
            _clean_text(
                article.get(
                    "pair_evidence_snippet"
                )
            ),

        "abstract":
            abstract,
    }


def _content_hash(
    *,
    ligand: str,
    receptor: str,
    article: dict[str, Any],
) -> str:
    payload = {
        "ligand":
            ligand,

        "receptor":
            receptor,

        "article":
            _article_payload(
                article
            ),

        "prompt_version":
            PROMPT_VERSION,
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


# --------------------------------------------------
# Local SQLite cache
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

    con.execute(
        """
        CREATE TABLE IF NOT EXISTS publication_evidence (
            content_hash TEXT NOT NULL,
            pmid TEXT NOT NULL,
            ligand TEXT NOT NULL,
            receptor TEXT NOT NULL,
            model TEXT NOT NULL,
            prompt_version TEXT NOT NULL,
            result_json TEXT NOT NULL,
            created_at TEXT NOT NULL,
            PRIMARY KEY (
                content_hash,
                model,
                prompt_version
            )
        )
        """
    )

    con.commit()
    con.close()


def _load_cached(
    *,
    content_hash: str,
    model: str,
) -> dict[str, Any] | None:
    init_cache()

    con = _connect()

    row = con.execute(
        """
        SELECT
            result_json,
            created_at
        FROM publication_evidence
        WHERE content_hash = ?
          AND model = ?
          AND prompt_version = ?
        LIMIT 1
        """,
        (
            content_hash,
            model,
            PROMPT_VERSION,
        ),
    ).fetchone()

    con.close()

    if row is None:
        return None

    result = json.loads(
        row[
            "result_json"
        ]
    )

    result[
        "cached"
    ] = True

    result[
        "cached_at"
    ] = row[
        "created_at"
    ]

    return result


def _save_cached(
    *,
    content_hash: str,
    pmid: str,
    ligand: str,
    receptor: str,
    model: str,
    result: dict[str, Any],
) -> None:
    init_cache()

    stored = dict(
        result
    )

    stored.pop(
        "cached",
        None,
    )

    stored.pop(
        "cached_at",
        None,
    )

    con = _connect()

    con.execute(
        """
        INSERT INTO publication_evidence (
            content_hash,
            pmid,
            ligand,
            receptor,
            model,
            prompt_version,
            result_json,
            created_at
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(
            content_hash,
            model,
            prompt_version
        )
        DO UPDATE SET
            result_json = excluded.result_json,
            created_at = excluded.created_at
        """,
        (
            content_hash,
            pmid,
            ligand,
            receptor,
            model,
            PROMPT_VERSION,
            _json_dumps(
                stored
            ),
            _utc_now(),
        ),
    )

    con.commit()
    con.close()


# --------------------------------------------------
# Scientific classification prompt
# --------------------------------------------------
SYSTEM_INSTRUCTIONS = """
You are classifying the type of biomedical evidence available
for a candidate ligand-receptor publication.

Use ONLY the publication metadata and abstract supplied in the
input. Do not use external biological knowledge. Do not infer
results from the title alone when the abstract does not support
them.

The PubMed records were found because they mention both pair
members. Co-mention alone is NOT evidence that the ligand binds
the receptor or that the pair mediates signaling.

Classify conservatively.

CONTROLLED EVIDENCE TYPES

direct_binding
- Use only when the available abstract explicitly supports a
  physical ligand-receptor interaction using direct molecular
  binding evidence, a binding/affinity assay, structural
  evidence, or equally explicit physical interaction evidence.
- Receptor dependence, pathway activation, co-expression,
  correlation, docking or computational prediction are NOT
  direct binding evidence by themselves.

perturbation_evidence
- The abstract describes a perturbation of the ligand, receptor,
  or their interaction (for example knockout, knockdown,
  blocking antibody, inhibitor, overexpression or rescue) and
  the result specifically supports a role for the candidate
  pair or one component in the reported response.
- Do not call this pair-specific unless the abstract actually
  links the perturbation to the ligand-receptor relationship.

functional_response
- Use only when the available abstract provides experimental
  evidence that explicitly links the candidate ligand to the
  candidate receptor in a functional response.
- Examples that can qualify include ligand manipulation causing
  a receptor-dependent response, receptor perturbation altering
  a ligand-induced response, or an experiment that directly
  tests the candidate ligand-receptor axis and reports a
  functional consequence.
- A computationally predicted ligand-receptor pair followed by
  an experiment showing a change in the receptor, the ligand,
  or the same cell population does NOT by itself qualify as
  pair-specific functional response evidence.
- General downstream effects, co-culture effects, expression
  changes, pathway changes or phenotypic changes should not be
  upgraded to functional_response unless the abstract explicitly
  connects those effects to the candidate ligand-receptor pair.

expression_or_coexpression
- The evidence is primarily expression, abundance,
  co-expression, localization or association of the two
  molecules/cell populations, without a demonstrated functional
  pair relationship.

computational_prediction
- The candidate relationship is inferred or predicted
  computationally and the supplied abstract does not provide
  experimental validation of the pair.

pair_co_mention_only
- Both pair members appear in the publication/abstract, but the
  supplied abstract does not explicitly establish a biological
  relationship between them.

insufficient_information
- The supplied information is too limited or ambiguous to make
  a defensible classification.

PAIR-SPECIFIC STATUS
- supported: the supplied abstract explicitly connects the
  ligand and receptor as a candidate pair in the same mechanism,
  experiment, computational inference or functional relationship.
- not_supported: both are present but no pair-specific
  relationship is established in the supplied abstract.
- unclear: the wording is too ambiguous to decide.

Important distinction:
- pair_specific_evidence = supported does NOT automatically mean
  experimental validation.
- A pair can be explicitly predicted computationally and still
  have pair_specific_evidence = supported while the primary
  evidence type remains computational_prediction.

DIRECT BINDING STATUS
- supported_in_abstract: the abstract explicitly provides direct
  physical binding evidence.
- not_supported_in_abstract: such evidence is not present in the
  supplied abstract. This does NOT mean the full paper proves
  that binding does not occur.
- unclear: the abstract wording is ambiguous.

EVIDENCE SCOPE
- specific_ligand_receptor_pair: evidence explicitly concerns
  the candidate ligand-receptor relationship.
- ligand_or_receptor_general: evidence concerns one component
  more generally, not the specific pair.
- co_mention_only: both are mentioned without demonstrated pair
  relationship.
- unclear: cannot determine.

STUDY MODE
- experimental, computational, mixed or unclear.

When several evidence types are supported, include all of them
in evidence_types and choose the most pair-informative type as
primary_evidence_type.

Primary evidence must reflect evidence for the candidate pair,
not merely the strongest experiment elsewhere in the paper.

If the pair itself is supported only by computational inference,
while a separate experiment validates a broader biological
context, expression change, cell-state change or one component
alone, keep computational_prediction as the primary evidence
type for the pair.

Use functional_response as primary only when the experimental
functional evidence itself is explicitly pair-specific.

Direct binding should only be primary when direct physical
evidence is explicit.

PAIR-SPECIFIC PATHWAY EXTRACTION
- Extract pair_specific_pathways only when the supplied abstract explicitly
  links the candidate ligand-receptor relationship to a named intracellular
  pathway, signaling cascade, downstream gene program, or functional response.
- Do NOT assign every receptor-general pathway to the candidate ligand.
- If the abstract only says that the receptor is generally associated with
  a pathway, use evidence_scope = receptor_general and include it only if the
  candidate pair is also explicitly discussed; otherwise leave
  pair_specific_pathways empty.
- If the effect is described, choose activation, inhibition, modulation, or
  association. Use unclear when the direction is not supported.
- The pathway extraction is abstract-limited. Do not infer activation in the
  user's dataset from a publication, CellChat category, SIGNOR, Reactome, or
  a CCI score.

Never say that the interaction is validated in the user's
dataset. Never convert CCI scores into biological truth.

The finding_summary must be a cautious paraphrase, not a quote.
If only an abstract is available, limitations should state that
the classification is limited to PubMed metadata/abstract and
may not capture evidence present only in the full text.
""".strip()


def _build_prompt(
    *,
    ligand: str,
    receptor: str,
    articles: list[dict[str, Any]],
) -> str:
    payload = {
        "candidate_pair": {
            "ligand":
                ligand,

            "receptor":
                receptor,
        },

        "publications": [
            _article_payload(
                article
            )
            for article
            in articles
        ],
    }

    return (
        SYSTEM_INSTRUCTIONS
        + "\n\n"
        + "Return exactly one classification for every supplied "
          "PMID and preserve each PMID exactly.\n\n"
        + "INPUT DATA:\n"
        + _json_dumps(
            payload
        )
    )


# --------------------------------------------------
# Gemini request
# --------------------------------------------------
def _gemini_request(
    *,
    ligand: str,
    receptor: str,
    articles: list[dict[str, Any]],
    model: str,
) -> list[dict[str, Any]]:
    """Classify the publications with Gemini (JSON output, validated against
    the ClassificationResponse schema). Uses the shared client in
    utils/gemini_client.py; the key comes from GEMINI_API_KEY (.env)."""
    from .gemini_client import generate_json

    system = (
        "You classify the evidence that biomedical publications provide for a "
        "ligand-receptor pair. Use only the supplied PubMed metadata and abstract. "
        "Return valid JSON only, matching this JSON schema exactly:\n"
        + _json_dumps(ClassificationResponse.model_json_schema())
    )
    data, _model_used = generate_json(
        system,
        _build_prompt(ligand=ligand, receptor=receptor, articles=articles),
        temperature=0.0,
    )
    try:
        parsed = ClassificationResponse.model_validate(data)
    except Exception as error:
        raise RuntimeError(
            "Gemini returned output that did not match "
            "the publication-evidence schema."
        ) from error

    return [
        item.model_dump()
        for item
        in parsed.classifications
    ]


# --------------------------------------------------
# Local rule-based fallback classifier
# --------------------------------------------------
DIRECT_BINDING_KEYWORDS = [
    "bind", "binds", "binding", "bound", "affinity", "ligand binding",
    "receptor binding", "surface plasmon resonance", "spr", "kd",
    "co-immunoprecipitation", "coimmunoprecipitation", "crystal structure",
]

PERTURBATION_KEYWORDS = [
    "knockout", "knock-down", "knockdown", "silencing", "sirna",
    "blocking antibody", "neutralizing antibody", "inhibitor", "inhibited",
    "overexpression", "overexpressed", "rescue", "antagonist", "agonist",
]

FUNCTIONAL_RESPONSE_KEYWORDS = [
    "activates", "activated", "activation", "induces", "induced", "stimulates",
    "stimulated", "promotes", "increases", "decreases", "suppresses",
    "response", "signaling", "phosphorylation", "migration", "proliferation",
    "apoptosis", "differentiation", "chemotaxis", "secretion", "cytokine",
]

EXPRESSION_KEYWORDS = [
    "expression", "expressed", "coexpression", "co-expression", "abundance",
    "localization", "localized", "transcript", "mrna", "rna-seq", "single-cell",
]

COMPUTATIONAL_KEYWORDS = [
    "predicted", "prediction", "computational", "in silico", "network",
    "cellchat", "cellphonedb", "natmi", "ligand-receptor analysis",
]

PATHWAY_TERMS = [
    "MAPK", "ERK", "JNK", "p38", "PI3K", "AKT", "mTOR", "NF-kB",
    "NFκB", "JAK", "STAT", "STAT3", "STAT5", "SMAD", "TGF-beta",
    "TGFβ", "WNT", "beta-catenin", "Notch", "Hedgehog", "Hippo", "YAP",
    "TAZ", "p53", "apoptosis", "proliferation", "migration", "invasion",
    "EMT", "angiogenesis", "inflammation", "chemotaxis", "cytokine",
    "interferon", "integrin", "VEGF", "hypoxia", "fibrosis",
]


def _contains_any(text: str, keywords: list[str]) -> bool:
    lowered = text.lower()
    return any(keyword.lower() in lowered for keyword in keywords)


def _first_matching_terms(text: str, terms: list[str], limit: int = 6) -> list[str]:
    lowered = text.lower()
    seen: list[str] = []
    for term in terms:
        if term.lower() in lowered and term not in seen:
            seen.append(term)
        if len(seen) >= limit:
            break
    return seen


def _short_context(text: str, term: str, max_chars: int = 260) -> str:
    if not text:
        return ""
    lowered = text.lower()
    index = lowered.find(term.lower()) if term else -1
    if index < 0:
        return text[:max_chars].strip()
    start = max(0, index - max_chars // 2)
    end = min(len(text), index + max_chars // 2)
    return text[start:end].strip()


def _pair_terms_present(text: str, ligand: str, receptor: str) -> bool:
    upper = text.upper()
    ligand_clean = _clean_text(ligand).upper()
    receptor_clean = _clean_text(receptor).upper()
    return bool(ligand_clean and receptor_clean and ligand_clean in upper and receptor_clean in upper)


def _local_classify_one_article(
    *,
    ligand: str,
    receptor: str,
    article: dict[str, Any],
) -> dict[str, Any]:
    pmid = str(article.get("pmid") or "").strip()
    title = _clean_text(article.get("title"))
    abstract = _clean_text(article.get("abstract"))
    snippet = _clean_text(article.get("pair_evidence_snippet"))
    text = " ".join(part for part in [title, snippet, abstract] if part)

    has_pair_terms = _pair_terms_present(text, ligand, receptor)
    has_pair_snippet = bool(snippet)

    direct = _contains_any(text, DIRECT_BINDING_KEYWORDS)
    perturbation = _contains_any(text, PERTURBATION_KEYWORDS)
    functional = _contains_any(text, FUNCTIONAL_RESPONSE_KEYWORDS)
    expression = _contains_any(text, EXPRESSION_KEYWORDS)
    computational = _contains_any(text, COMPUTATIONAL_KEYWORDS)

    if direct and (has_pair_terms or has_pair_snippet):
        primary = "direct_binding"
    elif perturbation and (has_pair_terms or has_pair_snippet):
        primary = "perturbation_evidence"
    elif functional and (has_pair_terms or has_pair_snippet):
        primary = "functional_response"
    elif expression:
        primary = "expression_or_coexpression"
    elif computational:
        primary = "computational_prediction"
    elif has_pair_terms or has_pair_snippet:
        primary = "pair_co_mention_only"
    else:
        primary = "insufficient_information"

    evidence_types: list[str] = []
    for flag, name in [
        (direct and (has_pair_terms or has_pair_snippet), "direct_binding"),
        (perturbation and (has_pair_terms or has_pair_snippet), "perturbation_evidence"),
        (functional and (has_pair_terms or has_pair_snippet), "functional_response"),
        (expression, "expression_or_coexpression"),
        (computational, "computational_prediction"),
    ]:
        if flag and name not in evidence_types:
            evidence_types.append(name)
    if not evidence_types:
        evidence_types = [primary]

    pair_specific = "supported" if has_pair_snippet and (direct or perturbation or functional or computational) else ("unclear" if has_pair_terms or has_pair_snippet else "not_supported")
    scope = "specific_ligand_receptor_pair" if pair_specific == "supported" else ("co_mention_only" if has_pair_terms or has_pair_snippet else "unclear")
    direct_status = "supported_in_abstract" if direct and (has_pair_terms or has_pair_snippet) else "not_supported_in_abstract"

    if computational and (direct or perturbation or functional or expression):
        study_mode = "mixed"
    elif computational:
        study_mode = "computational"
    elif direct or perturbation or functional or expression:
        study_mode = "experimental"
    else:
        study_mode = "unclear"

    context_bits = []
    for key in ["organisms", "tissue_context", "cell_types", "experimental_context"]:
        value = article.get(key) or []
        if isinstance(value, list):
            context_bits.extend(str(v).strip() for v in value if str(v).strip())
        elif str(value).strip():
            context_bits.append(str(value).strip())
    experimental_context = "; ".join(context_bits[:8]) or "Not specified in PubMed metadata/abstract."

    pathway_items = []
    pathway_terms = _first_matching_terms(text, PATHWAY_TERMS)
    if pathway_terms and pair_specific == "supported":
        for term in pathway_terms:
            if _contains_any(_short_context(text, term), ["inhibit", "inhibited", "suppress", "blocked", "decrease"]):
                direction = "inhibition"
            elif _contains_any(_short_context(text, term), ["activate", "activated", "induce", "stimulate", "increase", "phosphorylation"]):
                direction = "activation"
            elif _contains_any(_short_context(text, term), ["modulate", "regulate", "alter"]):
                direction = "modulation"
            else:
                direction = "association"
            pathway_items.append({
                "pathway_name": term,
                "pathway_description": f"The abstract/metadata mentions {term} in the same publication context as the selected pair.",
                "relationship_to_pair": "Rule-based extraction: pathway term appears in an abstract/metadata record that also contains pair-specific evidence signals. This is not proof of pathway activation in the user's dataset.",
                "effect_direction": direction,
                "evidence_scope": "specific_ligand_receptor_pair",
                "biological_context": experimental_context,
                "source_sentence": _short_context(text, term),
                "limitations": [
                    "Rule-based extraction from PubMed metadata/abstract only.",
                    "The full text may contain stronger or weaker evidence than the abstract.",
                    "This does not validate pathway activity in the user's dataset.",
                ],
            })

    return {
        "pmid": pmid,
        "primary_evidence_type": primary,
        "evidence_types": evidence_types,
        "pair_specific_evidence": pair_specific,
        "direct_binding_evidence": direct_status,
        "study_mode": study_mode,
        "evidence_scope": scope,
        "what_was_studied": snippet or title or "The PubMed record does not provide enough detail for a stronger local classification.",
        "finding_summary": (
            "Local rule-based classification found pair-specific evidence signals in the PubMed abstract/metadata."
            if pair_specific == "supported"
            else "Local rule-based classification did not find enough abstract-level evidence to treat this as a specific validated ligand-receptor mechanism."
        ),
        "experimental_context": experimental_context,
        "classification_reason": "Generated by deterministic local keyword/context rules, without Gemini. Treat as screening support, not as a curated manual review.",
        "limitations": [
            "Rule-based fallback; no LLM or full-text review used.",
            "Limited to PubMed metadata/abstract available to the platform.",
            "Co-mention and keyword proximity can overestimate or underestimate evidence.",
        ],
        "extraction_certainty": "moderate" if pair_specific == "supported" else "low",
        "pair_specific_pathways": pathway_items,
        "cached": False,
    }


def _local_classify_articles(
    *,
    ligand: str,
    receptor: str,
    articles: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    return [
        _local_classify_one_article(
            ligand=ligand,
            receptor=receptor,
            article=article,
        )
        for article in articles
        if str(article.get("pmid") or "").strip()
    ]


# --------------------------------------------------
# Public interface
# --------------------------------------------------
def classify_publication_evidence(
    ligand: str,
    receptor: str,
    articles: list[dict[str, Any]],
    *,
    model: str | None = None,
    force_refresh: bool = False,
    use_llm: bool = False,
) -> dict[str, Any]:
    """
    Classify what type of evidence each PubMed candidate
    actually provides for a ligand-receptor pair.

    The default path is deterministic and local, so the platform remains
    usable without Gemini, API keys or quota. Set use_llm=True explicitly
    to request the Gemini classifier.
    """

    ligand = _clean_text(
        ligand
    ).upper()

    receptor = _clean_text(
        receptor
    ).upper()

    if not ligand:
        raise ValueError(
            "A ligand identifier is required."
        )

    if not receptor:
        raise ValueError(
            "A receptor identifier is required."
        )

    resolved_model = (
        _clean_text(
            model
        )
        or DEFAULT_MODEL
    )

    if not articles:
        return {
            "ligand":
                ligand,

            "receptor":
                receptor,

            "provider":
                "Local rule-based",

            "model":
                "deterministic-fallback-v1",

            "classifications":
                [],

            "llm_called":
                False,

            "local_rule_based":
                True,
        }

    if not use_llm:
        return {
            "ligand":
                ligand,

            "receptor":
                receptor,

            "provider":
                "Local rule-based",

            "model":
                "deterministic-fallback-v1",

            "classifications":
                _local_classify_articles(
                    ligand=ligand,
                    receptor=receptor,
                    articles=articles,
                ),

            "llm_called":
                False,

            "local_rule_based":
                True,
        }

    cached_by_pmid: dict[
        str,
        dict[str, Any],
    ] = {}

    missing_articles: list[
        dict[str, Any]
    ] = []

    article_hash_by_pmid: dict[
        str,
        str,
    ] = {}

    for article in articles:
        pmid = str(
            article.get(
                "pmid"
            )
            or ""
        ).strip()

        if not pmid:
            continue

        content_hash = _content_hash(
            ligand=ligand,
            receptor=receptor,
            article=article,
        )

        article_hash_by_pmid[
            pmid
        ] = content_hash

        if not force_refresh:
            cached = _load_cached(
                content_hash=content_hash,
                model=resolved_model,
            )

            if cached is not None:
                cached_by_pmid[
                    pmid
                ] = cached

                continue

        missing_articles.append(
            article
        )

    fresh_by_pmid: dict[
        str,
        dict[str, Any],
    ] = {}

    if missing_articles:
        fresh_results = _gemini_request(
            ligand=ligand,
            receptor=receptor,
            articles=missing_articles,
            model=resolved_model,
        )

        requested_pmids = {
            str(
                article.get(
                    "pmid"
                )
                or ""
            ).strip()
            for article
            in missing_articles
        }

        for result in fresh_results:
            pmid = str(
                result.get(
                    "pmid"
                )
                or ""
            ).strip()

            if (
                not pmid
                or pmid not in requested_pmids
            ):
                continue

            result = dict(
                result
            )

            result[
                "cached"
            ] = False

            fresh_by_pmid[
                pmid
            ] = result

            _save_cached(
                content_hash=(
                    article_hash_by_pmid[
                        pmid
                    ]
                ),
                pmid=pmid,
                ligand=ligand,
                receptor=receptor,
                model=resolved_model,
                result=result,
            )

    ordered = []

    for article in articles:
        pmid = str(
            article.get(
                "pmid"
            )
            or ""
        ).strip()

        if not pmid:
            continue

        if pmid in cached_by_pmid:
            ordered.append(
                cached_by_pmid[
                    pmid
                ]
            )

        elif pmid in fresh_by_pmid:
            ordered.append(
                fresh_by_pmid[
                    pmid
                ]
            )

        else:
            ordered.append(
                {
                    "pmid":
                        pmid,

                    "primary_evidence_type":
                        "insufficient_information",

                    "evidence_types":
                        [
                            "insufficient_information"
                        ],

                    "pair_specific_evidence":
                        "unclear",

                    "direct_binding_evidence":
                        "unclear",

                    "study_mode":
                        "unclear",

                    "evidence_scope":
                        "unclear",

                    "what_was_studied":
                        (
                            "No structured classification "
                            "was returned."
                        ),

                    "finding_summary":
                        (
                            "The available publication "
                            "evidence could not be classified."
                        ),

                    "experimental_context":
                        (
                            "Insufficient information."
                        ),

                    "classification_reason":
                        (
                            "No structured Gemini result "
                            "was returned for this PMID."
                        ),

                    "limitations":
                        [
                            "No structured classification "
                            "was returned."
                        ],

                    "extraction_certainty":
                        "low",

                    "pair_specific_pathways":
                        [],

                    "cached":
                        False,
                }
            )

    return {
        "ligand":
            ligand,

        "receptor":
            receptor,

        "provider":
            "Gemini",

        "model":
            resolved_model,

        "classifications":
            ordered,

        "llm_called":
            bool(
                missing_articles
            ),
    }
