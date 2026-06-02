import argparse, os, sys, gzip, json, re, logging, inspect, time
from venv import logger
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
import pandas as pd
from typing import Dict, List, Optional, Tuple
from google import genai
from google.genai import types


# -----------------------------
# Env var names (user sets these in terminal)
# -----------------------------
OPENAI_API_KEY_ENV = "OPENAI_API_KEY"
GEMINI_API_KEY_ENV = "GEMINI_API_KEY"
ANTHROPIC_API_KEY_ENV = "ANTHROPIC_API_KEY"
UMGPT_API_KEY_ENV = "UMGPT_API_KEY"

# Optional overrides
OPENAI_MODEL_ENV = "OPENAI_MODEL"
GOOGLE_MODEL_ENV = "GOOGLE_MODEL"
ANTHROPIC_MODEL_ENV = "ANTHROPIC_MODEL"
UMGPT_MODEL_ENV = "UMGPT_MODEL"

# -----------------------------
# Defaults
# -----------------------------
DEFAULT_OPENAI_MODEL = "gpt-5.4-mini"
DEFAULT_GOOGLE_MODEL = "gemini-3.5-flash"
DEFAULT_ANTHROPIC_MODEL = "claude-opus-4-8"
DEFAULT_UMGPT_MODEL = "gpt-5-mini"

# U-M GPT Toolkit gateway (OpenAI-compatible). The working API root ends in /v1;
# the /responses endpoint serves all model families (GPT, Claude, Gemini) with
# Bearer auth and a uniform request/response shape.
DEFAULT_UMGPT_BASE_URL = "https://api.toolkit.umgpt.umich.edu/v1"

# Substrings identifying reasoning-capable models on the UMGPT gateway. Only
# these accept a 'reasoning.effort' parameter; others (e.g. gpt-4o, gpt-4.1,
# Llama) reject it with HTTP 400.
UMGPT_REASONING_HINTS = ("gpt-5", "o1", "o3", "claude-sonnet", "claude-opus", "gemini")

# Reasoning models (especially the GPT-5 series) spend tokens on hidden
# reasoning BEFORE emitting visible output. For this task the answer is a tiny
# JSON alias, so we use LOW effort and a large output budget so reasoning can
# finish AND still leave room for the JSON. If a response still comes back
# truncated (incomplete / max_output_tokens), the budget is escalated and the
# request retried once more.
UMGPT_REASONING_EFFORT = "low"
UMGPT_MAX_OUTPUT_TOKENS = 2048
UMGPT_MAX_OUTPUT_TOKENS_ESCALATED = 8192

# TOP_N_GENES = 10
# REQUEST_TIMEOUT_S = 60

# -----------------------------
# Helpers: formatting + parsing
# -----------------------------
def _normalize_alias(raw: str) -> str:
    """
    Normalize an alias into a compact UpperCamelCase-ish token while preserving
    common immunology symbols (e.g., CD4+T, NKCell, BCell, Macrophage).

    If the model already returns something good, we keep it mostly intact
    (just stripping whitespace and disallowed punctuation).
    """
    if raw is None:
        return "Unknown"

    s = raw.strip()

    # If JSON accidentally returned as text, try to extract alias from it.
    # (Caller also parses JSON, but this is a last-ditch fallback.)
    if s.startswith("{") and s.endswith("}"):
        try:
            obj = json.loads(s)
            if isinstance(obj, dict) and "alias" in obj:
                s = str(obj["alias"]).strip()
        except Exception:
            pass

    # Remove surrounding quotes/backticks
    s = s.strip("`\"'")

    # Truncate at first newline or sentence if model rambles
    s = s.splitlines()[0].strip()
    s = re.split(r"[.;]", s, maxsplit=1)[0].strip()

    # If there are spaces, convert to UpperCamelCase but keep + - / digits.
    if " " in s:
        parts = re.split(r"\s+", s)
        s = "".join(p[:1].upper() + p[1:] for p in parts if p)

    # Strip characters that tend to break "alias" formatting, but keep +-/0-9
    s = re.sub(r"[^A-Za-z0-9+\-/]", "", s)

    # Empty fallback
    return s if s else "Unknown"


# def _make_prompt(tissue: str, organism: str, genes: List[str], infer_type: str) -> str:
#     """
#     Ask for a single most likely {infer_type}. Force JSON output with alias only.
#     """
#     gene_list = ", ".join(genes)
#     return (
#         "You are annotating latent factors from bulk differential expression.\n"
#         f"Task: Identify the single most likely {infer_type} represented by these ordered, top marker genes.\n\n"
#         f"Organism: {organism}\n"
#         f"Tissue: {tissue}\n"
#         f"Top marker genes (comma-separated): {gene_list}\n\n"
#         "Return ONLY a JSON object with this exact schema:\n"
#         "{\"alias\": \"UpperCamelCaseInferredCellTypeOrProgramName\"}\n\n"
#         "Rules:\n"
#         "- alias must be terse and singular.\n"
#         "- Use informative shorthand when appropriate (e.g., CD4+T, CD8+T, NKCell, BCell, G2MPhaseCellCycle, CapillaryCaveolarTransportProgram, LipidTransportingCapillaryEndothelium, ImprintedGrowthMetabolicProgram).\n"
#         "- Avoid long phrases, parentheses, or multi-sentence explanations.\n"
#     )

def _make_prompt(template_str:str, template_keyvals: Dict[str, str], genes: List[str]) -> str:
    """
    Create a prompt using a template string and key-value pairs.
    """
    gene_list = ", ".join(genes)
    keyvals = template_keyvals.copy()
    keyvals["top_marker_gene_list"] = gene_list
    prompt = template_str
    for k, v in keyvals.items():
        placeholder = "{" + k + "}"
        ## make sure the placeholder exists in the template
        if placeholder not in prompt:
            raise ValueError(f"Argument {placeholder} not found in template. Available args: {re.findall(r'{(.*?)}', template_str)}")
        prompt = prompt.replace(placeholder, v)
    ## make sure that all placeholders were replaced
    if re.search(r"\{\S*?\}", prompt):
        unreplaced = re.findall(r"\{\S*?\}", prompt)
        raise ValueError(f"Unreplaced arguments {unreplaced} remain in the prompt. Please provide values for these arguments using --template-args or remove them from the template.")
    return prompt

def _extract_json_alias(text: str) -> str:
    """
    Parse the model's response for a JSON object and extract alias.
    Accepts either pure JSON or JSON embedded in text.
    """
    if not text:
        return "Unknown"

    t = text.strip()

    # Fast path: pure JSON
    if t.startswith("{") and t.endswith("}"):
        try:
            obj = json.loads(t)
            if isinstance(obj, dict) and "alias" in obj:
                return _normalize_alias(str(obj["alias"]))
        except Exception:
            pass

    # Embedded JSON: find first {...}
    m = re.search(r"\{.*\}", t, flags=re.DOTALL)
    if m:
        try:
            obj = json.loads(m.group(0))
            if isinstance(obj, dict) and "alias" in obj:
                return _normalize_alias(str(obj["alias"]))
        except Exception:
            pass

    # Fallback: treat whole text as alias
    return _normalize_alias(t)

# -----------------------------
# API Calls
# -----------------------------
def call_openai(prompt: str, model_name: str, request_timeout: int, max_retries: int, api_base_url: Optional[str]) -> str:
    """
    OpenAI Responses API via REST.
    Env: OPENAI_API_KEY, optional OPENAI_MODEL
    """
    api_key = os.environ.get(OPENAI_API_KEY_ENV, "")
    model = os.environ.get(OPENAI_MODEL_ENV, DEFAULT_OPENAI_MODEL) if model_name is None else model_name

    # Updated implementation for gpt-5 series
    url = api_base_url or "https://api.openai.com/v1/responses"  # New endpoint
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
    }
    payload = {
        "model": model,
        "input": [  # 'messages' is now 'input'
            {
                "role": "user", 
                "content": [{"type": "input_text", "text": prompt}] # Itemized content
            }
        ],
        # REMOVED: temperature, top_p, etc.
        "reasoning": { "effort": "medium" } # Optional: Controls thinking depth
    }

    for attempt in range(max_retries):
        try:
            r = requests.post(url, headers=headers, json=payload, timeout=request_timeout)
            r.raise_for_status()
            data = r.json()
            break
        except Exception as e:
            if attempt == max_retries - 1:
                raise e
            time.sleep(2 ** attempt + 5) ## wait for 7, 9, 13, ... seconds

    # Try common fields for Responses API
    # - output_text is present in many SDKs; in raw JSON, text is often in output[].content[].text
    if isinstance(data, dict) and "output_text" in data and data["output_text"]:
        return str(data["output_text"])

    # Fallback parse
    try:
        out = data.get("output", [])
        texts = []
        for item in out:
            for c in item.get("content", []):
                if c.get("type") == "output_text" and "text" in c:
                    texts.append(c["text"])
        return "\n".join(texts).strip()
    except Exception:
        return ""


def call_google(prompt: str, model_name: str, request_timeout: int, max_retries: int, api_base_url: Optional[str]) -> str:
    """
    Google Gemini 3 SDK call.
    Env: GEMINI_API_KEY, optional GOOGLE_MODEL
    """
    # The SDK automatically looks for the GEMINI_API_KEY environment variable
    client = genai.Client()

    if api_base_url:
        client._client._base_url = api_base_url  # Override base URL for custom/proxy endpoints
    
    model_id = os.environ.get(GOOGLE_MODEL_ENV, DEFAULT_GOOGLE_MODEL) if model_name is None else model_name
    
    # Configure Gemini 3 reasoning behavior
    # Options for thinking_level: LOW, MEDIUM, HIGH, or MINIMAL (Flash only)
    config = types.GenerateContentConfig(
        thinking_config=types.ThinkingConfig(
            thinking_level=types.ThinkingLevel.HIGH 
        ),
        # temperature=1.0 # Recommended to leave at default for Gemini 3
    )

    for attempt in range(max_retries):
        try:
            response = client.models.generate_content(
                model=model_id,
                contents=prompt,
                config=config
            )
            return response.text.strip()
        except Exception as e:
            if attempt == max_retries - 1:
                raise e
            time.sleep(2 ** attempt + 5) ## wait for 7, 9, 13, ... seconds
    return ""

def call_claude(prompt: str, model_name: str, request_timeout: int, max_retries: int, api_base_url: Optional[str]) -> str:
    """
    Anthropic Messages API via REST.
    Env: ANTHROPIC_API_KEY, optional ANTHROPIC_MODEL

    When api_base_url is provided (e.g. the U-M GPT Toolkit gateway), the API key
    is sent via the 'Authorization: Bearer' header. Otherwise the standard
    Anthropic 'x-api-key' header is used against api.anthropic.com.
    """
    api_key = os.environ.get(ANTHROPIC_API_KEY_ENV, "")
    model = os.environ.get(ANTHROPIC_MODEL_ENV, DEFAULT_ANTHROPIC_MODEL) if model_name is None else model_name

    if api_base_url:
        url = api_base_url
        headers = {
            "Authorization": f"Bearer {api_key}",
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        }
    else:
        url = "https://api.anthropic.com/v1/messages"
        headers = {
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        }

    payload = {
        "model": model,
        "max_tokens": 256,
        "messages": [{"role": "user", "content": prompt}],
    }

    data = None
    for attempt in range(max_retries):
        try:
            r = requests.post(url, headers=headers, json=payload, timeout=request_timeout)
            r.raise_for_status()
            data = r.json()
            break
        except Exception as e:
            if attempt == max_retries - 1:
                raise e
            time.sleep(2 ** attempt + 5)  ## wait for 7, 9, 13, ... seconds

    # content is a list of blocks; extract text blocks
    try:
        blocks = data.get("content", [])
        texts = [b.get("text", "") for b in blocks if b.get("type") == "text"]
        return "\n".join(texts).strip()
    except Exception:
        return ""


def _umgpt_supports_reasoning(model: str) -> bool:
    low = model.lower()
    return any(h in low for h in UMGPT_REASONING_HINTS)


def _umgpt_extract_text(data: dict) -> str:
    """Extract assistant text from a /responses payload, tolerating shapes."""
    if not isinstance(data, dict):
        return ""
    # Convenience field present in some responses.
    if data.get("output_text"):
        return str(data["output_text"]).strip()
    # Walk output[].content[] for output_text blocks.
    try:
        texts = []
        for item in data.get("output", []) or []:
            for c in item.get("content", []) or []:
                if c.get("type") == "output_text" and "text" in c:
                    texts.append(c["text"])
        return "\n".join(texts).strip()
    except Exception:
        return ""


def _umgpt_is_truncated(data: dict) -> bool:
    """True when the response was cut off before finishing (reasoning ate the
    whole output budget), which yields empty/partial text and downstream
    'Unknown'. Detected via status/incomplete_details, or empty output where
    usage shows the cap was reached."""
    if not isinstance(data, dict):
        return False
    if data.get("status") == "incomplete":
        return True
    inc = data.get("incomplete_details")
    if isinstance(inc, dict) and inc.get("reason") == "max_output_tokens":
        return True
    return False


def call_umgpt(prompt: str, model_name: str, request_timeout: int, max_retries: int, api_base_url: Optional[str]) -> str:
    """
    U-M GPT Toolkit gateway via REST (OpenAI-compatible /responses endpoint).
    Env: UMGPT_API_KEY, optional UMGPT_MODEL

    The Toolkit is a single OpenAI-compatible gateway that serves all model
    families (GPT, Claude, Gemini) through /responses with Bearer auth. The
    model name (e.g. 'gpt-5-mini', 'claude-opus-4-7', 'gemini-3-flash-preview')
    selects the underlying model. 'reasoning.effort' is sent only for
    reasoning-capable models, since others reject it with HTTP 400.

    GPT-5-series models spend heavily on hidden reasoning before producing
    output, which can exhaust the token budget and return truncated/empty text
    (downstream this becomes 'Unknown'). To avoid that we use low reasoning
    effort with a large output budget, and if a response still comes back
    truncated we escalate the budget and retry.
    """
    api_key = os.environ.get(UMGPT_API_KEY_ENV, "")
    model = os.environ.get(UMGPT_MODEL_ENV, DEFAULT_UMGPT_MODEL) if model_name is None else model_name

    # Resolve the /responses endpoint, tolerating various forms of base URL:
    #   .../v1                -> append /responses
    #   .../v1/               -> append responses
    #   .../v1/responses      -> use as-is
    #   https://host          -> assume /v1/responses
    base = (api_base_url or DEFAULT_UMGPT_BASE_URL).strip().rstrip("/")
    if base.endswith("/responses"):
        url = base
    elif base.endswith("/v1"):
        url = base + "/responses"
    else:
        url = base + "/v1/responses"

    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
    }
    supports_reasoning = _umgpt_supports_reasoning(model)

    def _build_payload(max_out: int) -> dict:
        p = {
            "model": model,
            "input": prompt,
            "max_output_tokens": max_out,
        }
        if supports_reasoning:
            p["reasoning"] = {"effort": UMGPT_REASONING_EFFORT}
        return p

    # Two budget tiers: normal, then escalated if the first comes back truncated.
    budgets = [UMGPT_MAX_OUTPUT_TOKENS, UMGPT_MAX_OUTPUT_TOKENS_ESCALATED]

    data = None
    text = ""
    for budget in budgets:
        payload = _build_payload(budget)
        data = None
        for attempt in range(max_retries):
            try:
                r = requests.post(url, headers=headers, json=payload, timeout=request_timeout)
                r.raise_for_status()
                data = r.json()
                break
            except Exception as e:
                if attempt == max_retries - 1:
                    raise e
                time.sleep(2 ** attempt + 5)  ## wait for 7, 9, 13, ... seconds

        text = _umgpt_extract_text(data)
        # If we got usable text and the response wasn't truncated, we're done.
        if text and not _umgpt_is_truncated(data):
            return text
        # Truncated or empty: escalate the budget and try once more.
        if budget != budgets[-1]:
            logging.getLogger(__name__).warning(
                f"UMGPT response truncated or empty for model '{model}' at "
                f"max_output_tokens={budget}; retrying with larger budget."
            )

    # Return whatever text we managed to extract (may be empty -> 'Unknown').
    return text

# -----------------------------
# Core logic
# -----------------------------
def read_bulk_de(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"gene", "factor"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input missing required columns: {sorted(required)}")

    # Allow either "Chi2" or "FoldChange" columns as in your example
    return df


def top_genes_per_factor(
    df: pd.DataFrame,
    primary_rank: str,
    secondary_rank: str,
    top_n: int,
) -> Dict[int, List[str]]:
    # Validate columns exist
    for col in (primary_rank, secondary_rank):
        if col not in df.columns:
            raise ValueError(f"Ranking column '{col}' not found in input. Available: {list(df.columns)}")

    # Ensure factors are sortable ints (your example uses 0..k-1)
    try:
        factors = sorted(df["factor"].astype(int).unique().tolist())
        df2 = df.copy()
        df2["factor"] = df2["factor"].astype(int)
        df = df2
    except Exception:
        factors = sorted(df["factor"].unique().tolist())

    out: Dict[int, List[str]] = {}

    for f in factors:
        sub = df[df["factor"] == f].copy()

        # Get top N by each ranking column
        prim_list = (
            sub.sort_values(primary_rank, ascending=False)["gene"]
            .astype(str)
            .head(top_n)
            .tolist()
        )
        sec_list = (
            sub.sort_values(secondary_rank, ascending=False)["gene"]
            .astype(str)
            .head(top_n)
            .tolist()
        )

        # Interleave: P1, S1, P2, S2, ... skipping duplicates
        seen = set()
        merged: List[str] = []
        for i in range(top_n):
            if i < len(prim_list):
                g = prim_list[i]
                if g not in seen:
                    merged.append(g)
                    seen.add(g)
            if i < len(sec_list):
                g = sec_list[i]
                if g not in seen:
                    merged.append(g)
                    seen.add(g)

        out[int(f)] = merged

    return out


def annotate_factors(
    factor2genes: Dict[int, List[str]],
    template_str: str,
    template_keyvals: Dict[str, str],
    api_type: str,
    model_name: str,
    request_timeout: int,
    max_retries: int,
    threads: int,
    api_base_url: Optional[str],
    logger: logging.Logger
) -> Dict[str, List[Tuple[int, str]]]:
    """
    Returns dict keyed by engine name -> list of (factor_index, alias)
    """
    if api_type == "openai":
        api_fn = call_openai
    elif api_type == "google":
        api_fn = call_google
    elif api_type == "claude":
        api_fn = call_claude
    elif api_type == "umgpt":
        api_fn = call_umgpt
    else:
        raise ValueError(f"Unknown API type: {api_type}")

    def _annotate_one(idx: int) -> Tuple[int, str]:
        logger.info(f"Annotating factor {idx} with {api_type} API...")
        genes = factor2genes[idx]
        prompt = _make_prompt(template_str=template_str, template_keyvals=template_keyvals, genes=genes)
        text = api_fn(prompt, model_name, request_timeout, max_retries, api_base_url)
        alias = _extract_json_alias(text)
        return idx, alias

    sorted_indices = sorted(factor2genes.keys())
    idx2alias: Dict[int, str] = {}

    n_workers = max(1, int(threads))
    logger.info(f"Annotating {len(sorted_indices)} factors with {n_workers} thread(s)...")

    ## fail-fast: cancel pending tasks and don't wait on running ones if any factor errors
    executor = ThreadPoolExecutor(max_workers=n_workers)
    try:
        future2idx = {executor.submit(_annotate_one, idx): idx for idx in sorted_indices}
        for future in as_completed(future2idx):
            idx = future2idx[future]
            try:
                _, alias = future.result()
            except Exception as e:
                logger.error(f"Failed to annotate factor {idx}: {e}")
                raise
            idx2alias[idx] = alias
    finally:
        executor.shutdown(wait=False, cancel_futures=True)

    ## preserve factor-index order so duplicate suffixing is deterministic
    results: List[List] = [[idx, idx2alias[idx]] for idx in sorted_indices]

    ## keep track of duplicates
    alias2cnts: Dict[str, int] = {}
    for _, alias in results:
        alias2cnts[alias] = alias2cnts.get(alias, 0) + 1

    ## rename duplicate factors
    logger.info(f"Resolving duplicate aliases...")
    alias2iter: Dict[str, int] = {}
    for i in range(len(results)):
        idx, alias = results[i]
        if alias2cnts[alias] > 1:
            alias2iter[alias] = alias2iter.get(alias, 0) + 1
            results[i][1] = f"{alias}_{alias2iter[alias]}"

    return results

def write_output(out: str, rows: List[Tuple[int, str]]) -> str:
    df_out = pd.DataFrame(rows, columns=["index", "alias"])
    df_out = df_out.sort_values("index")
    df_out.to_csv(out, sep="\t", index=False)
    return out

def annotate_bulk_de_with_ai(_args):
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="""
                                     Annotate Bulk DE test results with Generative AI APIs.
                                     """)
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--de', required= True, type=str, help='TSV file containing bulk DE test results.')
    inout_params.add_argument('--out', required=True, type=str, help='Output file name')
    inout_params.add_argument('--template', type=str, default=f"{repo_dir}/assets/template.ai_anno.txt", help='Path to prompt template file (default: assets/template.ai_anno.txt)')
    inout_params.add_argument('--template-args', type=str, nargs='+', metavar="KEY=VALUE", help="Template variables as key=value pairs (e.g., tissue='lung')")
    inout_params.add_argument('--template-argfile', type=str, help='Path to a TSV or JSON file containing template variables. TSV should have two collumns: key and value. JSON should be a flat object with string keys and values.')
    # inout_params.add_argument('--tissue', type=str, required=True, help='Tissue name used in the prompt')
    # inout_params.add_argument('--organism', default="human", help='Organism name for gene interpretation (e.g., human, mouse)')
    inout_params.add_argument('--api-type', type=str, required=True, choices=['openai', 'google', 'claude', 'umgpt'], help='API for generative AI model')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Other parameters")
    aux_params.add_argument('--primary-rank', type=str, default="Chi2", help='Primary ranking column name in the input TSV (e.g., Chi2)')
    aux_params.add_argument('--secondary-rank', type=str, default="FoldChange", help='Secondary ranking column name in the input TSV (e.g., FoldChange)')
    aux_params.add_argument('--top-n', type=int, default=10, help='Number of top genes to use for annotation (default: 10)')
    aux_params.add_argument('--model-name', type=str, help='Model name for generative AI API. Default will be used otherwise')
    # aux_params.add_argument('--infer-type', type=str, default="cell type, subcellular transcriptional program, or known biological pathway", help='Type of entity to infer (e.g., cell type, transcriptional program, etc)')
    aux_params.add_argument('--request-timeout', type=int, default=60, help='Request timeout (in seconds) for generative AI API (default: 60)')
    aux_params.add_argument('--max-retries', type=int, default=3, help='Maximum number of retries for failed requests (default: 3)')
    aux_params.add_argument('--threads', type=int, default=1, help='Number of threads to use for parallel API calls (default: 1)')
    aux_params.add_argument('--api-base-url', type=str, help=f'Base URL for the API endpoint if using a custom or proxy service. For --api-type umgpt, defaults to {DEFAULT_UMGPT_BASE_URL}')

    args = parser.parse_args(_args)

    # For the U-M GPT Toolkit, default the base URL when the user didn't supply
    # one, so the gateway endpoint is used out of the box.
    if args.api_type == "umgpt" and not args.api_base_url:
        args.api_base_url = DEFAULT_UMGPT_BASE_URL

    log_format = "[%(asctime)s - %(levelname)s - %(message)s]"
    date_format = "[%Y-%m-%d %H:%M:%S]"  # Clean timestamp without milliseconds

    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        datefmt=date_format
    )

    logger = logging.getLogger(__name__)

    logger.info(f"Reading DE test results from {args.de}")

    # Read + collect top genes
    df = read_bulk_de(args.de)

    logger.info(f"Collecting top {args.top_n} genes per factor, first by {args.primary_rank}, then by {args.secondary_rank}")
    factor2genes = top_genes_per_factor(
        df,
        primary_rank=args.primary_rank,
        secondary_rank=args.secondary_rank,
        top_n=args.top_n
    )

    ## make sure that the template file exists
    if not os.path.isfile(args.template):
        logger.error(f"Template file not found: {args.template}")
        sys.exit(1)
    ## read the template file and store it in a string
    with open(args.template, "r") as f:
        template_str = f.read()

    ## parse the template arguments
    template_keyvals = {}
    if args.template_args is not None and len(args.template_args) > 0:
        if args.template_argfile:
            logger.error("Cannot specify both --template-args and --template-argfile")
            sys.exit(1)
        for kv in args.template_args:
            if "=" not in kv:
                logger.error(f"Invalid template argument: {kv}. Must be in KEY=VALUE format.")
                sys.exit(1)
            key, value = kv.split("=", 1)
            key = key.strip("'\"")
            value = value.strip("'\"") ## allow values to be optionally enclosed in quotes
            template_keyvals[key] = value
    elif args.template_argfile:
        if not os.path.isfile(args.template_argfile):
            logger.error(f"Template argument file not found: {args.template_argfile}")
            sys.exit(1)
        if args.template_argfile.endswith(".json"):
            with open(args.template_argfile, "r") as f:
                try:
                    obj = json.load(f)
                    if not isinstance(obj, dict):
                        logger.error("Template argument JSON file must contain a flat object with string keys and values")
                        sys.exit(1)
                    for k, v in obj.items():
                        template_keyvals[k] = str(v)
                except Exception as e:
                    logger.error(f"Failed to parse template argument JSON file: {e}")
                    sys.exit(1)
        else: ## assume TSV with key in first column and value in second column
            with open(args.template_argfile, "r") as f:
                for line in f:
                    if "\t" not in line:
                        logger.error(f"Invalid line in template argument TSV file (missing tab): {line.strip()}")
                        sys.exit(1)
                    key, value = line.rstrip().split("\t", maxsplit=1)
                    key = key.rstrip(":") ## allow keys to optionally end with a colon, which is common in templates
                    key = key.strip("'\"")
                    value = value.strip("'\"") ## allow values to be optionally enclosed in quotes
                    template_keyvals[key] = value

    logger.info(f"Annotating {len(factor2genes)} factors with {args.api_type} API...")

    # Annotate
    results = annotate_factors(
        factor2genes=factor2genes,
        template_str=template_str,
        template_keyvals=template_keyvals,
        api_type=args.api_type,
        model_name=args.model_name,
        request_timeout=args.request_timeout,
        max_retries=args.max_retries,
        threads=args.threads,
        api_base_url=args.api_base_url,
        logger=logger
    )

    # results = annotate_factors(
    #     factor2genes=factor2genes,
    #     tissue=args.tissue,
    #     organism=args.organism,
    #     api_type=args.api_type,
    #     model_name=args.model_name,
    #     infer_type=args.infer_type,
    #     request_timeout=args.request_timeout,
    #     max_retries=args.max_retries,
    #     threads=args.threads,
    #     api_base_url=args.api_base_url,
    #     logger=logger
    # )

    logger.info(f"Writing results to {args.out}")

    write_output(args.out, results)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])