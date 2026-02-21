import argparse, os, sys, gzip, json, re, logging, inspect, time, base64
from venv import logger
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

# Optional overrides
OPENAI_MODEL_ENV = "OPENAI_MODEL"
GOOGLE_MODEL_ENV = "GOOGLE_MODEL"
ANTHROPIC_MODEL_ENV = "ANTHROPIC_MODEL"

# -----------------------------
# Defaults
# -----------------------------
DEFAULT_OPENAI_MODEL = "gpt-5-mini"
DEFAULT_GOOGLE_MODEL = "gemini-3-flash-preview"
DEFAULT_ANTHROPIC_MODEL = "claude-opus-4-6"

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


def _make_prompt(tissue: str, organism: str, genes: List[str], has_thumbnail: bool) -> str:
    """
    Ask for a single most likely cell type. Force JSON output with alias only.
    If has_thumbnail=True, the model will also receive a tissue thumbnail image
    representing spatial gene expression across ALL factors.
    """
    gene_list = ", ".join(genes)

    thumbnail_block = ""
    if has_thumbnail:
        thumbnail_block = (
            "\nAdditional input:\n"
            "- A thumbnail image is attached. It is a tissue overview showing spatial gene-expression patterns "
            "for all factors (composite/overview). Use it ONLY as supporting evidence to disambiguate the cell type "
            "suggested by the marker genes.\n"
            "- Prioritize marker genes first; use spatial pattern second.\n"
        )

    return (
        "You are annotating latent factors from bulk differential expression.\n"
        "Task: Identify the single most likely cell type represented by these ordered, top marker genes.\n\n"
        f"Organism: {organism}\n"
        f"Tissue: {tissue}\n"
        f"Top marker genes (comma-separated): {gene_list}\n"
        f"{thumbnail_block}\n"
        "Return ONLY a JSON object with this exact schema:\n"
        "{\"alias\": \"UpperCamelCaseCellType\"}\n\n"
        "Rules:\n"
        "- alias must be terse and singular.\n"
        "- Use informative shorthand when appropriate (e.g., CD4+T, CD8+T, NKCell, BCell, PlasmaCell).\n"
        "- Avoid long phrases, parentheses, or multi-sentence explanations.\n"
        "- Do not include any extra keys besides 'alias'.\n"
    )



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
  
def _load_thumbnail_b64(path: Optional[str]) -> Optional[Tuple[str, str]]:
    """
    Returns (mime_type, base64_string) for an image file, or None if no path.
    Currently expects PNG/JPG/JPEG by extension.
    """
    if not path:
        return None

    if not os.path.exists(path):
        raise FileNotFoundError(f"--thumbnail not found: {path}")

    ext = os.path.splitext(path)[1].lower()
    if ext == ".png":
        mime = "image/png"
    elif ext in (".jpg", ".jpeg"):
        mime = "image/jpeg"
    else:
        raise ValueError(f"--thumbnail must be .png/.jpg/.jpeg (got {ext})")

    with open(path, "rb") as f:
        b = f.read()

    b64 = base64.b64encode(b).decode("utf-8")
    return mime, b64

# -----------------------------
# API Calls
# -----------------------------
def call_openai(prompt: str, thumbnail: Optional[Tuple[str, str]], model_name: str, request_timeout: int, max_retries: int) -> str:
    """
    OpenAI Responses API via REST.
    Env: OPENAI_API_KEY, optional OPENAI_MODEL
    """
    api_key = os.environ.get(OPENAI_API_KEY_ENV, "")
    model = os.environ.get(OPENAI_MODEL_ENV, DEFAULT_OPENAI_MODEL) if model_name is None else model_name

    # Updated implementation for gpt-5 series
    url = "https://api.openai.com/v1/responses"  # New endpoint
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
    }
    content = [{"type": "input_text", "text": prompt}]
    if thumbnail is not None:
        mime, b64 = thumbnail
        content.append({
            "type": "input_image",
            "image_url": f"data:{mime};base64,{b64}",
        })

    payload = {
        "model": model,
        "input": [
            {
                "role": "user",
                "content": content
            }
        ],
        "reasoning": { "effort": "medium" }
    }

    for attempt in range(max_retries):
        try:
            r = requests.post(url, headers=headers, json=payload, timeout=request_timeout)
            r.raise_for_status()
            data = r.json()
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


def call_google(prompt: str, thumbnail: Optional[Tuple[str, str]], model_name: str, request_timeout: int, max_retries: int) -> str:
    """
    Google Gemini 3 SDK call.
    Env: GEMINI_API_KEY, optional GOOGLE_MODEL
    """
    # The SDK automatically looks for the GEMINI_API_KEY environment variable
    client = genai.Client()
    model_id = os.environ.get(GOOGLE_MODEL_ENV, DEFAULT_GOOGLE_MODEL) if model_name is None else model_name

    config = types.GenerateContentConfig(
        thinking_config=types.ThinkingConfig(thinking_level=types.ThinkingLevel.HIGH),
    )

    # Build multimodal contents
    if thumbnail is None:
        contents = prompt
    else:
        mime, b64 = thumbnail
        img_bytes = base64.b64decode(b64)
        contents = [
            prompt,
            types.Part.from_bytes(data=img_bytes, mime_type=mime),
        ]

    for attempt in range(max_retries):
        try:
            response = client.models.generate_content(
                model=model_id,
                contents=contents,
                config=config
            )
            return response.text.strip()
        except Exception as e:
            if attempt == max_retries - 1:
                raise e
            time.sleep(2 ** attempt + 5) ## wait for 7, 9, 13, ... seconds
    return ""


def call_claude(prompt: str, thumbnail: Optional[Tuple[str, str]], model_name: str, request_timeout: int, max_retries: int) -> str:
    """
    Anthropic Messages API via REST.
    Env: ANTHROPIC_API_KEY, optional ANTHROPIC_MODEL
    """
    api_key = os.environ.get(ANTHROPIC_API_KEY_ENV, "")
    model = os.environ.get(ANTHROPIC_MODEL_ENV, DEFAULT_ANTHROPIC_MODEL) if model_name is None else model_name

    url = "https://api.anthropic.com/v1/messages"
    headers = {
        "x-api-key": api_key,
        "anthropic-version": "2023-06-01",
        "content-type": "application/json",
    }
    content_blocks = []
    if thumbnail is not None:
        mime, b64 = thumbnail
        content_blocks.append({
            "type": "image",
            "source": {
                "type": "base64",
                "media_type": mime,
                "data": b64
            }
        })
    content_blocks.append({"type": "text", "text": prompt})

    payload = {
        "model": model,
        "max_tokens": 256,
        "temperature": 0.2,
        "messages": [{"role": "user", "content": content_blocks}],
    }

    for attempt in range(max_retries):
        try:
            r = requests.post(url, headers=headers, json=payload, timeout=request_timeout)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            if attempt == max_retries - 1:
                raise e
            time.sleep(2 ** attempt + 5) ## wait for 7, 9, 13, ... seconds

    # content is a list of blocks; extract text blocks
    try:
        blocks = data.get("content", [])
        texts = [b.get("text", "") for b in blocks if b.get("type") == "text"]
        return "\n".join(texts).strip()
    except Exception:
        return ""


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
    tissue: str,
    organism: str,
    api_type: str,
    model_name: str,
    request_timeout: int,
    max_retries: int,
    thumbnail: Optional[Tuple[str, str]],
    logger: logging.Logger
) -> Dict[str, List[Tuple[int, str]]]:
    """
    Returns dict keyed by engine name -> list of (factor_index, alias)
    """
    results: List[Tuple[int, str]] = []

    ## keep track of duplicates
    alias2cnts = {}

    for idx in sorted(factor2genes.keys()):
        logger.info(f"Annotating factor {idx} with {api_type} API...")
        genes = factor2genes[idx]
        prompt = _make_prompt(tissue=tissue, organism=organism, genes=genes, has_thumbnail=(thumbnail is not None))

        if api_type == "openai":
            text = call_openai(prompt, thumbnail, model_name, request_timeout, max_retries)
            alias = _extract_json_alias(text)
            results.append([idx, alias])
        elif api_type == "google":
            text = call_google(prompt, thumbnail, model_name, request_timeout, max_retries)
            alias = _extract_json_alias(text)
            results.append([idx, alias])
        elif api_type == "claude":
            text = call_claude(prompt, thumbnail, model_name, request_timeout, max_retries)
            alias = _extract_json_alias(text)
            results.append([idx, alias])
        else:
            raise ValueError(f"Unknown API type: {api_type}")
        alias2cnts[alias] = alias2cnts.get(alias, 0) + 1

    ## rename duplicate factors
    logger.info(f"Resolving duplicate aliases...")
    alias2iter = {}
    for i in range(len(results)):
        idx, alias = results[i]
        if alias2cnts[alias] > 1:
            alias2iter[alias] = alias2iter.get(alias,0) + 1
            results[i][1] = f"{alias}_{alias2iter[alias]}"

    return results

def write_output(out: str, rows: List[Tuple[int, str]]) -> str:
    df_out = pd.DataFrame(rows, columns=["index", "alias"])
    df_out = df_out.sort_values("index")
    df_out.to_csv(out, sep="\t", index=False)
    return out

def annotate_bulk_de_with_ai(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="""
                                     Annotate Bulk DE test results with Generative AI APIs.
                                     """)
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--de', required= True, type=str, help='TSV file containing bulk DE test results.')
    inout_params.add_argument('--out', required=True, type=str, help='Output file name')
    inout_params.add_argument('--tissue', type=str, required=True, help='Tissue name used in the prompt')
    inout_params.add_argument('--organism', default="human", help='Organism name for gene interpretation (e.g., human, mouse)')
    inout_params.add_argument('--api-type', type=str, required=True, choices=['openai', 'google', 'claude'], help='API for generative AI model')
    inout_params.add_argument('--thumbnail', type=str, default=None, help='Optional PNG/JPG thumbnail file location of tissue overview with gene expression of all factors; will be included in the prompt as an image')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Other parameters")
    aux_params.add_argument('--primary-rank', type=str, default="Chi2", help='Primary ranking column name in the input TSV (e.g., Chi2)')
    aux_params.add_argument('--secondary-rank', type=str, default="FoldChange", help='Secondary ranking column name in the input TSV (e.g., FoldChange)')
    aux_params.add_argument('--top-n', type=int, default=10, help='Number of top genes to use for annotation (default: 10)')
    aux_params.add_argument('--model-name', type=str, help='Model name for generative AI API. Default will be used otherwise')
    aux_params.add_argument('--request-timeout', type=int, default=60, help='Request timeout (in seconds) for generative AI API (default: 60)')
    aux_params.add_argument('--max-retries', type=int, default=3, help='Maximum number of retries for failed requests (default: 3)')

    args = parser.parse_args(_args)
    
    thumbnail = _load_thumbnail_b64(args.thumbnail)

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

    logger.info(f"Annotating {len(factor2genes)} factors with {args.api_type} API...")

    # Annotate
    results = annotate_factors(
        factor2genes=factor2genes,
        tissue=args.tissue,
        organism=args.organism,
        api_type=args.api_type,
        model_name=args.model_name,
        request_timeout=args.request_timeout,
        max_retries=args.max_retries,
        thumbnail=thumbnail,
        logger=logger
    )

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
