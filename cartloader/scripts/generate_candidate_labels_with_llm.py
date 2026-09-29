"""
generate_candidate_labels_with_llm.py -- build candidate label lists (cell types and transcriptional programs)
for one or more tissue contexts, for use with annotate_factors_with_llm.py --candidates.

One LLM call per tissue context (no marker genes are sent), so a list can be reused for every dataset from
that tissue. Lists are cached per model and reasoning effort:
    <cache-dir>/<model>/<effort>/candidates/<organism>/<mode>/<tissue>__candidates.v1.json
A request reuses a list made at the same effort or higher, never lower. The layout matches jevanno's cache,
so lists generated there are reused here and vice versa.

Self-contained (requests; the `anthropic` package only for --api-type claude).

Example:
  python generate_candidate_labels_with_llm.py --organism human --tissue "lung cancer" "lymph node" \
      --api-type claude --out candidates.json
"""

import argparse
import datetime as dt
import inspect
import json
import logging
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from string import Template

import requests

API_KEY_ENV = {"claude": "ANTHROPIC_API_KEY", "openai": "OPENAI_API_KEY", "umgpt": "UMGPT_API_KEY",
               "google": "GEMINI_API_KEY"}
DEFAULT_MODEL = {"claude": "claude-opus-5-5", "openai": "gpt-5.4-mini", "umgpt": "gpt-5-mini",
                 "google": "gemini-3.5-flash"}
DEFAULT_BASE_URL = {"openai": "https://api.openai.com/v1", "umgpt": "https://api.toolkit.umgpt.umich.edu/v1"}
ANTHROPIC_FALLBACK_BETA = "server-side-fallback-2026-07-01"
PROMPT_VERSION = "candidates.v1"
EFFORT_ORDER = ("low", "medium", "high", "xhigh", "max")

# Prompt of jevanno's candidate generator (prompts/candidates.v1.txt), verbatim.
CANDIDATE_PROMPT = """# SYSTEM
You are an expert in single-cell and spatial transcriptomics and in ${organism} tissue biology. You build reference label lists that a classifier will later use to annotate latent factors (topics) learned from spatial transcriptomics of a tissue. Labels must be biologically precise, mutually distinguishable, and useful at the resolution such factors usually resolve.

# USER
Build a Candidate Set for this tissue.

Organism: ${organism}
Tissue context: ${tissue_context}
Annotation mode: ${mode_description}

What to include:
${mode_guidance}
- Aim for ${min_candidates}-${max_candidates} candidates. Cover the tissue thoroughly; prefer distinct, well-supported entities over near-duplicates.
- If the tissue is a tumour, include malignant/tumour cell states typical of this cancer (lineage, differentiation, proliferative, EMT-like, hypoxic, stress, oncofetal / cancer-testis antigen programs) as well as the non-malignant compartments.
- Include a small number of non-biological or low-information labels that latent factors often capture (e.g. low-count background / ambient signal, mitochondrial or ribosomal content), so they are not forced onto real cell types.

For every candidate give:
- label: UpperCamelCase, at most 40 characters, no spaces (e.g. AlveolarType2Cell, CD8TExhausted, G2MCellCycleProgram).
- alias: a terse display form of the label, at most 24 characters (may equal label).
- kind: "cell_type" or "program".
- description: one sentence stating what the entity is AND what distinguishes it from the most similar other candidates in this list.
- canonical_markers: 5-10 characteristic genes, most specific first, as official ${organism} gene symbols (${symbol_case}). Prefer genes that are specific to this entity within this tissue over ubiquitous ones.

Return only the JSON object required by the schema.
"""

MODE_DESCRIPTIONS = {
    "cell_type": "cell types (including cell states/subtypes)",
    "program": "transcriptional programs (pathways, cell-cycle phases, stress/response and differentiation programs)",
    "either": "whichever best describes a factor: a cell type/state or a transcriptional program",
}
_CELL_TYPE_GUIDANCE = (
    "- Cell types and cell states expected in this tissue (epithelial, stromal, vascular, immune, neural, "
    "tissue-specific), at the resolution latent factors usually resolve (major types plus well-known subtypes/states)."
)
_PROGRAM_GUIDANCE = (
    "- Transcriptional programs that commonly appear as latent factors: cell-cycle phases, proliferation, "
    "hypoxia/glycolysis, interferon and inflammatory responses, stress and immediate-early responses, EMT, "
    "ribosome biogenesis/translation, RNA processing, metabolic programs, tissue-specific differentiation or "
    "secretory programs."
)
MODE_GUIDANCE = {"cell_type": _CELL_TYPE_GUIDANCE, "program": _PROGRAM_GUIDANCE,
                 "either": _CELL_TYPE_GUIDANCE + "\n" + _PROGRAM_GUIDANCE}

CANDIDATES_SCHEMA = {
    "type": "object",
    "properties": {"candidates": {"type": "array", "items": {
        "type": "object",
        "properties": {
            "label": {"type": "string"},
            "alias": {"type": "string"},
            "kind": {"type": "string", "enum": ["cell_type", "program"]},
            "description": {"type": "string"},
            "canonical_markers": {"type": "array", "items": {"type": "string"}},
        },
        "required": ["label", "alias", "kind", "description", "canonical_markers"],
        "additionalProperties": False,
    }}},
    "required": ["candidates"],
    "additionalProperties": False,
}

log = logging.getLogger("generate_candidate_labels_with_llm")


def symbol_case(organism):
    o = organism.lower()
    if o in ("human", "homo sapiens", "hs"):
        return "HGNC, upper case, e.g. EPCAM"
    if o in ("mouse", "mus musculus", "mm", "rat"):
        return "MGI, capitalised, e.g. Epcam"
    return "official symbols in the organism's convention"


def normalize_tissue(t):
    return re.sub(r"\s+", " ", t.strip().lower())


def normalize_alias(raw):
    if not raw:
        return "Unknown"
    s = raw.strip().strip("`\"'").splitlines()[0].strip()
    if " " in s:
        s = "".join(p[:1].upper() + p[1:] for p in re.split(r"\s+", s) if p)
    s = re.sub(r"[^A-Za-z0-9+\-/]", "", s)
    return s or "Unknown"


def _safe(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "-", name).strip("-") or "unknown"


def cache_key(organism, tissue, mode):
    return {"organism": organism.strip().lower(), "tissue_context": normalize_tissue(tissue),
            "annotation_mode": mode, "prompt_version": PROMPT_VERSION}


def cache_path(cache_dir, key, model, effort):
    slug = re.sub(r"[^a-z0-9]+", "-", key["tissue_context"]).strip("-")
    return os.path.join(cache_dir, _safe(model), _safe(effort), "candidates", key["organism"],
                        key["annotation_mode"], f"{slug}__{key['prompt_version']}.json")


def find_cached(cache_dir, key, model, effort):
    efforts = list(EFFORT_ORDER[EFFORT_ORDER.index(effort):]) if effort in EFFORT_ORDER else [effort]
    for e in efforts:
        path = cache_path(cache_dir, key, model, e)
        if os.path.exists(path):
            return path
    return None


def render(template, **values):
    system, user = template.split("# USER", 1)
    system = system.replace("# SYSTEM", "", 1).strip()
    return Template(system).substitute(values), Template(user.strip()).substitute(values)


def complete_json(args, system, prompt, schema):
    if args.api_type == "claude":
        import anthropic  # only needed for --api-type claude

        client = anthropic.Anthropic()
        kwargs = {"model": args.model_name, "max_tokens": args.max_output_tokens, "system": system,
                  "messages": [{"role": "user", "content": prompt}],
                  "output_config": {"effort": args.effort, "format": {"type": "json_schema", "schema": schema}}}
        stream_cm = (client.messages.stream(**kwargs) if args.no_fallbacks else client.beta.messages.stream(
            betas=[ANTHROPIC_FALLBACK_BETA], extra_body={"fallbacks": "default"}, **kwargs))
        with stream_cm as stream:
            message = stream.get_final_message()
        if message.stop_reason in ("refusal", "max_tokens"):
            raise RuntimeError(f"Claude stopped with {message.stop_reason}")
        text = next(b.text for b in message.content if b.type == "text")
        return json.loads(text), message.model, {"input_tokens": message.usage.input_tokens,
                                                 "output_tokens": message.usage.output_tokens}
    key = os.environ.get(API_KEY_ENV[args.api_type], "")
    if args.api_type in ("openai", "umgpt"):
        url = (args.api_base_url or DEFAULT_BASE_URL[args.api_type]).rstrip("/")
        url = url if url.endswith("/responses") else url + "/responses"
        payload = {"model": args.model_name, "instructions": system, "input": prompt,
                   "max_output_tokens": args.max_output_tokens,
                   "text": {"format": {"type": "json_schema", "name": "candidates", "schema": schema, "strict": True}}}
        if args.effort in ("low", "medium", "high"):
            payload["reasoning"] = {"effort": args.effort}
        r = requests.post(url, headers={"Authorization": f"Bearer {key}"}, json=payload, timeout=args.request_timeout)
        if r.status_code == 400 and "reasoning" in r.text:
            payload.pop("reasoning", None)
            r = requests.post(url, headers={"Authorization": f"Bearer {key}"}, json=payload,
                              timeout=args.request_timeout)
        r.raise_for_status()
        data = r.json()
        text = data.get("output_text") or "".join(c.get("text", "") for item in data.get("output", []) or []
                                                   for c in item.get("content", []) or [] if c.get("type") == "output_text")
        return json.loads(text), data.get("model", args.model_name), data.get("usage", {})
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{args.model_name}:generateContent"
    payload = {"systemInstruction": {"parts": [{"text": system}]},
               "contents": [{"role": "user", "parts": [{"text": prompt}]}],
               "generationConfig": {"responseMimeType": "application/json", "responseJsonSchema": schema}}
    r = requests.post(url, headers={"x-goog-api-key": key}, json=payload, timeout=args.request_timeout)
    r.raise_for_status()
    data = r.json()
    return (json.loads("".join(p.get("text", "") for p in data["candidates"][0]["content"]["parts"])),
            args.model_name, data.get("usageMetadata", {}))


def generate(args, tissue):
    key = cache_key(args.organism, tissue, args.mode)
    found = find_cached(args.cache_dir, key, args.model_name, args.effort)
    if found:
        log.info("%s: cached list %s", tissue, found)
        with open(found, encoding="utf-8") as fh:
            return json.load(fh)
    path = cache_path(args.cache_dir, key, args.model_name, args.effort)
    if args.no_generate:
        raise FileNotFoundError(f"no cached list for {tissue!r} at effort {args.effort} or higher ({path})")
    system, user = render(CANDIDATE_PROMPT, organism=args.organism, tissue_context=tissue,
                          min_candidates=args.min_candidates, max_candidates=args.max_candidates,
                          mode_description=MODE_DESCRIPTIONS[args.mode], mode_guidance=MODE_GUIDANCE[args.mode],
                          symbol_case=symbol_case(args.organism))
    log.info("%s: generating with %s (%s, effort %s)", tissue, args.api_type, args.model_name, args.effort)
    data, model, usage = complete_json(args, system, user, CANDIDATES_SCHEMA)
    seen, cands = set(), []
    for it in data["candidates"]:
        label = normalize_alias(it["label"])
        if label in seen or label == "Unknown":
            continue
        seen.add(label)
        cands.append({"label": label, "alias": normalize_alias(it.get("alias") or label),
                      "description": it["description"].strip(),
                      "canonical_markers": [g.strip() for g in it["canonical_markers"] if g.strip()],
                      "kind": it["kind"] if args.mode == "either" else args.mode,
                      "tissues": [normalize_tissue(tissue)], "source": "tissue"})
    out = {"schema_version": 1, "key": key, "candidates": cands,
           "generated_by": {"provider": args.api_type, "model": model, "effort": args.effort, "usage": usage,
                            "created_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),  # noqa: UP017 (py<3.11)
                            "tissue_context_as_given": tissue}}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    log.info("%s: %d candidates (%s) -> %s", tissue, len(cands), usage, path)
    return out


def generate_candidate_labels_with_llm(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Generate candidate label lists per tissue context (cached per model and effort).")
    io = parser.add_argument_group("Input/Output Parameters")
    io.add_argument("--organism", required=True, help="e.g. human, mouse")
    io.add_argument("--tissue", required=True, nargs="+", help="One or more tissue contexts")
    io.add_argument("--out", required=True, help="Output JSON with the candidates of all tissues (for --candidates)")
    io.add_argument("--api-type", required=True, choices=list(API_KEY_ENV), help="LLM provider")
    aux = parser.add_argument_group("Auxiliary Parameters")
    aux.add_argument("--mode", default="either", choices=list(MODE_DESCRIPTIONS))
    aux.add_argument("--cache-dir", default="jevanno_cache", help="Cache root (default: ./jevanno_cache)")
    aux.add_argument("--no-generate", action="store_true", help="Use cached lists only; fail if one is missing")
    aux.add_argument("--min-candidates", type=int, default=80)
    aux.add_argument("--max-candidates", type=int, default=200)
    aux.add_argument("--model-name", help=f"LLM model (defaults: {DEFAULT_MODEL})")
    aux.add_argument("--effort", default="low", help="Reasoning effort (default: low)")
    aux.add_argument("--max-output-tokens", type=int, default=64000)
    aux.add_argument("--no-fallbacks", action="store_true", help="claude: disable server-side refusal fallbacks")
    aux.add_argument("--api-base-url", help="Base URL for openai/umgpt")
    aux.add_argument("--request-timeout", type=int, default=900)
    aux.add_argument("--threads", type=int, default=8)
    aux.add_argument("--verbose", action="store_true")
    args = parser.parse_args(_args)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="[%(asctime)s - %(levelname)s - %(message)s]", datefmt="%Y-%m-%d %H:%M:%S")
    args.model_name = args.model_name or DEFAULT_MODEL[args.api_type]

    tissues = list(dict.fromkeys(args.tissue))
    with ThreadPoolExecutor(max_workers=max(1, min(args.threads, len(tissues)))) as ex:
        sets = list(ex.map(lambda t: generate(args, t), tissues))
    merged = {}
    for s in sets:
        for c in s["candidates"]:
            if c["label"] in merged:
                m = merged[c["label"]]
                m["tissues"] = list(dict.fromkeys(m["tissues"] + c["tissues"]))
                m["canonical_markers"] = list(dict.fromkeys(m["canonical_markers"] + c["canonical_markers"]))[:12]
            else:
                merged[c["label"]] = {**c, "tissues": list(c["tissues"])}
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump({"schema_version": 1, "key": {"union_of": [s["key"] for s in sets]},
                   "generated_by": {"union_of": [s["generated_by"] for s in sets]},
                   "candidates": list(merged.values())}, fh, indent=2, ensure_ascii=False)
    log.info("Wrote %d candidate labels from %d tissue context(s) to %s", len(merged), len(tissues), args.out)


if __name__ == "__main__":
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], script_name)
    func(sys.argv[1:])
