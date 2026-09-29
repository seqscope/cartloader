"""
annotate_factors_with_llm.py -- label the latent factors of a cartloader run (cell types or transcriptional
programs) with a generative LLM that sees all factors at once, and write an interactive report.

Point it at a factor prefix:
  single sample       .../cartl/<sample>/t18-f192   reads <prefix>-bulk-de.tsv and <prefix>-model.tsv
  multi-sample root   .../cartl/t12-f192            reads <prefix>-de.tsv and <prefix>-model.tsv (joint factors),
                                                    and each sample's <multi_id>-<sample>/t12-f192-pixel-pseudobulk.tsv.gz
<prefix>-rgb.tsv, when present, gives each factor the colour it has in CartoScope.
Samples of a multi-sample run come from multi-catalog.yaml next to the prefix or, without it, from the sample
sub-directories that hold a pseudobulk for the prefix. The tissue of every sample is --tissue, or the `tissue`
column of --sample-sheet (the run_together sample sheet, matched by `id`) for pan-tissue runs.

Evidence per factor:
  * marker genes: union of the top-N genes by chi-square and by fold enrichment, interleaved (as in
    annotate_bulk_de_with_ai.py)
  * abundance (share of counts), DE support (number of enriched genes), per-gene specificity
  * multi-sample: which samples/tissues hold the factor and its top genes within each sample

All factors are labelled together in one deep-reasoning call (split in halves automatically only if the answer
does not fit the output budget), so the model contrasts factors and may answer "Unresolved". For each factor it
returns a best-guess label, kind (cell type / program), compartment, confidence, key genes and rationale, plus
alternative interpretations when the markers support more than one reading (e.g. a cell type and a
transcriptional program), and a summary of the dataset. Candidate labels (--candidates, from
generate_candidate_labels_with_llm.py, or "auto") are optional suggestions, or a closed list with --closed.

Outputs (default stem: <prefix>-alias-ai):
  <stem>.tsv               index<TAB>alias, the format of annotate_bulk_de_with_ai.py
  <stem>.annotations.tsv   one row per factor: status, label, kind, compartment, confidence, alternatives, evidence
  <stem>.html              self-contained report: overview, factor table, interactive dotplots, run details
  <stem>.llm/              request and response of every LLM call (reused on re-runs unless --no-reuse)

Examples:
  cartloader annotate_factors_with_llm --prefix cartl/sample-rep1/t18-f192 --organism mouse \
      --tissue "whole mouse pup" --api-type claude
  cartloader annotate_factors_with_llm --prefix cartl/t12-f192 --organism human --sample-sheet samples.tsv \
      --api-type claude
"""

import argparse
import base64
import csv
import datetime as dt
import gzip
import hashlib
import inspect
import json
import logging
import math
import os
import re
import shlex
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from string import Template
from types import SimpleNamespace

import numpy as np
import pandas as pd
import requests

# -----------------------------
# Defaults
# -----------------------------
API_KEY_ENV = {"claude": "ANTHROPIC_API_KEY", "openai": "OPENAI_API_KEY", "umgpt": "UMGPT_API_KEY",
               "google": "GEMINI_API_KEY"}
DEFAULT_MODEL = {"claude": "claude-opus-5-5", "openai": "gpt-5.4-mini", "umgpt": "gpt-5-mini",
                 "google": "gemini-3.5-flash"}
DEFAULT_BASE_URL = {"openai": "https://api.openai.com/v1", "umgpt": "https://api.toolkit.umgpt.umich.edu/v1"}
DEFAULT_MAX_OUTPUT_TOKENS = {"claude": 128000, "openai": 128000, "umgpt": 64000, "google": 64000}
ANTHROPIC_FALLBACK_BETA = "server-side-fallback-2026-07-01"
PROMPT_VERSION = "joint.v2"

GENE_COLUMN_ALIASES = ("gene", "feature", "gene_id", "geneid")
PSEUDOCOUNT = 0.5
CHI2_PVAL_MAX = 1e-3
UNRESOLVED = "Unresolved"
CONFIDENCE_LEVELS = {"high": 0.9, "medium": 0.6, "low": 0.3}
MODE_DESCRIPTIONS = {
    "cell_type": "cell types (including cell states/subtypes)",
    "program": "transcriptional programs (pathways, cell-cycle phases, stress/response and differentiation programs)",
    "either": "whichever best describes a factor: a cell type/state or a transcriptional program",
}
# Fixed order: it is also the colour order of the report (validated adjacent-pair palette).
COMPARTMENTS = ("epithelial", "immune", "stromal", "vascular", "neural", "program", "other", "malignant")

# File layout under a factor prefix (first existing wins)
DE_SUFFIXES = ("-bulk-de.tsv", "-de.tsv", "-bulk-de.tsv.gz", "-de.tsv.gz")
MODEL_SUFFIXES = ("-model.tsv", "-model.tsv.gz", "-pseudobulk.tsv.gz")
PSEUDOBULK_SUFFIXES = ("-pixel-pseudobulk.tsv.gz", "-pseudobulk.tsv.gz")
RGB_SUFFIX = "-rgb.tsv"
ALIAS_SUFFIX = "-alias-ai.tsv"
SHEET_UNSET = {"", "-", ".", "NA"}
SHEET_ID_COLUMNS = ("id", "sample", "sample_id")
SHEET_TISSUE_COLUMNS = ("tissue", "tissue_context")

# Joint annotator prompt. v1 is jevanno's evaluated prompt (prompts/joint.v1.txt); v2 keeps its evidence and
# rules and adds alternative interpretations, compartment, key genes and a dataset summary.
JOINT_PROMPT = """# SYSTEM
You are an expert in single-cell and spatial transcriptomics and in ${organism} tissue biology. You annotate the latent factors of a factor model (topic model) fitted to spatial transcriptomics. Each factor is a gene program: it may correspond to a cell type or cell state, or to a transcriptional program shared by several cell types. You annotate all factors of a dataset together, using the contrasts between factors to decide what each one is.

# USER
Organism: ${organism}
Tissue context(s) of the dataset: ${tissue_summary}
Annotation mode: ${mode_description}

How to read the evidence
- marker_genes: the factor's most characteristic genes, most characteristic first. They interleave the top genes by chi-square significance and the top genes by fold enrichment against all other factors, so the list mixes abundant markers with highly specific rare ones.
${evidence_notes}
Candidate labels (${candidate_rule}). Format: label [kind]: canonical markers
${candidates_block}

Factors to annotate in this request (JSON, one per line):
${factors_block}
${digest_block}
Instructions
1. Annotate every factor listed under "Factors to annotate" exactly once, using its id.
2. ${label_rule}
3. Compare factors with each other. When two factors share a biological identity, give them labels that state what distinguishes them (state, program, subtype or sample of origin), not identical labels.
4. A factor can be a program rather than a cell type (e.g. cell cycle, interferon response, hypoxia, stress). Choose the kind that best explains its marker genes.
5. If the marker genes do not support any coherent identity, use the label "Unresolved" with kind "unresolved" and low confidence. Do not guess.
6. confidence: "high" when several specific markers agree, "medium" when the call is likely but partly ambiguous, "low" otherwise.
7. alias: UpperCamelCase display name without spaces, at most 40 characters. Spell words out rather than abbreviating them (OvarianCarcinomaSERPINA5+ not OvCaSERPINA5, NaiveCentralMemoryT not NaiveCMT); keep gene symbols as written and only standard abbreviations (NK, DC, CD8, ECM).
8. alternatives: other interpretations the marker genes also support, best first, at most 3, each with its own kind, confidence and the marker genes that support it. List one whenever a second reading is genuinely plausible: the factor read as a transcriptional program instead of a cell type or vice versa (e.g. proliferating T cells vs a G2/M cell-cycle program), two cell types that share the factor, or a cell state vs its parent type. Leave the list empty when the call is unambiguous.
9. compartment of the label: "epithelial", "immune", "stromal" (fibroblasts, smooth muscle, pericytes, adipocytes, mesenchyme), "vascular" (blood and lymphatic endothelium), "neural" (neurons, glia), "malignant" (tumour cells), "program" (a program not tied to one compartment) or "other".
10. key_genes: the 3-6 marker genes that most support the label.
11. rationale: at most 200 characters naming the genes that drove the call.
12. in_candidate_set: true only if label is copied exactly from the candidate labels.
13. dataset_summary: 3-6 sentences on what these factors show about the dataset as a whole: the compartments and major populations represented, notable programs, sample- or tissue-specific factors, and factors that look technical or unresolved.

Return only the JSON object required by the schema.
"""

EVIDENCE_NOTES_SUMMARY = (
    "- abundance: the factor's share of all counts (\"major\" >= 0.5 %, \"minor\", \"trace\" < 0.1 %). "
    "de_support: how many genes are significantly enriched in the factor (\"strong\" > 150, \"moderate\", "
    "\"weak\" <= 30). A trace factor with weak support and marker genes that share no coherent biology is "
    "usually noise or a very rare population: say so with low confidence rather than forcing a label."
)
EVIDENCE_NOTES_MULTI = (
    "- found_in / distribution: which Samples hold the factor. \"specific\" = essentially one Sample, "
    "\"restricted\" = Samples of one or two tissue contexts, \"shared\" = many. A factor specific to one tumour "
    "Sample often represents that tumour's malignant cells or a Sample-specific program (occasionally a technical "
    "artefact); a factor shared across many tissues usually represents immune, stromal or vascular cells or a "
    "general program. Sample names are given with their tissue context."
)
EVIDENCE_NOTES_SAMPLE_MARKERS = (
    "- markers_by_sample: the factor's top genes recomputed within each Sample that holds it. Differences "
    "between Samples show tissue-specific variants of a shared factor; the label must still fit all of them."
)
EVIDENCE_NOTES_SPECIFICITY = (
    "- marker_specificity: for each marker gene, \"exclusive\" (at least half of its expression is in this factor), "
    "\"highest\" (this factor expresses it most), \"top3\" (among its three highest factors) or \"shared\" "
    "(other factors express it more). Exclusive and highest genes define the factor; shared genes are context."
)

log = logging.getLogger("annotate_factors_with_llm")


class OutputTruncated(RuntimeError):
    """The model ran out of output tokens before finishing the JSON answer."""


# -----------------------------
# Input
# -----------------------------
def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def first_existing(prefix, suffixes):
    for s in suffixes:
        if os.path.exists(prefix + s):
            return prefix + s
    return None


def read_gene_by_factor(path):
    """Gene x factor matrix (model matrix or pixel pseudobulk): first column genes, one column per factor."""
    with open_text(path) as fh:
        df = pd.read_csv(fh, sep="\t")
    df = df.set_index(df.columns[0])
    df.index = df.index.astype(str)
    df.columns = [str(c) for c in df.columns]
    return df.astype(float)


def read_rgb(path):
    """Factor colours from a cartloader rgb.tsv (Name, [Color_index,] R, G, B in 0-1 or 0-255) -> {factor: '#rrggbb'}."""
    with open_text(path) as fh:
        df = pd.read_csv(fh, sep="\t", dtype={"Name": str})
    missing = {"Name", "R", "G", "B"} - set(df.columns)
    if missing:
        raise ValueError(f"RGB TSV {path} is missing columns {sorted(missing)}; has {list(df.columns)}")
    rgb = df[["R", "G", "B"]].astype(float)
    if (rgb.to_numpy() > 1).any():
        rgb = rgb / 255.0
    rgb = (rgb.clip(0, 1) * 255).round().astype(int)
    return {str(n): "#%02x%02x%02x" % tuple(v) for n, v in zip(df["Name"], rgb.to_numpy())}


def read_de(path):
    with open_text(path) as fh:
        df = pd.read_csv(fh, sep="\t")
    if "gene" not in df.columns:
        lower2col = {str(c).lower(): c for c in df.columns}
        for alias in GENE_COLUMN_ALIASES:
            if alias in lower2col:
                df = df.rename(columns={lower2col[alias]: "gene"})
                break
    missing = {"gene", "factor"} - set(df.columns)
    if missing:
        raise ValueError(f"DE TSV {path} is missing columns {sorted(missing)}; has {list(df.columns)}")
    df["gene"] = df["gene"].astype(str)
    df["factor"] = df["factor"].astype(str)
    return df


def read_sample_sheet(path):
    """run_together-style sample sheet (TSV). Sample id from `id` (or `sample`/`sample_id`, else the basename of
    `in_dir`/`in_prefix`); optional `tissue` (or `tissue_context`) and `pseudobulk` (path relative to the sheet).
    Other columns (input roles) are ignored. '', '-', '.', 'NA' mean unset."""
    base = os.path.dirname(os.path.abspath(path))
    out = []
    with open_text(path) as fh:
        for i, row in enumerate(csv.DictReader(fh, delimiter="\t"), 2):
            row = {k.strip(): (v or "").strip() for k, v in row.items() if k is not None}
            row = {k: v for k, v in row.items() if v not in SHEET_UNSET}
            if not row:
                continue
            sid = next((row[c] for c in SHEET_ID_COLUMNS if c in row), None)
            if not sid and (row.get("in_dir") or row.get("in_prefix")):
                sid = os.path.basename(os.path.normpath(row.get("in_dir") or row["in_prefix"]))
            if not sid:
                raise ValueError(f"{path}: line {i} has no sample id (columns {', '.join(SHEET_ID_COLUMNS)})")
            pb = row.get("pseudobulk")
            out.append({"name": sid, "tissue": next((row[c] for c in SHEET_TISSUE_COLUMNS if c in row), None),
                        "pseudobulk": (pb if os.path.isabs(pb) else os.path.join(base, pb)) if pb else None})
    if len({s["name"] for s in out}) != len(out):
        raise ValueError(f"Sample ids in {path} are not unique")
    return out


def strip_common_prefix(names):
    """<multi_id>-<sample> directory names -> sample ids (the shared '<multi_id>-' is dropped)."""
    if len(names) < 2:
        return list(names)
    cut = os.path.commonprefix(names).rfind("-") + 1
    return [n[cut:] for n in names] if cut > 0 and all(len(n) > cut for n in names) else list(names)


def discover_sample_dirs(prefix, multi_catalog):
    """Sample id -> directory for a multi-sample prefix, from multi-catalog.yaml or the sub-directories that hold
    a pseudobulk for the prefix. Empty for a single-sample prefix."""
    root, base = os.path.dirname(os.path.abspath(prefix)), os.path.basename(prefix)
    mc = os.path.join(root, multi_catalog)
    if os.path.exists(mc):
        import yaml
        with open(mc, encoding="utf-8") as fh:
            cat = yaml.safe_load(fh) or {}
        found = {str(sid): os.path.join(root, os.path.dirname(rel)) for sid, rel in (cat.get("samples") or {}).items()}
        if found:
            return found
    dirs = sorted(d for d in os.listdir(root)
                  if os.path.isdir(os.path.join(root, d)) and first_existing(os.path.join(root, d, base), PSEUDOBULK_SUFFIXES))
    return dict(zip(strip_common_prefix(dirs), (os.path.join(root, d) for d in dirs)))


def resolve_samples(args):
    """[{name, tissue, pseudobulk}] for a multi-sample run, [] for a single-sample run."""
    sheet = read_sample_sheet(args.sample_sheet) if args.sample_sheet else None
    if args.single_sample:
        if sheet:
            raise ValueError("--single-sample and --sample-sheet are mutually exclusive")
        return []
    dirs = discover_sample_dirs(args.prefix, args.multi_catalog) if args.prefix else {}
    if sheet is None:
        sheet = [{"name": n, "tissue": None, "pseudobulk": None} for n in dirs]
    base = os.path.basename(args.prefix) if args.prefix else None
    samples = []
    for s in sheet:
        pb = s["pseudobulk"]
        if not pb:
            hits = [d for n, d in dirs.items() if s["name"] in (n, os.path.basename(d))
                    or os.path.basename(d).endswith("-" + s["name"])]
            if len(hits) != 1 or base is None:
                raise ValueError(f"Sample {s['name']!r}: no `pseudobulk` column in the sheet and "
                                 f"{'no' if not hits else 'several'} matching sample directories "
                                 f"(found: {', '.join(sorted(dirs)) or 'none'}; needs --prefix)")
            pb = first_existing(os.path.join(hits[0], base), PSEUDOBULK_SUFFIXES)
        if not pb or not os.path.exists(pb):
            raise FileNotFoundError(f"Sample {s['name']!r}: pseudobulk not found ({pb})")
        tissue = s["tissue"] or args.tissue
        if not tissue:
            raise ValueError(f"No tissue for sample {s['name']!r}: add a `tissue` column to --sample-sheet or pass --tissue")
        samples.append({"name": s["name"], "tissue": tissue, "pseudobulk": pb})
    return samples


def factor_sort_key(f):
    return (0, int(f), "") if f.lstrip("-").isdigit() else (1, 0, f)


# -----------------------------
# Statistics (cartloader bulk DE: 2x2 table + 0.5 per cell; Chi2 and odds-ratio fold change)
# -----------------------------
def _chi2_sf_1df(x):
    return math.erfc(math.sqrt(max(x, 0.0) / 2.0))


def _chi2_threshold(p):
    lo, hi = 0.0, 1000.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if _chi2_sf_1df(mid) > p else (lo, mid)
    return hi


def de_stats(values):
    """Chi2 and odds-ratio fold change of every row (gene) in every column (factor) against all other columns."""
    gene_total = values.sum(axis=1, keepdims=True)
    factor_total = values.sum(axis=0, keepdims=True)
    grand_total = values.sum()
    a = values + PSEUDOCOUNT
    b = (gene_total - values) + PSEUDOCOUNT
    c = (factor_total - values) + PSEUDOCOUNT
    d = (grand_total - factor_total - gene_total + values) + PSEUDOCOUNT
    n = a + b + c + d
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = n * (a * d - b * c) ** 2 / ((a + b) * (c + d) * (a + c) * (b + d))
        fold_change = (a * d) / (b * c)
    return chi2, fold_change, gene_total


def compute_de(matrix, min_fold_change=1.5, max_pval=CHI2_PVAL_MAX, min_count=1.0):
    values = matrix.to_numpy(dtype=float)
    chi2, fold_change, gene_total = de_stats(values)
    keep = (fold_change >= min_fold_change) & (chi2 > _chi2_threshold(max_pval)) & (values >= min_count)
    gi, fi = np.nonzero(keep)
    de = pd.DataFrame({
        "gene": matrix.index.to_numpy()[gi],
        "factor": np.asarray(matrix.columns)[fi].astype(str),
        "Chi2": chi2[gi, fi],
        "pval": [_chi2_sf_1df(x) for x in chi2[gi, fi]],
        "FoldChange": fold_change[gi, fi],
        "gene_total": gene_total[gi, 0],
    })
    order = pd.Categorical(de["factor"], categories=[str(c) for c in matrix.columns], ordered=True)
    return de.assign(_o=order).sort_values(["_o", "Chi2"], ascending=[True, False]).drop(columns="_o")


def marker_genes(de, primary="Chi2", secondary="FoldChange", top_n=10, factors=None):
    """Top N by primary and top N by secondary, interleaved P1, S1, P2, S2, ... without duplicates."""
    for col in (primary, secondary):
        if col not in de.columns:
            raise ValueError(f"Ranking column {col!r} not in DE table; has {list(de.columns)}")
    grouped = {str(f): sub for f, sub in de.groupby("factor", sort=False)}
    out = {}
    for f in factors if factors is not None else sorted(grouped, key=factor_sort_key):
        sub = grouped.get(str(f))
        if sub is None:
            out[str(f)] = []
            continue
        prim = sub.sort_values(primary, ascending=False, kind="stable")["gene"].astype(str).head(top_n).tolist()
        sec = sub.sort_values(secondary, ascending=False, kind="stable")["gene"].astype(str).head(top_n).tolist()
        merged = []
        for i in range(top_n):
            for lst in (prim, sec):
                if i < len(lst) and lst[i] not in merged:
                    merged.append(lst[i])
        out[str(f)] = merged
    return out


# -----------------------------
# Evidence
# -----------------------------
def factor_summaries(de, factors, matrix=None, trace=0.001, minor=0.005, weak=30, moderate=150):
    n_de = de.groupby("factor").size()
    totals = matrix.sum(axis=0) if matrix is not None else None
    share = (totals / totals.sum()) if matrix is not None else None
    out = {}
    for f in factors:
        n = int(n_de.get(f, 0))
        d = {"de_support": "weak" if n <= weak else "moderate" if n <= moderate else "strong", "n_de": n}
        if share is not None:
            s = float(share[f])
            d["abundance"] = "trace" if s < trace else "minor" if s < minor else "major"
            d["share"] = s
            d["total"] = float(totals[f])
        out[f] = d
    return out


def specificity_categories(matrix, genes_by_factor, exclusive=0.5):
    freq = matrix / matrix.sum(axis=0).replace(0, np.nan)
    share = freq.div(freq.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    rank = freq.rank(axis=1, ascending=False, method="min")
    out = {}
    for f, genes in genes_by_factor.items():
        cats = {}
        for g in genes:
            if g not in share.index or f not in share.columns:
                cats[g] = "shared"
                continue
            s, r = float(share.at[g, f]), float(rank.at[g, f])
            cats[g] = "exclusive" if s >= exclusive else "highest" if r <= 1 else "top3" if r <= 3 else "shared"
        out[f] = cats
    return out


def read_sample_pseudobulks(samples, factors):
    """Each sample's pixel pseudobulk, aligned to the factor columns of the model."""
    matrices = {}
    for s in samples:
        pb = read_gene_by_factor(s["pseudobulk"])
        missing = [f for f in factors if f not in pb.columns]
        if missing:
            log.warning("%s: %d factor column(s) missing from %s (treated as 0)", s["name"], len(missing), s["pseudobulk"])
        matrices[s["name"]] = pb.reindex(columns=factors, fill_value=0.0)
    return matrices


def factor_distributions(samples, matrices, present_share=0.1, specific_share=0.8):
    """Where each factor lives: share of each sample's pixel pseudobulk, normalised across samples."""
    weights = pd.DataFrame({name: pb.sum(axis=0) / pb.to_numpy().sum() for name, pb in matrices.items()})
    tissue_of = {s["name"]: s["tissue"] for s in samples}
    norm = weights.div(weights.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    out = {}
    for f, row in norm.iterrows():
        ranked = row.sort_values(ascending=False)
        holders = [s for s, w in ranked.items() if w >= present_share]
        tissues = list(dict.fromkeys(tissue_of[s] for s in holders))
        if ranked.iloc[0] <= 0:
            cat = "absent"
        elif ranked.iloc[0] >= specific_share:
            cat = "specific"
        elif len(tissues) <= 2:
            cat = "restricted"
        else:
            cat = "shared"
        out[str(f)] = {"category": cat, "samples": holders, "tissues": tissues}
    return out


def per_sample_markers(matrices, dists, top_n=4, max_samples=3):
    wanted = {}
    for f, d in dists.items():
        for s in d["samples"][:max_samples]:
            wanted.setdefault(s, set()).add(f)
    found = {}
    for name, pb in matrices.items():
        if name not in wanted:
            continue
        mk = marker_genes(compute_de(pb), top_n=top_n)
        for f in wanted[name]:
            if mk.get(f):
                found.setdefault(f, {})[name] = mk[f]
    return {f: {s: d[s] for s in dists[f]["samples"] if s in d} for f, d in found.items()}


def build_profiles(args, samples):
    matrix = read_gene_by_factor(args.model) if args.model else None
    de = read_de(args.de) if args.de else compute_de(matrix)
    factors = [str(c) for c in matrix.columns] if matrix is not None else sorted(de["factor"].unique(),
                                                                                 key=factor_sort_key)
    markers = marker_genes(de, args.primary_rank, args.secondary_rank, args.top_n, factors)
    dists, by_sample, sample_totals = {}, {}, None
    if samples:
        matrices = read_sample_pseudobulks(samples, factors)
        sample_totals = pd.DataFrame({name: pb.sum(axis=0) for name, pb in matrices.items()})  # factor x sample
        dists = factor_distributions(samples, matrices)
        if not args.no_sample_markers:
            by_sample = per_sample_markers(matrices, dists)
        tissues = list(dict.fromkeys(s["tissue"] for s in samples))
    else:
        tissues = [args.tissue]
    spec = specificity_categories(matrix, markers) if (matrix is not None and not args.no_specificity) else {}
    summaries = factor_summaries(de, factors, matrix)
    profiles = []
    for f in factors:
        ev = dict(summaries.get(f, {}))
        if f in dists:
            d = dists[f]
            ev.update({"distribution": d["category"], "found_in_samples": d["samples"], "found_in_tissues": d["tissues"]})
        if f in spec:
            ev["marker_specificity"] = spec[f]
        if by_sample.get(f):
            ev["markers_by_sample"] = by_sample[f]
        profiles.append({"factor": f, "marker_genes": markers.get(f, []), "evidence": ev})
    genes = set(matrix.index) if matrix is not None else set(de["gene"])
    return SimpleNamespace(profiles=profiles, tissues=tissues, genes=genes, matrix=matrix, sample_totals=sample_totals)


# -----------------------------
# Candidate labels (optional; files written by generate_candidate_labels_with_llm.py)
# -----------------------------
def load_candidates(paths, panel_genes):
    """Union of candidate-set JSON files by label, markers restricted to genes measured in this dataset."""
    merged = {}
    for p in paths:
        with open(p, encoding="utf-8") as fh:
            for c in json.load(fh)["candidates"]:
                if c["label"] in merged:
                    m = merged[c["label"]]
                    m["tissues"] = list(dict.fromkeys(m["tissues"] + c.get("tissues", [])))
                    m["canonical_markers"] = list(dict.fromkeys(m["canonical_markers"] + c["canonical_markers"]))[:12]
                else:
                    merged[c["label"]] = {**c, "tissues": list(c.get("tissues", []))}
    for c in merged.values():
        c["canonical_markers"] = [g for g in c["canonical_markers"] if g in panel_genes] or c["canonical_markers"]
    return list(merged.values())


def generate_candidates(args, tissues, out_path):
    """--candidates auto: one cached candidate list per tissue via generate_candidate_labels_with_llm."""
    from cartloader.scripts.generate_candidate_labels_with_llm import generate
    gargs = SimpleNamespace(
        organism=args.organism, mode=args.mode, cache_dir=os.path.expanduser(args.candidate_cache_dir),
        no_generate=False, min_candidates=80, max_candidates=200, api_type=args.api_type,
        model_name=args.model_name, effort=args.candidate_effort, max_output_tokens=64000,
        no_fallbacks=args.no_fallbacks, api_base_url=args.api_base_url, request_timeout=args.request_timeout)
    with ThreadPoolExecutor(max_workers=max(1, min(args.threads, len(tissues)))) as ex:
        sets = list(ex.map(lambda t: generate(gargs, t), tissues))
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump({"schema_version": 1, "candidates": [c for s in sets for c in s["candidates"]]}, fh, indent=1)
    return out_path


# -----------------------------
# Prompt
# -----------------------------
def render(template, **values):
    system, user = template.split("# USER", 1)
    system = system.replace("# SYSTEM", "", 1).strip()
    return Template(system).substitute(values), Template(user.strip()).substitute(values)


def factor_line(p, sample_tissue):
    ev = p["evidence"]
    d = {"id": p["factor"], "marker_genes": p["marker_genes"]}
    for k in ("abundance", "de_support"):
        if k in ev:
            d[k] = ev[k]
    if "distribution" in ev:
        d["distribution"] = ev["distribution"]
        d["found_in"] = [f"{s} ({sample_tissue.get(s, '?')})" for s in ev.get("found_in_samples", [])]
    if "marker_specificity" in ev:
        d["marker_specificity"] = ev["marker_specificity"]
    if "markers_by_sample" in ev:
        d["markers_by_sample"] = ev["markers_by_sample"]
    return json.dumps(d, ensure_ascii=False)


def digest_line(p):
    ev = p["evidence"]
    tags = [ev[k] for k in ("abundance", "de_support") if k in ev]
    if "distribution" in ev:
        tags.append(f"{ev['distribution']}: {', '.join(ev.get('found_in_samples', [])[:3])}")
    return f"- {p['factor']}: {', '.join(p['marker_genes'][:6])}" + (f" [{'; '.join(tags)}]" if tags else "")


def make_batches(profiles, size):
    """Group factors by their strongest sample (then id) so a batch holds related factors. size <= 0: one batch."""
    if size <= 0 or size >= len(profiles):
        return [profiles] if profiles else []

    def key(p):
        samples = p["evidence"].get("found_in_samples") or [""]
        return (samples[0], len(samples))
    ordered = sorted(profiles, key=key) if any("found_in_samples" in p["evidence"] for p in profiles) else profiles
    return [ordered[i:i + size] for i in range(0, len(ordered), size)]


def response_schema(labels=None):
    label = {"type": "string"} if labels is None else {"type": "string", "enum": labels + [UNRESOLVED]}
    confidence = {"type": "string", "enum": ["high", "medium", "low"]}
    genes = {"type": "array", "items": {"type": "string"}}
    alternative = {
        "type": "object",
        "properties": {"label": label, "kind": {"type": "string", "enum": ["cell_type", "program"]},
                       "confidence": confidence, "genes": genes},
        "required": ["label", "kind", "confidence", "genes"],
        "additionalProperties": False,
    }
    item = {
        "type": "object",
        "properties": {
            "factor": {"type": "string"},
            "label": label,
            "alias": {"type": "string"},
            "kind": {"type": "string", "enum": ["cell_type", "program", "unresolved"]},
            "compartment": {"type": "string", "enum": list(COMPARTMENTS)},
            "confidence": confidence,
            "key_genes": genes,
            "alternatives": {"type": "array", "items": alternative},
            "rationale": {"type": "string"},
            "in_candidate_set": {"type": "boolean"},
        },
        "required": ["factor", "label", "alias", "kind", "compartment", "confidence", "key_genes", "alternatives",
                     "rationale", "in_candidate_set"],
        "additionalProperties": False,
    }
    return {"type": "object",
            "properties": {"dataset_summary": {"type": "string"}, "annotations": {"type": "array", "items": item}},
            "required": ["dataset_summary", "annotations"], "additionalProperties": False}


def prompt_context(usable, samples, tissues, candidates, organism, mode, closed):
    notes = []
    if any("de_support" in p["evidence"] for p in usable):
        notes.append(EVIDENCE_NOTES_SUMMARY)
    if samples:
        notes.append(EVIDENCE_NOTES_MULTI)
    if any("marker_specificity" in p["evidence"] for p in usable):
        notes.append(EVIDENCE_NOTES_SPECIFICITY)
    if any("markers_by_sample" in p["evidence"] for p in usable):
        notes.append(EVIDENCE_NOTES_SAMPLE_MARKERS)
    labels = [c["label"] for c in candidates]
    cand_block = ("\n".join(f"- {c['label']} [{c['kind']}]: {', '.join(c['canonical_markers'][:8])}" for c in candidates)
                  if candidates else "(none provided)")
    if closed:
        if not candidates:
            raise ValueError("--closed requires --candidates")
        candidate_rule, label_rule = "a closed list", "Choose label only from the candidate labels (or Unresolved)."
    else:
        candidate_rule = "suggestions, not a closed list"
        label_rule = ("Use a candidate label when it fits. When none fits well, write a new, more accurate label in "
                      "the same style and set in_candidate_set to false.")
    tissue_summary = "; ".join(tissues)
    if samples and len(tissues) > 1:
        tissue_summary += " (Samples: " + "; ".join(f"{s['name']} = {s['tissue']}" for s in samples) + ")"
    elif samples:
        tissue_summary += " (Samples: " + ", ".join(s["name"] for s in samples) + ")"
    return {"organism": organism, "tissue_summary": tissue_summary, "mode_description": MODE_DESCRIPTIONS[mode],
            "evidence_notes": "\n".join(notes), "candidate_rule": candidate_rule, "candidates_block": cand_block,
            "label_rule": label_rule, "schema": response_schema(labels if closed else None),
            "sample_tissue": {s["name"]: s["tissue"] for s in samples}, "usable": usable}


def build_request(ctx, batch, name):
    """(name, system, user, schema, factor ids) for one batch; factors outside the batch are listed for contrast."""
    ids = {p["factor"] for p in batch}
    others = [p for p in ctx["usable"] if p["factor"] not in ids]
    digest = ""
    if others:
        digest = ("\nOther factors in the same dataset, for contrast only (do not annotate these):\n"
                  + "\n".join(digest_line(p) for p in others) + "\n")
    system, user = render(
        JOINT_PROMPT, organism=ctx["organism"], tissue_summary=ctx["tissue_summary"],
        mode_description=ctx["mode_description"], evidence_notes=ctx["evidence_notes"],
        candidate_rule=ctx["candidate_rule"], candidates_block=ctx["candidates_block"], label_rule=ctx["label_rule"],
        factors_block="\n".join(factor_line(p, ctx["sample_tissue"]) for p in batch), digest_block=digest)
    return (name, system, user, ctx["schema"], [p["factor"] for p in batch])


# -----------------------------
# LLM calls (JSON-schema constrained output)
# -----------------------------
def call_claude(system, prompt, schema, model, effort, max_tokens, fallbacks):
    import anthropic  # only needed for --api-type claude

    client = anthropic.Anthropic()
    kwargs = {"model": model, "max_tokens": max_tokens, "messages": [{"role": "user", "content": prompt}],
              "output_config": {"effort": effort, "format": {"type": "json_schema", "schema": schema}}}
    if "haiku" not in model:
        kwargs["thinking"] = {"type": "adaptive"}  # explicit: models before Opus 5 do not think unless asked
    if system:
        kwargs["system"] = system
    stream_cm = (client.beta.messages.stream(betas=[ANTHROPIC_FALLBACK_BETA], extra_body={"fallbacks": "default"},
                                             **kwargs) if fallbacks else client.messages.stream(**kwargs))
    start = last = time.monotonic()
    with stream_cm as stream:
        for _ in stream:
            now = time.monotonic()
            if now - last >= 30:
                last = now
                log.info("  still generating (%.0fs)", now - start)
        message = stream.get_final_message()
    if message.stop_reason == "refusal":
        raise RuntimeError(f"Claude declined the request: {getattr(message, 'stop_details', None)}")
    if message.stop_reason == "max_tokens":
        raise OutputTruncated(f"Claude hit max_tokens={max_tokens}")
    text = next(b.text for b in message.content if b.type == "text")
    usage = {"input_tokens": message.usage.input_tokens, "output_tokens": message.usage.output_tokens}
    return json.loads(text), message.model, usage


def call_openai_responses(system, prompt, schema, model, effort, max_tokens, base_url, api_key, timeout, retries):
    url = base_url.rstrip("/")
    url = url if url.endswith("/responses") else url + "/responses"
    payload = {"model": model, "input": prompt, "max_output_tokens": max_tokens,
               "text": {"format": {"type": "json_schema", "name": "annotations", "schema": schema, "strict": True}}}
    if system:
        payload["instructions"] = system
    if effort in ("low", "medium", "high"):
        payload["reasoning"] = {"effort": effort}
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    for attempt in range(retries):
        r = requests.post(url, headers=headers, json=payload, timeout=timeout)
        if r.status_code == 400 and "reasoning" in payload and "reasoning" in r.text:
            payload.pop("reasoning")  # model does not accept reasoning.effort
            continue
        if r.status_code in (429, 500, 502, 503, 504, 529) and attempt < retries - 1:
            time.sleep(2 ** attempt * 5)
            continue
        r.raise_for_status()
        break
    data = r.json()
    if data.get("status") == "incomplete":
        details = data.get("incomplete_details") or {}
        if details.get("reason") == "max_output_tokens":
            raise OutputTruncated(f"response hit max_output_tokens={max_tokens}")
        raise RuntimeError(f"response incomplete ({details})")
    text = data.get("output_text") or "".join(c.get("text", "") for item in data.get("output", []) or []
                                               for c in item.get("content", []) or [] if c.get("type") == "output_text")
    return json.loads(text), data.get("model", model), data.get("usage", {})


def call_gemini(system, prompt, schema, model, api_key, timeout, retries):
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent"
    payload = {"contents": [{"role": "user", "parts": [{"text": prompt}]}],
               "generationConfig": {"responseMimeType": "application/json", "responseJsonSchema": schema}}
    if system:
        payload["systemInstruction"] = {"parts": [{"text": system}]}
    for attempt in range(retries):
        r = requests.post(url, headers={"x-goog-api-key": api_key}, json=payload, timeout=timeout)
        if r.status_code in (429, 500, 503) and attempt < retries - 1:
            time.sleep(2 ** attempt * 5)
            continue
        r.raise_for_status()
        break
    data = r.json()
    cand = data["candidates"][0]
    if cand.get("finishReason") == "MAX_TOKENS":
        raise OutputTruncated("Gemini hit its output token limit")
    text = "".join(p.get("text", "") for p in cand["content"]["parts"])
    return json.loads(text), model, data.get("usageMetadata", {})


def complete_json(args, system, prompt, schema):
    if args.api_type == "claude":
        return call_claude(system, prompt, schema, args.model_name, args.effort, args.max_output_tokens,
                           not args.no_fallbacks)
    key = os.environ.get(API_KEY_ENV[args.api_type], "")
    if args.api_type in ("openai", "umgpt"):
        return call_openai_responses(system, prompt, schema, args.model_name, args.effort, args.max_output_tokens,
                                     args.api_base_url or DEFAULT_BASE_URL[args.api_type], key,
                                     args.request_timeout, args.max_retries)
    return call_gemini(system, prompt, schema, args.model_name, key, args.request_timeout, args.max_retries)


def run_request(args, llm_dir, ctx, req):
    """Send one request (or reuse its saved response). If the answer overflows the output budget, split the batch
    in halves and retry each. Returns {"annotations", "summaries", "calls"}."""
    name, system, user, schema, ids = req
    key = hashlib.sha256(json.dumps([args.api_type, args.model_name, args.effort, system, user, schema]).encode()
                         ).hexdigest()[:16]
    path = os.path.join(llm_dir, f"{name}.{key}.json")
    split_marker = os.path.join(llm_dir, f"{name}.{key}.split")
    call = {"name": name, "n_factors": len(ids), "file": os.path.basename(path), "reused": False}

    def split(reason):
        log.warning("%s: %s; splitting the %d factors into two requests", name, reason, len(ids))
        by_id = {p["factor"]: p for p in ctx["usable"]}
        batch = [by_id[i] for i in ids]
        half = len(batch) // 2
        parts = [run_request(args, llm_dir, ctx, build_request(ctx, b, f"{name}{s}"))
                 for s, b in (("a", batch[:half]), ("b", batch[half:]))]
        return {k: parts[0][k] + parts[1][k] for k in ("annotations", "summaries", "calls")}

    if os.path.exists(split_marker) and not args.no_reuse and len(ids) >= 2:
        return split("overflowed the output budget in an earlier run")
    if os.path.exists(path) and not args.no_reuse:
        log.info("%s: reusing saved response %s", name, path)
        with open(path, encoding="utf-8") as fh:
            saved = json.load(fh)
        call.update(model=saved.get("model"), usage=saved.get("usage"), reused=True)
        data = saved["data"]
    else:
        log.info("%s: sending %d factors to %s (%s, effort %s, max %d output tokens)", name, len(ids), args.api_type,
                 args.model_name, args.effort, args.max_output_tokens)
        try:
            data, model, usage = complete_json(args, system, user, schema)
        except OutputTruncated as e:
            if len(ids) < 2:
                raise
            with open(split_marker, "w", encoding="utf-8") as fh:
                fh.write(f"{e}\n")
            return split(str(e))
        with open(path, "w", encoding="utf-8") as fh:
            json.dump({"model": model, "api_type": args.api_type, "effort": args.effort, "usage": usage,
                       "prompt_version": PROMPT_VERSION, "system": system,
                       "prompt": user, "data": data}, fh, indent=1)
        call.update(model=model, usage=usage)
        log.info("%s: done (%s)", name, usage)
    got = {str(a["factor"]) for a in data["annotations"]}
    missing, extra = set(ids) - got, got - set(ids)
    if missing:
        log.warning("%s: model skipped factors %s; they are left Unresolved", name, sorted(missing, key=factor_sort_key))
    if extra:
        log.warning("%s: ignoring annotations for factors not in the request: %s", name, sorted(extra, key=factor_sort_key))
    anns = [a for a in data["annotations"] if str(a["factor"]) in set(ids)]
    summary = (data.get("dataset_summary") or "").strip()
    return {"annotations": anns, "summaries": [summary] if summary else [], "calls": [call]}


def load_saved_responses(llm_dir, factor_ids):
    """--report-only: annotations from the responses saved in llm_dir, newest file first, so each factor takes its
    most recent answer even if the prompt has changed since. Returns (raw, summaries, calls, models, efforts)."""
    if not os.path.isdir(llm_dir):
        raise FileNotFoundError(f"No saved LLM responses: {llm_dir} does not exist (run once without --report-only)")
    files = sorted((os.path.join(llm_dir, f) for f in os.listdir(llm_dir) if f.endswith(".json")),
                   key=os.path.getmtime, reverse=True)
    raw, summaries, calls, models, efforts = {}, [], [], [], []
    for path in files:
        try:
            with open(path, encoding="utf-8") as fh:
                saved = json.load(fh)
            anns = saved["data"]["annotations"]
        except (OSError, ValueError, KeyError, TypeError) as e:
            log.warning("Skipping %s: not a saved response (%s)", path, e)
            continue
        new = [a for a in anns if str(a.get("factor")) in factor_ids and str(a["factor"]) not in raw]
        if not new:
            continue
        raw.update({str(a["factor"]): a for a in new})
        summary = (saved["data"].get("dataset_summary") or "").strip()
        if summary:
            summaries.append(summary)
        calls.append({"name": os.path.basename(path).split(".")[0], "n_factors": len(new), "file": os.path.basename(path),
                      "reused": True, "model": saved.get("model"), "usage": saved.get("usage")})
        models.append(saved.get("model"))
        efforts.append(saved.get("effort"))
        log.info("%s: %d factor annotation(s)", path, len(new))
    if not raw:
        raise FileNotFoundError(f"No saved LLM responses for these factors in {llm_dir}")
    return raw, summaries, calls, [m for m in dict.fromkeys(models) if m], [e for e in dict.fromkeys(efforts) if e]


# -----------------------------
# Output
# -----------------------------
def normalize_alias(raw):
    if not raw:
        return "Unknown"
    s = raw.strip().strip("`\"'").splitlines()[0].strip()
    if " " in s:
        s = "".join(p[:1].upper() + p[1:] for p in re.split(r"\s+", s) if p)
    s = re.sub(r"[^A-Za-z0-9+\-/]", "", s)
    return s or "Unknown"


def suffix_duplicates(pairs):
    counts = {}
    for _, a in pairs:
        counts[a] = counts.get(a, 0) + 1
    seen, out = {}, []
    for f, a in pairs:
        if counts[a] > 1:
            seen[a] = seen.get(a, 0) + 1
            a = f"{a}_{seen[a]}"
        out.append((f, a))
    return out


def is_unresolved(a):
    return a["label"] == UNRESOLVED or a["kind"] == "unresolved"


def to_annotations(profiles, usable, raw, candidate_labels, pending=None):
    usable_ids = {p["factor"] for p in usable}
    rows = []
    for p in profiles:
        f, ev = p["factor"], p["evidence"]
        base = {"factor": f, "marker_genes": ",".join(p["marker_genes"]), "abundance": ev.get("abundance", ""),
                "share": round(ev["share"], 6) if "share" in ev else "", "de_support": ev.get("de_support", ""),
                "n_de": ev.get("n_de", ""), "distribution": ev.get("distribution", ""),
                "found_in": ",".join(ev.get("found_in_samples", []))}
        a = raw.get(f)
        if a is None:
            why = (pending if pending and f in usable_ids else
                   "too few marker genes" if f not in usable_ids else "not returned by the model")
            rows.append({**base, "status": "pending" if pending and f in usable_ids else "unfitted", "label": "",
                         "alias": UNRESOLVED, "kind": "", "compartment": "", "confidence": "", "key_genes": "",
                         "alternatives": "", "rationale": why, "in_candidate_set": ""})
            continue
        unresolved = is_unresolved(a)
        status = "unfitted" if unresolved else ("ambiguous" if a["confidence"] == "low" else "assigned")
        alts = "; ".join(f"{x['label']} [{x['kind']}, {x['confidence']}]" for x in a.get("alternatives", []))
        rows.append({**base, "status": status, "label": "" if unresolved else a["label"],
                     "alias": UNRESOLVED if unresolved else normalize_alias(a["alias"]),
                     "kind": "" if unresolved else a["kind"],
                     "compartment": "" if unresolved else a.get("compartment", ""),
                     "confidence": CONFIDENCE_LEVELS[a["confidence"]], "key_genes": ",".join(a.get("key_genes", [])),
                     "alternatives": alts, "rationale": a["rationale"],
                     "in_candidate_set": (a["label"] in candidate_labels) if candidate_labels else ""})
    cols = ["factor", "status", "alias", "label", "kind", "compartment", "confidence", "key_genes", "alternatives",
            "rationale", "in_candidate_set", "abundance", "share", "de_support", "n_de", "distribution", "found_in",
            "marker_genes"]
    return pd.DataFrame(rows)[cols]


# -----------------------------
# Report
# -----------------------------
def _b64(arr, dtype):
    return base64.b64encode(np.ascontiguousarray(arr, dtype=dtype).tobytes()).decode("ascii")


def _lfc_i2(fold_change):
    with np.errstate(divide="ignore", invalid="ignore"):
        lfc = np.nan_to_num(np.log2(fold_change), nan=0.0, posinf=327.0, neginf=-327.0)
    return np.clip(np.round(lfc * 100), -32767, 32767)


def dotplot_payload(matrix, genes):
    """Count and log2 fold change (odds ratio vs all other factors, as bulk DE) for genes x factors."""
    values = matrix.to_numpy(dtype=float)
    _, fold_change, _ = de_stats(values)
    idx = matrix.index.get_indexer(genes)
    return {"genes": list(genes), "counts": _b64(values[idx], "<f4"), "lfc": _b64(_lfc_i2(fold_change[idx]), "<i2")}


def report_payload(bundle, table, raw, samples, summaries, calls, candidates, meta, gene_mode, colors=None):
    profiles = bundle.profiles
    colors = colors or {}
    rows = {r["factor"]: r for r in table.to_dict("records")}
    factors = []
    for p in profiles:
        f, ev, r, a = p["factor"], p["evidence"], rows[p["factor"]], raw.get(p["factor"]) or {}
        factors.append({
            "id": f, "color": colors.get(f), "total": ev.get("total"), "status": r["status"], "label": r["label"], "alias": r["alias"], "kind": r["kind"],
            "compartment": r["compartment"], "confidence": a.get("confidence", ""),
            "key_genes": list(a.get("key_genes", [])),
            "alternatives": [{"label": x["label"], "kind": x["kind"], "confidence": x["confidence"],
                              "genes": list(x.get("genes", []))} for x in a.get("alternatives", [])],
            "rationale": r["rationale"], "in_candidate_set": r["in_candidate_set"], "markers": p["marker_genes"],
            "spec": ev.get("marker_specificity", {}), "share": ev.get("share"), "abundance": ev.get("abundance", ""),
            "de_support": ev.get("de_support", ""), "n_de": ev.get("n_de"), "distribution": ev.get("distribution", ""),
            "found_in": ev.get("found_in_samples", []), "markers_by_sample": ev.get("markers_by_sample", {}),
        })
    dot = None
    if bundle.matrix is not None:
        if gene_mode == "all":
            genes = list(bundle.matrix.index)
        else:
            wanted = []
            for fa in factors:
                wanted += fa["markers"] + fa["key_genes"] + [g for x in fa["alternatives"] for g in x["genes"]]
            wanted += [g for c in candidates for g in c.get("canonical_markers", [])]
            genes = [g for g in dict.fromkeys(wanted) if g in bundle.genes]
        dot = dotplot_payload(bundle.matrix, genes)
    sdot = None
    if bundle.sample_totals is not None:
        st = bundle.sample_totals.reindex(index=[fa["id"] for fa in factors], columns=[s["name"] for s in samples])
        values = st.fillna(0.0).to_numpy(dtype=float)
        _, fold_change, _ = de_stats(values)
        sdot = {"counts": _b64(values, "<f4"), "lfc": _b64(_lfc_i2(fold_change), "<i2")}
    return {"meta": meta, "compartments": list(COMPARTMENTS), "factors": factors, "summaries": summaries,
            "samples": [{"name": s["name"], "tissue": s["tissue"]} for s in samples], "dot": dot, "sdot": sdot,
            "calls": calls}


def render_html(payload):
    data = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    title = payload["meta"]["title"].replace("&", "&amp;").replace("<", "&lt;")
    return HTML_TEMPLATE.replace("__TITLE__", title).replace("__DATA__", data)


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  *, *::before, *::after { box-sizing: border-box; }
  [hidden] { display: none !important; }
  :root {
    color-scheme: light;
    --bg: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --muted: #6f6d68; --grid: #e1e0d9;
    --axis: #c3c2b7; --border: rgba(11,11,11,.12); --row-alt: #f3f2ee; --hover: rgba(42,120,214,.10);
    --accent: #2a78d6; --accent-ink: #1c5cab; --accent-soft: #e1ecfa;
    --good: #006300; --good-soft: #e2f0e2; --warn: #875800; --warn-soft: #faefd6; --none: #66655f; --none-soft: #ebeae5;
    --c1: #2a78d6; --c2: #eb6834; --c3: #1baf7a; --c4: #eda100; --c5: #e87ba4; --c6: #008300; --c7: #4a3aa7; --c8: #e34948;
    --c-unres: #b9b7b0; --div-neg: #1c5cab; --div-mid: #f0efec; --div-pos: #c43838;
    --bar: #b3b1aa; --bar-on: #52514e;
    --sans: system-ui, -apple-system, "Segoe UI", sans-serif;
    --mono: ui-monospace, "SFMono-Regular", Menlo, Consolas, monospace;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      color-scheme: dark;
      --bg: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink2: #c3c2b7; --muted: #9a988f; --grid: #2c2c2a;
      --axis: #383835; --border: rgba(255,255,255,.12); --row-alt: #202020; --hover: rgba(57,135,229,.18);
      --accent: #3987e5; --accent-ink: #86b6ef; --accent-soft: #1a2a40;
      --good: #0ca30c; --good-soft: #15301a; --warn: #e0ad52; --warn-soft: #3a2e17; --none: #a3a198; --none-soft: #2a2a28;
      --c1: #3987e5; --c2: #d95926; --c3: #199e70; --c4: #c98500; --c5: #d55181; --c6: #008300; --c7: #9085e9; --c8: #e66767;
      --c-unres: #5d5c58; --div-neg: #3987e5; --div-mid: #383835; --div-pos: #e66767;
      --bar: #5f5e59; --bar-on: #c3c2b7;
    }
  }
  :root[data-theme="dark"] {
    color-scheme: dark;
    --bg: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink2: #c3c2b7; --muted: #9a988f; --grid: #2c2c2a;
    --axis: #383835; --border: rgba(255,255,255,.12); --row-alt: #202020; --hover: rgba(57,135,229,.18);
    --accent: #3987e5; --accent-ink: #86b6ef; --accent-soft: #1a2a40;
    --good: #0ca30c; --good-soft: #15301a; --warn: #e0ad52; --warn-soft: #3a2e17; --none: #a3a198; --none-soft: #2a2a28;
    --c1: #3987e5; --c2: #d95926; --c3: #199e70; --c4: #c98500; --c5: #d55181; --c6: #008300; --c7: #9085e9; --c8: #e66767;
    --c-unres: #5d5c58; --div-neg: #3987e5; --div-mid: #383835; --div-pos: #e66767;
    --bar: #5f5e59; --bar-on: #c3c2b7;
  }
  body { margin: 0; background: var(--bg); color: var(--ink); font-family: var(--sans); font-size: 14.5px; line-height: 1.5; }
  .page { max-width: 1440px; margin: 0 auto; padding: 24px 16px 60px; display: grid; grid-template-columns: minmax(0, 1fr); gap: 18px; }
  header { display: flex; flex-wrap: wrap; gap: 8px 24px; align-items: flex-start; justify-content: space-between; }
  header .hd { display: grid; gap: 4px; min-width: 0; flex: 1 1 320px; }
  h1, .meta { overflow-wrap: anywhere; }
  .eyebrow { font-size: 12px; letter-spacing: .06em; text-transform: uppercase; color: var(--muted); }
  h1 { font-size: clamp(22px, 3.2vw, 30px); line-height: 1.15; margin: 0; font-weight: 650; word-break: break-word; }
  h2 { font-size: 16px; margin: 0 0 8px; font-weight: 600; }
  p { margin: 0; max-width: 110ch; }
  .meta, .muted { color: var(--muted); }
  .meta { font-size: 13.5px; }
  button, select, input, textarea { font: inherit; font-size: 13.5px; color: var(--ink); }
  button, select, input[type=search], input[type=number], textarea { background: var(--surface); border: 1px solid var(--border); border-radius: 4px; padding: 5px 9px; }
  button { cursor: pointer; }
  button:hover { border-color: var(--axis); }
  :focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }
  .tabs { display: flex; gap: 2px; border-bottom: 1px solid var(--border); overflow-x: auto; }
  .tabs [role=tab] { border: 0; border-bottom: 2px solid transparent; border-radius: 0; background: none; padding: 8px 14px; color: var(--ink2); white-space: nowrap; }
  .tabs [role=tab][aria-selected=true] { color: var(--ink); border-bottom-color: var(--accent); font-weight: 600; }
  section[role=tabpanel] { display: grid; grid-template-columns: minmax(0, 1fr); gap: 18px; min-width: 0; }
  .card { background: var(--surface); border: 1px solid var(--border); border-radius: 6px; padding: 16px; min-width: 0; }
  .stats { display: flex; flex-wrap: wrap; gap: 10px; }
  .stat { background: var(--surface); border: 1px solid var(--border); border-radius: 6px; padding: 10px 14px; min-width: 120px; }
  .stat b { display: block; font-size: 22px; font-weight: 600; }
  .stat span { font-size: 13px; color: var(--ink2); }
  .summary p + p { margin-top: 8px; }
  .sbar { display: flex; gap: 2px; height: 22px; border-radius: 4px; overflow: hidden; background: var(--surface); }
  .sbar .seg { flex: 1 1 0; min-width: 0; height: 100%; }
  .sbar .seg:hover { filter: brightness(1.12); }
  .legend { display: flex; flex-wrap: wrap; gap: 6px 16px; font-size: 13px; color: var(--ink2); margin-top: 10px; align-items: center; }
  .sw { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 6px; vertical-align: -1px; }
  .legend b { color: var(--ink); font-weight: 600; margin-left: 4px; }
  .srows { display: grid; grid-template-columns: minmax(120px, max-content) minmax(0, 1fr); gap: 8px 14px; align-items: center; }
  .srows .sl { font-size: 13px; min-width: 0; }
  .srows .sl span { display: block; color: var(--muted); font-size: 12px; }
  .tablewrap { overflow-x: auto; border: 1px solid var(--border); border-radius: 6px; background: var(--surface); }
  table { border-collapse: collapse; width: 100%; font-size: 13.5px; }
  th, td { text-align: left; padding: 8px 12px; border-bottom: 1px solid var(--grid); vertical-align: top; }
  thead th { font-size: 12px; font-weight: 500; letter-spacing: .04em; text-transform: uppercase; color: var(--muted); background: var(--surface); position: sticky; top: 0; z-index: 1; }
  tbody tr:nth-child(even) td { background: var(--row-alt); }
  td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
  #ftable { min-width: 1050px; }
  #ftable td.c-f { width: 150px; } #ftable td.c-an { width: 30%; } #ftable td.c-alt { width: 24%; }
  #ftable tr.flash td { background: var(--accent-soft); transition: background .6s; }
  .fid { font-family: var(--mono); font-weight: 600; font-size: 16px; display: flex; align-items: center; gap: 7px; }
  .fsw { display: inline-block; width: 13px; height: 13px; border-radius: 3px; box-shadow: inset 0 0 0 1px var(--border); flex: none; }
  .morebtn { font-size: 12px; padding: 1px 8px; border-radius: 12px; }
  .tags { display: flex; flex-wrap: wrap; gap: 4px; margin-top: 4px; }
  .tag { font-size: 11.5px; border: 1px solid var(--border); border-radius: 3px; padding: 0 5px; color: var(--ink2); white-space: nowrap; }
  .tag.d-specific { border-color: var(--accent); color: var(--accent-ink); }
  .where { font-size: 12px; color: var(--muted); margin-top: 4px; word-break: break-word; }
  .pill { display: inline-block; font-size: 11.5px; border-radius: 3px; padding: 1px 7px; white-space: nowrap; }
  .p-assigned { background: var(--good-soft); color: var(--good); }
  .p-ambiguous { background: var(--warn-soft); color: var(--warn); }
  .p-unfitted, .p-pending { background: var(--none-soft); color: var(--none); }
  .lab { font-weight: 600; font-size: 14px; margin-top: 4px; word-break: break-word; }
  .sub { display: flex; flex-wrap: wrap; gap: 2px 10px; margin-top: 2px; font-size: 12.5px; color: var(--ink2); align-items: center; }
  .alias { font-family: var(--mono); }
  .why { margin-top: 6px; font-size: 12.5px; color: var(--ink2); }
  .altitem + .altitem { margin-top: 8px; padding-top: 8px; border-top: 1px dashed var(--grid); }
  .altlab { font-weight: 500; word-break: break-word; }
  .genes { display: flex; flex-wrap: wrap; gap: 2px 3px; margin-top: 3px; }
  .g { font-family: var(--mono); font-style: italic; font-size: 12.5px; padding: 0 3px; border-radius: 3px; }
  .g.s-exclusive { background: var(--accent-soft); color: var(--accent-ink); font-weight: 600; }
  .g.s-highest { box-shadow: inset 0 0 0 1px var(--accent); }
  .g.s-top3 { box-shadow: inset 0 0 0 1px var(--axis); }
  .g.s-shared { color: var(--muted); }
  .g.key { text-decoration: underline; text-decoration-thickness: 2px; text-underline-offset: 3px; }
  .persample { margin-top: 8px; display: grid; gap: 2px; font-size: 12px; }
  .sname { font-size: 11.5px; color: var(--muted); font-family: var(--mono); margin-right: 4px; }
  .linkbtn { margin-top: 8px; font-size: 12px; padding: 2px 8px; }
  .chip { display: inline-flex; align-items: center; gap: 5px; font-size: 12.5px; border: 1px solid var(--border); border-radius: 12px; padding: 0 8px 0 6px; margin: 2px 4px 2px 0; cursor: pointer; background: var(--surface); }
  .chip:hover { border-color: var(--axis); }
  .chip .fidm { font-family: var(--mono); color: var(--muted); }
  .dot { display: inline-block; width: 8px; height: 8px; border-radius: 50%; flex: none; }
  .controls { display: flex; flex-wrap: wrap; gap: 8px 14px; align-items: center; }
  .controls label { display: inline-flex; gap: 6px; align-items: center; font-size: 13px; color: var(--ink2); }
  .controls input[type=search] { flex: 1 1 220px; max-width: 320px; min-width: 0; }
  .controls input[type=number] { width: 64px; }
  .sticky { position: sticky; top: 0; z-index: 3; background: var(--bg); padding-block: 10px; border-bottom: 1px solid var(--border); }
  .count { margin-left: auto; font-size: 13px; color: var(--muted); font-variant-numeric: tabular-nums; }
  .note { font-size: 12.5px; color: var(--muted); }
  textarea { width: 100%; min-height: 54px; font-family: var(--mono); font-size: 12.5px; resize: vertical; }
  .dp-scroll { position: relative; overflow: auto; border: 1px solid var(--border); border-radius: 6px; background: var(--surface); outline-offset: -2px; }
  .dp-inner { position: relative; }
  .dp-canvas { position: sticky; top: 0; left: 0; display: block; }
  .dplegend { display: flex; flex-wrap: wrap; gap: 10px 32px; align-items: center; font-size: 12.5px; color: var(--ink2); }
  .dplegend svg { overflow: visible; }
  .cbar { width: 160px; height: 10px; border-radius: 2px; background: linear-gradient(90deg, var(--div-neg), var(--div-mid), var(--div-pos)); border: 1px solid var(--border); }
  .cticks { display: flex; justify-content: space-between; width: 160px; font-variant-numeric: tabular-nums; color: var(--muted); font-size: 11.5px; }
  .kv { display: grid; grid-template-columns: max-content minmax(0, 1fr); gap: 4px 16px; font-size: 13.5px; }
  .kv dt { color: var(--muted); } .kv dd { margin: 0; word-break: break-word; }
  pre { white-space: pre-wrap; word-break: break-word; font-family: var(--mono); font-size: 12px; background: var(--row-alt); padding: 12px; border-radius: 4px; max-height: 480px; overflow: auto; }
  details summary { cursor: pointer; color: var(--ink2); }
  .tip { position: fixed; z-index: 10; pointer-events: none; background: var(--surface); color: var(--ink); border: 1px solid var(--border); border-radius: 6px; padding: 8px 10px; font-size: 12.5px; box-shadow: 0 4px 16px rgba(0,0,0,.14); max-width: 340px; }
  .tip .tt { font-weight: 600; margin-bottom: 4px; word-break: break-word; }
  .tip .tip-row { display: flex; gap: 8px; align-items: baseline; }
  .tip .tip-row b { font-variant-numeric: tabular-nums; }
  .tip .tip-row span { color: var(--ink2); }
  footer { color: var(--muted); font-size: 12.5px; border-top: 1px solid var(--border); padding-top: 12px; }
  @media (max-width: 640px) { .page { padding-top: 16px; } .stat { min-width: 0; flex: 1 1 40%; } }
</style>
</head>
<body>
<div class="page">
  <header>
    <div class="hd">
      <div class="eyebrow" id="eyebrow"></div>
      <h1 id="title"></h1>
      <p class="meta" id="metaline"></p>
      <p class="meta" id="samplesline"></p>
    </div>
    <button id="theme" type="button" aria-label="Colour theme">Theme: auto</button>
  </header>
  <nav class="tabs" role="tablist" aria-label="Report sections">
    <button role="tab" data-tab="overview" aria-selected="true">Overview</button>
    <button role="tab" data-tab="factors" aria-selected="false">Factors</button>
    <button role="tab" data-tab="dotplot" aria-selected="false">Dotplot</button>
    <button role="tab" data-tab="samples" aria-selected="false">Samples</button>
    <button role="tab" data-tab="run" aria-selected="false">Run details</button>
  </nav>

  <section id="tab-overview" role="tabpanel">
    <div class="stats" id="stats"></div>
    <div class="card summary"><h2>Dataset summary</h2><div id="summary"></div>
      <p class="note" style="margin-top:10px">Written by the LLM from the factor evidence; labels are best guesses to check against the marker genes.</p></div>
    <div class="card"><h2>Composition by compartment</h2>
      <p class="note" style="margin-bottom:10px">Share of all counts assigned to the factors of each compartment (the compartment of each factor's best-guess label).</p>
      <div id="comp-bar"></div><div class="legend" id="comp-legend"></div></div>
    <div class="tablewrap"><table id="comp-table"><thead><tr><th>Compartment</th><th class="num">Factors</th><th class="num">Share of counts</th><th>Largest factors</th></tr></thead><tbody></tbody></table></div>
  </section>

  <section id="tab-factors" role="tabpanel" hidden>
    <div class="controls sticky" role="search">
      <input type="search" id="fq" placeholder="Search factor, label or gene" aria-label="Search factor, label or gene">
      <select id="fstatus" aria-label="Status"><option value="">Any status</option><option value="assigned">labelled</option><option value="ambiguous">low confidence</option><option value="unfitted">unresolved</option><option value="pending">not annotated</option></select>
      <select id="fkind" aria-label="Kind"><option value="">Cell types and programs</option><option value="cell_type">cell types</option><option value="program">programs</option></select>
      <select id="fcomp" aria-label="Compartment"><option value="">Any compartment</option></select>
      <select id="fdist" aria-label="Distribution"><option value="">Any distribution</option><option value="specific">specific (one sample)</option><option value="restricted">restricted</option><option value="shared">shared</option></select>
      <label><input type="checkbox" id="falt"> with alternatives</label>
      <select id="fsort" aria-label="Sort"><option value="factor">Sort: factor</option><option value="share">Sort: abundance</option><option value="confidence">Sort: confidence</option><option value="compartment">Sort: compartment</option></select>
      <span class="count" id="fcount" aria-live="polite"></span>
    </div>
    <div class="legend" id="gene-legend">
      <span>Marker genes, most characteristic first:</span>
      <span><span class="g s-exclusive">exclusive</span> ≥50% of the gene's expression</span>
      <span><span class="g s-highest">highest</span> this factor expresses it most</span>
      <span><span class="g s-top3">top3</span></span>
      <span><span class="g s-shared">shared</span> expressed more elsewhere</span>
      <span><span class="g key">underlined</span> key genes named by the LLM</span>
    </div>
    <div class="tablewrap"><table id="ftable"><thead><tr><th>Factor</th><th>Best-guess annotation</th><th>Alternative interpretations</th><th>Marker genes</th></tr></thead><tbody></tbody></table></div>
  </section>

  <section id="tab-dotplot" role="tabpanel" hidden>
    <div class="card" style="display:grid;gap:10px">
      <div class="controls">
        <label>Genes <select id="dgenes"><option value="key">key genes (LLM)</option><option value="top">top marker genes</option><option value="custom">custom list</option></select></label>
        <label id="dklab">per factor <input type="number" id="dk" min="1" max="20" value="2"></label>
        <label>Gene order <select id="dgorder"><option value="listed">as selected</option><option value="peak">by peak factor</option></select></label>
        <label>Factor order <select id="dorder"><option value="factor">factor id</option><option value="compartment">compartment</option><option value="share">abundance</option><option value="match" title="Factors that express the shown genes most strongly come first (mean of each gene's count relative to its max)">expression of shown genes</option></select></label>
      </div>
      <div id="dcustomwrap" hidden><textarea id="dcustom" placeholder="Gene symbols separated by spaces, commas or new lines" aria-label="Custom gene list"></textarea>
        <div class="controls" style="margin-top:6px"><button type="button" id="dapply">Apply gene list</button><span class="note" id="dcustomnote"></span></div></div>
      <div class="controls">
        <input type="search" id="dq" placeholder="Filter factors by id, label or alias" aria-label="Filter factors">
        <select id="dcomp" aria-label="Compartment"><option value="">All compartments</option></select>
        <label><input type="checkbox" id="dhideunres"> hide unresolved</label>
        <label>Dot size <select id="dsize"><option value="col">relative to each gene's max</option><option value="global">one scale for all genes</option></select></label>
        <label>Scale <select id="dscale"><option value="sqrt">area ∝ √count</option><option value="log">area ∝ log(1+count)</option><option value="area">area ∝ count</option></select></label>
        <label>Colour range ± <select id="dclamp"><option>2</option><option>3</option><option selected>4</option><option>6</option><option>8</option></select></label>
        <button type="button" id="ddownload">Download view (TSV)</button>
        <span class="count" id="dcount"></span>
      </div>
      <div class="dplegend" id="dlegend"></div>
    </div>
    <div id="dhost"></div>
    <p class="note">Rows are factors, columns genes. Dot size: the gene's count in the factor (factor model). Dot colour: log2 fold change of the gene in the factor against all other factors (the odds ratio of cartloader's bulk DE). Grey bars (log scale): each gene's total count over all factors, and each factor's total count. The square beside each factor id is its colour in CartoScope; the thin stripe is its compartment. A cell without a dot holds a count too small to draw (under a pixel at the current scale, e.g. below ~0.1 on the log scale); it is not necessarily zero, and hovering shows its value. Hover for values; click a factor to open it in the table; click a gene to order factors by it. "Show markers in dotplot" in the factor table plots that factor's markers across all factors, ordered by expression of those genes, with the factor highlighted. Arrow keys move through the plot when it has focus.</p>
  </section>

  <section id="tab-samples" role="tabpanel" hidden>
    <div class="card"><h2>Composition of each sample</h2>
      <p class="note" style="margin-bottom:12px">Share of each sample's decoded counts (pixel pseudobulk) by compartment.</p>
      <div class="srows" id="sample-bars"></div><div class="legend" id="sample-legend"></div></div>
    <div class="card" style="display:grid;gap:10px">
      <div class="controls">
        <input type="search" id="sq" placeholder="Filter factors by id, label or alias" aria-label="Filter factors">
        <select id="scomp" aria-label="Compartment"><option value="">All compartments</option></select>
        <label>Factor order <select id="sorder"><option value="factor">factor id</option><option value="compartment">compartment</option><option value="share">abundance</option><option value="peak">by peak sample</option></select></label>
        <label>Dot size <select id="ssize"><option value="row">relative to each factor's max</option><option value="col">relative to each sample's max</option><option value="global">one scale</option></select></label>
        <label>Scale <select id="sscale"><option value="sqrt">area ∝ √count</option><option value="log">area ∝ log(1+count)</option><option value="area">area ∝ count</option></select></label>
        <label>Colour range ± <select id="sclamp"><option>2</option><option selected>4</option><option>6</option><option>8</option></select></label>
      </div>
      <div class="dplegend" id="slegend"></div>
    </div>
    <div id="shost"></div>
    <p class="note">Rows are factors, columns samples. Dot size: the factor's count in the sample's pixel pseudobulk. Dot colour: log2 fold change of the factor in the sample against the other samples. Grey bars (log scale): each sample's decoded counts and each factor's decoded counts over all samples.</p>
  </section>

  <section id="tab-run" role="tabpanel" hidden>
    <div class="card"><h2>Run</h2><dl class="kv" id="runkv"></dl></div>
    <div class="tablewrap"><table id="calls"><thead><tr><th>Call</th><th class="num">Factors</th><th>Model</th><th class="num">Input tokens</th><th class="num">Output tokens</th><th>Response file</th></tr></thead><tbody></tbody></table></div>
    <div class="card" id="cmdcard"><h2>Command</h2><pre id="cmd"></pre></div>
  </section>

  <footer>Generated by cartloader annotate_factors_with_llm. The alias TSV next to this file holds the best-guess aliases in cartloader's index/alias format; the .annotations.tsv holds every column shown here.</footer>
</div>
<div class="tip" id="tip" hidden></div>
<script id="report-data" type="application/json">__DATA__</script>
<script>
(function () {
  "use strict";
  const D = JSON.parse(document.getElementById("report-data").textContent);
  const $ = (s) => document.querySelector(s);
  function el(tag, cls, text) { const e = document.createElement(tag); if (cls) e.className = cls; if (text != null) e.textContent = text; return e; }
  function cssVar(n) { return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }
  function fmtPct(x) { if (x == null || isNaN(x)) return "–"; const p = x * 100; return (p >= 10 ? p.toFixed(0) : p >= 1 ? p.toFixed(1) : p.toFixed(2)) + "%"; }
  function fmtNum(x) {
    if (x == null || isNaN(x)) return "–"; const a = Math.abs(x);
    return a >= 1e9 ? (x / 1e9).toFixed(2) + "G" : a >= 1e6 ? (x / 1e6).toFixed(2) + "M" : a >= 1e4 ? (x / 1e3).toFixed(1) + "k"
      : a >= 100 ? x.toFixed(0) : a >= 1 ? x.toFixed(1) : a === 0 ? "0" : x.toPrecision(2);
  }
  function fmtLfc(v) { return (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v).toFixed(2); }
  function decode(s, T) { const bin = atob(s || ""); const u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); return new T(u.buffer); }
  const store = { get(k) { try { return localStorage.getItem(k); } catch (e) { return null; } }, set(k, v) { try { localStorage.setItem(k, v); } catch (e) { /* no storage */ } } };

  const F = D.factors, COMP = D.compartments, nF = F.length;
  const FI = new Map(F.map((f, i) => [f.id, i]));
  const STATUS = { assigned: "labelled", ambiguous: "low confidence", unfitted: "unresolved", pending: "not annotated" };
  const KIND = { cell_type: "cell type", program: "program" };
  const compOf = f => (f.label && f.compartment) ? f.compartment : "unresolved";
  const compVar = c => { const i = COMP.indexOf(c); return i < 0 ? "var(--c-unres)" : "var(--c" + (i + 1) + ")"; };
  const compRank = c => { const i = COMP.indexOf(c); return i < 0 ? COMP.length : i; };
  const nameOf = f => f.label || (f.status === "pending" ? "not annotated" : "Unresolved");
  const shortOf = f => f.label ? f.alias : (f.status === "pending" ? "" : "Unresolved");
  const fcolor = f => f.color || null;
  const searchText = f => [f.id, f.label, f.alias, f.alternatives.map(a => a.label).join(" ")].join(" ").toLowerCase();
  const COMP_PRESENT = COMP.concat(["unresolved"]).filter(c => F.some(f => compOf(f) === c));

  // ---------- tooltip ----------
  const tip = $("#tip");
  function placeTip(x, y) {
    const r = tip.getBoundingClientRect(); let left = x + 14, top = y + 14;
    if (left + r.width > innerWidth - 8) left = x - r.width - 14;
    if (top + r.height > innerHeight - 8) top = y - r.height - 14;
    tip.style.left = Math.max(8, left) + "px"; tip.style.top = Math.max(8, top) + "px";
  }
  function showTip(x, y, build) { tip.replaceChildren(); build(tip); tip.hidden = false; placeTip(x, y); }
  function hideTip() { tip.hidden = true; }
  function tipTitle(t, text) { t.append(el("div", "tt", text)); }
  function tipRow(t, value, label) { const d = el("div", "tip-row"); d.append(el("b", null, value), el("span", null, label)); t.append(d); }

  // ---------- header ----------
  const M = D.meta;
  $("#eyebrow").textContent = "factor annotation · " + M.mode + " · " + M.date;
  $("#title").textContent = M.title;
  $("#metaline").textContent = nF + " factors · organism: " + M.organism + " · " +
    (D.samples.length && M.tissues.length > 1 ? M.tissues.length + " tissue contexts" : "tissue: " + M.tissues.join("; ")) +
    " · " + (M.dry_run ? "dry run, no LLM call" : "labels by " + M.model + " (effort " + M.effort + ")");
  if (D.samples.length) $("#samplesline").textContent = D.samples.length + " samples: " + D.samples.map(s => s.name + " (" + s.tissue + ")").join("; ");

  // ---------- tabs ----------
  const tabButtons = Array.from(document.querySelectorAll(".tabs [role=tab]"));
  if (!D.samples.length) tabButtons.find(b => b.dataset.tab === "samples").remove();
  const tabs = tabButtons.filter(b => b.isConnected);
  const onShow = {};
  function showTab(name) {
    if (!tabs.some(b => b.dataset.tab === name)) name = "overview";
    tabs.forEach(b => { const on = b.dataset.tab === name; b.setAttribute("aria-selected", on ? "true" : "false"); b.tabIndex = on ? 0 : -1; document.getElementById("tab-" + b.dataset.tab).hidden = !on; });
    hideTip();
    if (onShow[name]) onShow[name]();
    try { history.replaceState(null, "", "#" + name); } catch (e) { /* file:// may refuse */ }
  }
  tabs.forEach((b, i) => {
    b.addEventListener("click", () => showTab(b.dataset.tab));
    b.addEventListener("keydown", e => {
      if (e.key !== "ArrowRight" && e.key !== "ArrowLeft") return;
      const n = tabs[(i + (e.key === "ArrowRight" ? 1 : tabs.length - 1)) % tabs.length]; n.focus(); showTab(n.dataset.tab);
    });
  });

  // ---------- shared bits ----------
  function geneChip(g, cls) { return el("span", "g" + (cls ? " " + cls : ""), g); }
  function compTag(c) { const s = el("span"); const d = el("span", "dot"); d.style.background = compVar(c); s.append(d, document.createTextNode(" " + c)); s.style.display = "inline-flex"; s.style.alignItems = "center"; s.style.gap = "4px"; return s; }
  function factorChip(fi) {
    const f = F[fi], b = el("button", "chip"); b.type = "button";
    const d = el("span", "dot"); d.style.background = fcolor(f) || compVar(compOf(f));
    if (fcolor(f)) d.style.boxShadow = "inset 0 0 0 1px var(--border)";
    b.append(d, el("span", "fidm", f.id), document.createTextNode(shortOf(f) || nameOf(f)));
    b.title = nameOf(f); b.addEventListener("click", () => goFactor(f.id)); return b;
  }
  function fillCompSelect(sel) { COMP_PRESENT.forEach(c => { const o = el("option", null, c); o.value = c; sel.append(o); }); }
  function aggregate(valueOf) {
    const m = new Map();
    F.forEach((f, i) => { const c = compOf(f); const e = m.get(c) || { key: c, value: 0, n: 0, factors: [] }; e.value += valueOf(i); e.n++; e.factors.push(i); m.set(c, e); });
    return COMP_PRESENT.filter(c => m.has(c)).map(c => m.get(c));
  }
  function stackedBar(host, segs, what) {
    const total = segs.reduce((a, s) => a + s.value, 0) || 1;
    const bar = el("div", "sbar"); bar.setAttribute("role", "img");
    bar.setAttribute("aria-label", segs.map(s => s.key + " " + fmtPct(s.value / total)).join(", "));
    segs.forEach(s => {
      if (!(s.value > 0)) return;
      const d = el("div", "seg"); d.style.flexGrow = String(s.value / total); d.style.background = compVar(s.key);
      d.addEventListener("pointermove", e => showTip(e.clientX, e.clientY, t => { tipTitle(t, s.key); tipRow(t, fmtPct(s.value / total), what); if (s.n != null) tipRow(t, String(s.n), "factors"); }));
      d.addEventListener("pointerleave", hideTip); bar.append(d);
    });
    host.append(bar);
  }
  function compLegend(host, segs) {
    const total = segs.reduce((a, s) => a + s.value, 0) || 1;
    segs.forEach(s => { const it = el("span"); const sw = el("span", "sw"); sw.style.background = compVar(s.key); it.append(sw, document.createTextNode(s.key)); if (segs[0].value != null) it.append(el("b", null, fmtPct(s.value / total))); host.append(it); });
  }

  // ---------- overview ----------
  function renderOverview() {
    const n = k => F.filter(k).length;
    const tiles = [[n(f => f.status === "assigned"), "labelled"], [n(f => f.status === "ambiguous"), "low confidence"],
      [n(f => f.status === "unfitted"), "unresolved"], [n(f => f.label && f.kind === "cell_type"), "cell types"],
      [n(f => f.label && f.kind === "program"), "programs"], [n(f => f.alternatives.length > 0), "with alternative readings"]];
    const pend = n(f => f.status === "pending"); if (pend) tiles.unshift([pend, "not annotated (dry run)"]);
    tiles.forEach(([v, l]) => { const d = el("div", "stat"); d.append(el("b", null, String(v)), el("span", null, l)); $("#stats").append(d); });
    const sum = $("#summary");
    if (D.summaries.length) D.summaries.forEach(s => sum.append(el("p", null, s)));
    else sum.append(el("p", "muted", M.dry_run ? "Dry run: no LLM call was made, so there are no labels or summary yet. The evidence, factor table and dotplots are complete." : "The model returned no summary."));
    const hasShare = F.some(f => f.share != null);
    const segs = aggregate(i => hasShare ? (F[i].share || 0) : 1);
    stackedBar($("#comp-bar"), segs, hasShare ? "of all counts" : "of factors");
    compLegend($("#comp-legend"), segs);
    const tb = $("#comp-table tbody"), total = segs.reduce((a, s) => a + s.value, 0) || 1;
    segs.forEach(s => {
      const tr = el("tr"), c0 = el("td"); const sw = el("span", "sw"); sw.style.background = compVar(s.key); c0.append(sw, document.createTextNode(s.key));
      const cf = el("td"), sorted = s.factors.slice().sort((a, b) => (F[b].share || 0) - (F[a].share || 0));
      sorted.forEach((fi, j) => { const ch = factorChip(fi); if (j >= 8) ch.hidden = true; cf.append(ch); });
      if (sorted.length > 8) {
        const more = el("button", "chip morebtn", "+" + (sorted.length - 8) + " more"); more.type = "button"; more.setAttribute("aria-expanded", "false");
        more.addEventListener("click", () => {
          const open = more.getAttribute("aria-expanded") !== "true";
          Array.from(cf.querySelectorAll(".chip:not(.morebtn)")).forEach((ch, j) => { if (j >= 8) ch.hidden = !open; });
          more.setAttribute("aria-expanded", open ? "true" : "false"); more.textContent = open ? "show fewer" : "+" + (sorted.length - 8) + " more";
          cf.append(more);
        });
        cf.append(more);
      }
      tr.append(c0, el("td", "num", String(s.n)), el("td", "num", hasShare ? fmtPct(s.value / total) : "–"), cf); tb.append(tr);
    });
  }

  // ---------- factor table ----------
  const rowEl = new Map();
  function buildRow(f) {
    const tr = el("tr"); tr.id = "factor-" + f.id;
    const td1 = el("td", "c-f"), fid = el("div", "fid");
    if (fcolor(f)) { const sw = el("span", "fsw"); sw.style.background = fcolor(f); sw.title = "colour in CartoScope"; fid.append(sw); }
    fid.append(document.createTextNode(f.id)); td1.append(fid);
    const tags = el("div", "tags");
    if (f.share != null) { const t = el("span", "tag", fmtPct(f.share)); t.title = "share of all counts (" + f.abundance + ")"; tags.append(t); }
    if (f.de_support) { const t = el("span", "tag", f.de_support + " DE"); t.title = f.n_de + " significantly enriched genes"; tags.append(t); }
    if (f.distribution) tags.append(el("span", "tag d-" + f.distribution, f.distribution));
    td1.append(tags);
    if (f.found_in.length) td1.append(el("div", "where", f.found_in.join(", ")));
    const td2 = el("td", "c-an"); td2.append(el("span", "pill p-" + f.status, STATUS[f.status] || f.status));
    const lab = el("div", "lab" + (f.label ? "" : " muted"), nameOf(f)); td2.append(lab);
    if (f.label) {
      const sub = el("div", "sub"); sub.append(el("span", "alias", f.alias), el("span", null, KIND[f.kind] || f.kind), compTag(compOf(f)));
      if (f.confidence) sub.append(el("span", null, f.confidence + " confidence")); td2.append(sub);
    }
    if (f.rationale) td2.append(el("div", "why", f.rationale));
    const td3 = el("td", "c-alt");
    if (f.alternatives.length) {
      if (!f.label) td3.append(el("div", "note", "possible readings:"));
      f.alternatives.forEach(a => {
        const d = el("div", "altitem"); d.append(el("div", "altlab", a.label));
        const m = el("div", "sub"); m.append(el("span", null, KIND[a.kind] || a.kind), el("span", null, a.confidence + " confidence")); d.append(m);
        if (a.genes.length) { const gs = el("div", "genes"); a.genes.forEach(g => gs.append(geneChip(g))); d.append(gs); }
        td3.append(d);
      });
    } else td3.append(el("span", "muted", f.status === "pending" ? "–" : "none"));
    const td4 = el("td", "c-mk"), gs = el("div", "genes"), key = new Set(f.key_genes);
    f.markers.forEach(g => gs.append(geneChip(g, "s-" + (f.spec[g] || "none") + (key.has(g) ? " key" : ""))));
    if (!f.markers.length) gs.append(el("span", "muted", "no marker genes"));
    td4.append(gs);
    const mbs = Object.entries(f.markers_by_sample || {});
    if (mbs.length) { const ps = el("div", "persample"); mbs.forEach(([s, g]) => { const d = el("div"); d.append(el("span", "sname", s)); g.forEach(x => d.append(geneChip(x), document.createTextNode(" "))); ps.append(d); }); td4.append(ps); }
    if (D.dot && f.markers.length) { const b = el("button", "linkbtn", "Show markers in dotplot"); b.type = "button"; b.addEventListener("click", () => dotForFactor(f)); td4.append(b); }
    tr.append(td1, td2, td3, td4); return tr;
  }
  const fctl = { q: $("#fq"), status: $("#fstatus"), kind: $("#fkind"), comp: $("#fcomp"), dist: $("#fdist"), alt: $("#falt"), sort: $("#fsort") };
  function renderFactors() {
    fillCompSelect(fctl.comp);
    if (!D.samples.length) fctl.dist.remove();
    if (!F.some(f => f.status === "pending")) fctl.status.querySelector("option[value=pending]").remove();
    const tb = $("#ftable tbody"); F.forEach(f => { const r = buildRow(f); rowEl.set(f.id, r); tb.append(r); });
    Object.values(fctl).forEach(c => { if (c.isConnected) { c.addEventListener("input", applyFactorFilter); c.addEventListener("change", applyFactorFilter); } });
    applyFactorFilter();
  }
  const CONF = { high: 0, medium: 1, low: 2 };
  function applyFactorFilter() {
    const q = fctl.q.value.trim().toLowerCase(), st = fctl.status.value, kd = fctl.kind.value, cp = fctl.comp.value;
    const ds = fctl.dist.isConnected ? fctl.dist.value : "", alt = fctl.alt.checked;
    let shown = 0;
    const order = F.map((f, i) => i);
    const s = fctl.sort.value;
    if (s === "share") order.sort((a, b) => (F[b].share || 0) - (F[a].share || 0));
    else if (s === "confidence") order.sort((a, b) => ((F[a].label ? CONF[F[a].confidence] ?? 3 : 4) - (F[b].label ? CONF[F[b].confidence] ?? 3 : 4)) || ((F[b].share || 0) - (F[a].share || 0)));
    else if (s === "compartment") order.sort((a, b) => (compRank(compOf(F[a])) - compRank(compOf(F[b]))) || (a - b));
    const tb = $("#ftable tbody");
    order.forEach(i => {
      const f = F[i], r = rowEl.get(f.id);
      const ok = (!st || f.status === st) && (!kd || (f.label && f.kind === kd)) && (!cp || compOf(f) === cp) && (!ds || f.distribution === ds)
        && (!alt || f.alternatives.length) && (!q || (searchText(f) + " " + f.markers.join(" ").toLowerCase()).indexOf(q) >= 0);
      r.hidden = !ok; if (ok) shown++; tb.append(r);
    });
    $("#fcount").textContent = "Showing " + shown + " of " + nF + " factors";
  }
  function goFactor(id) {
    showTab("factors");
    const r = rowEl.get(id); if (!r) return;
    if (r.hidden) { fctl.q.value = ""; fctl.status.value = ""; fctl.kind.value = ""; fctl.comp.value = ""; if (fctl.dist.isConnected) fctl.dist.value = ""; fctl.alt.checked = false; applyFactorFilter(); }
    r.scrollIntoView({ block: "center" }); r.classList.add("flash"); setTimeout(() => r.classList.remove("flash"), 1400);
  }

  // ---------- dotplot component ----------
  function Dotplot(host, opt) {
    const scroll = el("div", "dp-scroll"), inner = el("div", "dp-inner"), cv = el("canvas", "dp-canvas");
    scroll.tabIndex = 0; scroll.setAttribute("role", "img"); scroll.setAttribute("aria-label", opt.ariaLabel || "dotplot");
    inner.append(cv); scroll.append(inner); host.append(scroll);
    const ctx = cv.getContext("2d"), cell = opt.cell || 16, PAD = 10, SW = 22;
    const barW = opt.rowTotal ? 64 : 0, barH = opt.colTotal ? 46 : 0;
    // labelW/labelH are the plot origin: row label text + row total bars, column total bars + column label text
    let nr = 0, nc = 0, textW = 120, textH = 80, labelW = 120, labelH = 80, vw = 0, vh = 0, colMax = null, rowMax = null, gMax = 0;
    let hover = null, cursor = null, lut = null, raf = 0, rowScale = null, colScale = null;
    function font(px, mono, italic, weight) { return (italic ? "italic " : "") + (weight || 400) + " " + px + "px " + (mono ? cssVar("--mono") : cssVar("--sans")); }
    function measure() {
      ctx.font = font(12, false); let w = 0;
      for (let r = 0; r < nr; r++) { const [a, b] = opt.rowLabel(r); ctx.font = font(11.5, true); let x = ctx.measureText(a).width + 8; ctx.font = font(12, false); x += ctx.measureText(b).width; if (x > w) w = x; }
      textW = Math.min(340, Math.ceil(w) + SW + 8); labelW = textW + barW;
      ctx.font = font(12, !!opt.colMono, !!opt.colItalic); let h = 0;
      for (let c = 0; c < nc; c++) { const x = ctx.measureText(opt.colLabel(c)).width; if (x > h) h = x; }
      textH = Math.min(170, Math.ceil(h) + 16); labelH = barH + textH;
    }
    function logScale(n, get, maxTicks) {
      let lo = Infinity, hi = 0;
      for (let i = 0; i < n; i++) { const v = get(i); if (v > 0) { if (v < lo) lo = v; if (v > hi) hi = v; } }
      if (!(hi > 0)) return null;
      const a = Math.floor(Math.log10(lo)), b = Math.max(a + 1, Math.log10(hi));
      const top = Math.floor(b), ticks = [];
      if (maxTicks <= 2) { ticks.push(a); if (top > a) ticks.push(top); }
      else { const step = Math.ceil((top - a + 1) / maxTicks); for (let k = a; k <= b; k += step) ticks.push(k); }
      return { frac: v => v > 0 ? Math.max(0, Math.min(1, (Math.log10(v) - a) / (b - a))) : 0, ticks, tickFrac: k => (k - a) / (b - a) };
    }
    function fmtPow(k) { const n = ["1", "10", "100", "1k", "10k", "100k", "1M", "10M", "100M", "1G", "10G"]; return k >= 0 && k < n.length ? n[k] : "1e" + k; }
    function computeNorm() {
      colMax = new Float64Array(nc); rowMax = new Float64Array(nr); gMax = 0;
      for (let r = 0; r < nr; r++) for (let c = 0; c < nc; c++) { const v = opt.count(r, c); if (v > colMax[c]) colMax[c] = v; if (v > rowMax[r]) rowMax[r] = v; if (v > gMax) gMax = v; }
      rowScale = opt.rowTotal ? logScale(nr, opt.rowTotal, 2) : null; colScale = opt.colTotal ? logScale(nc, opt.colTotal, 3) : null;
    }
    function sizeOf(v, r, c) {
      const mode = opt.sizeMode(), base = mode === "col" ? colMax[c] : mode === "row" ? rowMax[r] : gMax;
      if (!(v > 0) || !(base > 0)) return 0;
      const sc = opt.scale ? opt.scale() : "area";
      return sc === "log" ? Math.log1p(v) / Math.log1p(base) : sc === "sqrt" ? Math.sqrt(v / base) : v / base;
    }
    function hexRgb(h) { h = h.replace("#", ""); if (h.length === 3) h = h.split("").map(x => x + x).join(""); const n = parseInt(h, 16); return [n >> 16 & 255, n >> 8 & 255, n & 255]; }
    function buildLut() {
      const neg = hexRgb(cssVar("--div-neg")), mid = hexRgb(cssVar("--div-mid")), pos = hexRgb(cssVar("--div-pos")); lut = [];
      for (let i = -64; i <= 64; i++) { const t = Math.abs(i) / 64, p = i < 0 ? neg : pos; lut.push("rgb(" + [0, 1, 2].map(k => Math.round(mid[k] + (p[k] - mid[k]) * t)).join(",") + ")"); }
    }
    function colorOf(l) { const L = opt.clamp(); const t = Math.max(-1, Math.min(1, l / L)); return lut[Math.round(t * 64) + 64]; }
    function layout() {
      const W = labelW + nc * cell + PAD, H = labelH + nr * cell + PAD;
      inner.style.width = W + "px"; inner.style.height = H + "px";
      scroll.style.height = Math.max(160, Math.min(H + 2, Math.round(innerHeight * 0.74))) + "px";
      resize();
    }
    function resize() {
      const W = labelW + nc * cell + PAD, H = labelH + nr * cell + PAD;
      vw = Math.min(scroll.clientWidth, W); vh = Math.min(scroll.clientHeight, H);
      const dpr = window.devicePixelRatio || 1;
      cv.width = Math.max(1, Math.round(vw * dpr)); cv.height = Math.max(1, Math.round(vh * dpr));
      cv.style.width = vw + "px"; cv.style.height = vh + "px"; ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      draw();
    }
    function draw() {
      raf = 0; if (!vw || !vh) return;
      if (!lut) buildLut();
      const sx = scroll.scrollLeft, sy = scroll.scrollTop;
      const surface = cssVar("--surface"), ink = cssVar("--ink"), ink2 = cssVar("--ink2"), muted = cssVar("--muted");
      const rowAlt = cssVar("--row-alt"), hov = cssVar("--hover"), axis = cssVar("--axis"), ring = cssVar("--border");
      ctx.fillStyle = surface; ctx.fillRect(0, 0, vw, vh);
      const c0 = Math.max(0, Math.floor(sx / cell)), c1 = Math.min(nc, Math.ceil((sx + vw - labelW) / cell) + 1);
      const r0 = Math.max(0, Math.floor(sy / cell)), r1 = Math.min(nr, Math.ceil((sy + vh - labelH) / cell) + 1);
      const X = c => labelW + c * cell - sx, Y = r => labelH + r * cell - sy;
      const focus = cursor || hover;
      ctx.save(); ctx.beginPath(); ctx.rect(labelW, labelH, vw - labelW, vh - labelH); ctx.clip();
      for (let r = r0; r < r1; r++) if (r % 2) { ctx.fillStyle = rowAlt; ctx.fillRect(labelW, Y(r), vw, cell); }
      if (opt.markRow != null && opt.markRow() != null) { const r = opt.markRow(); ctx.fillStyle = hov; ctx.fillRect(labelW, Y(r), vw, cell); }
      if (focus && focus.r != null) { ctx.fillStyle = hov; ctx.fillRect(labelW, Y(focus.r), vw, cell); }
      if (focus && focus.c != null) { ctx.fillStyle = hov; ctx.fillRect(X(focus.c), labelH, cell, vh); }
      const maxR = cell * 0.46; ctx.lineWidth = 0.75; ctx.strokeStyle = ring;
      for (let r = r0; r < r1; r++) for (let c = c0; c < c1; c++) {
        const v = opt.count(r, c), s = sizeOf(v, r, c); if (s <= 0) continue;
        const rad = maxR * Math.sqrt(s); if (rad < 0.6) continue;
        ctx.beginPath(); ctx.arc(X(c) + cell / 2, Y(r) + cell / 2, rad, 0, 6.2832);
        ctx.fillStyle = colorOf(opt.lfc(r, c)); ctx.fill(); if (rad > 2) ctx.stroke();
      }
      ctx.restore();
      // row labels
      const barCol = cssVar("--bar"), barOn = cssVar("--bar-on");
      ctx.save(); ctx.fillStyle = surface; ctx.fillRect(0, labelH, labelW, vh); ctx.beginPath(); ctx.rect(0, labelH, textW - 4, vh - labelH); ctx.clip();
      ctx.textBaseline = "middle";
      for (let r = r0; r < r1; r++) {
        const [a, b, sw, fc] = opt.rowLabel(r), y = Y(r) + cell / 2, on = focus && focus.r === r;
        if (sw) { ctx.fillStyle = sw; ctx.fillRect(3, y - 5, 3, 10); }
        if (fc) { ctx.fillStyle = fc; ctx.fillRect(9, y - 5, 10, 10); ctx.strokeStyle = ring; ctx.lineWidth = 1; ctx.strokeRect(9.5, y - 4.5, 9, 9); }
        ctx.font = font(11.5, true); ctx.fillStyle = muted; ctx.fillText(a, SW, y);
        const aw = ctx.measureText(a).width; ctx.font = font(12, false, false, on ? 650 : 400); ctx.fillStyle = on ? ink : ink2; ctx.fillText(b, SW + aw + 7, y);
      }
      ctx.restore();
      if (rowScale) {  // row totals: horizontal log-scale bars between the labels and the plot
        const x0 = textW + 2, len = barW - 10;
        ctx.save(); ctx.beginPath(); ctx.rect(textW, labelH, barW, vh - labelH); ctx.clip();
        ctx.strokeStyle = cssVar("--grid"); ctx.lineWidth = 1; ctx.beginPath();
        rowScale.ticks.forEach(k => { const x = Math.round(x0 + len * rowScale.tickFrac(k)) + 0.5; ctx.moveTo(x, labelH); ctx.lineTo(x, vh); }); ctx.stroke();
        for (let r = r0; r < r1; r++) {
          const w = len * rowScale.frac(opt.rowTotal(r)); if (w <= 0) continue;
          ctx.fillStyle = focus && focus.r === r ? barOn : barCol; ctx.fillRect(x0, Y(r) + cell * 0.2, Math.max(1, w), cell * 0.6);
        }
        ctx.restore();
      }
      // column labels
      ctx.save(); ctx.fillStyle = surface; ctx.fillRect(labelW, 0, vw, labelH); ctx.beginPath(); ctx.rect(labelW, barH, vw - labelW, textH); ctx.clip();
      ctx.textBaseline = "middle";
      for (let c = c0; c < c1; c++) {
        const on = focus && focus.c === c;
        ctx.save(); ctx.translate(X(c) + cell / 2, labelH - 8); ctx.rotate(-Math.PI / 2);
        ctx.font = font(12, !!opt.colMono, !!opt.colItalic, on ? 650 : 400); ctx.fillStyle = on ? ink : ink2; ctx.fillText(opt.colLabel(c), 0, 0); ctx.restore();
      }
      ctx.restore();
      if (colScale) {  // column totals: vertical log-scale bars above the column labels
        const y0 = barH - 4, len = barH - 8;
        ctx.save(); ctx.beginPath(); ctx.rect(labelW, 0, vw - labelW, barH); ctx.clip();
        ctx.strokeStyle = cssVar("--grid"); ctx.lineWidth = 1; ctx.beginPath();
        colScale.ticks.forEach(k => { const y = Math.round(y0 - len * colScale.tickFrac(k)) + 0.5; ctx.moveTo(labelW, y); ctx.lineTo(vw, y); }); ctx.stroke();
        for (let c = c0; c < c1; c++) {
          const h = len * colScale.frac(opt.colTotal(c)); if (h <= 0) continue;
          ctx.fillStyle = focus && focus.c === c ? barOn : barCol; ctx.fillRect(X(c) + cell * 0.2, y0 - h, cell * 0.6, Math.max(1, h));
        }
        ctx.restore();
      }
      ctx.fillStyle = surface; ctx.fillRect(0, 0, labelW, labelH);
      ctx.fillStyle = muted; ctx.font = font(11.5, false); ctx.textBaseline = "alphabetic";
      ctx.fillText((opt.cornerRows || "factors") + " ↓", 8, labelH - 8); ctx.fillText((opt.cornerCols || "genes") + " →", 8, labelH - 24);
      ctx.font = font(10.5, false);
      if (colScale) {
        ctx.textAlign = "right"; ctx.textBaseline = "middle";
        colScale.ticks.forEach(k => ctx.fillText(fmtPow(k), labelW - 6, barH - 4 - (barH - 8) * colScale.tickFrac(k)));
        ctx.textAlign = "left"; ctx.textBaseline = "top"; ctx.fillText("total count (log)", 8, 4);
      }
      if (rowScale) {
        ctx.textBaseline = "alphabetic";
        rowScale.ticks.forEach((k, i) => { ctx.textAlign = i === 0 ? "left" : "right"; ctx.fillText(fmtPow(k), textW + 2 + (barW - 10) * rowScale.tickFrac(k) + (i === 0 ? -1 : 3), labelH - 5); });
      }
      ctx.textAlign = "left";
      ctx.strokeStyle = axis; ctx.lineWidth = 1; ctx.beginPath();
      ctx.moveTo(labelW - 0.5, labelH); ctx.lineTo(labelW - 0.5, vh); ctx.moveTo(labelW, labelH - 0.5); ctx.lineTo(vw, labelH - 0.5); ctx.stroke();
    }
    function schedule() { if (!raf) raf = requestAnimationFrame(draw); }
    function locate(e) {
      const b = cv.getBoundingClientRect(), x = e.clientX - b.left, y = e.clientY - b.top;
      if (x < 0 || y < 0 || x > vw || y > vh) return null;
      const c = Math.floor((x - labelW + scroll.scrollLeft) / cell), r = Math.floor((y - labelH + scroll.scrollTop) / cell);
      const inR = r >= 0 && r < nr && y >= labelH, inC = c >= 0 && c < nc && x >= labelW;
      if (inR && inC) return { r, c }; if (inR && x < labelW) return { r, c: null }; if (inC && y < labelH) return { r: null, c };
      return null;
    }
    function tipFor(h, x, y) {
      if (!h) { hideTip(); return; }
      showTip(x, y, t => { if (h.r != null && h.c != null) opt.cellTip(t, h.r, h.c); else if (h.r != null) opt.rowTip(t, h.r); else opt.colTip(t, h.c); });
    }
    scroll.addEventListener("scroll", () => { hideTip(); schedule(); }, { passive: true });
    scroll.addEventListener("pointermove", e => { const h = locate(e); cursor = null; const same = h && hover && h.r === hover.r && h.c === hover.c; hover = h; if (!same) schedule(); tipFor(h, e.clientX, e.clientY); scroll.style.cursor = h ? "pointer" : "default"; });
    scroll.addEventListener("pointerleave", () => { hover = null; hideTip(); schedule(); });
    scroll.addEventListener("click", e => { const h = locate(e); if (!h) return; if (h.r != null && opt.onRow) opt.onRow(h.r); else if (h.c != null && opt.onCol) opt.onCol(h.c); });
    scroll.addEventListener("keydown", e => {
      const k = e.key; if (!nr || !nc) return;
      if (k === "Escape") { cursor = null; hideTip(); schedule(); return; }
      if (k === "Enter" && cursor && opt.onRow) { opt.onRow(cursor.r); return; }
      const d = { ArrowUp: [-1, 0], ArrowDown: [1, 0], ArrowLeft: [0, -1], ArrowRight: [0, 1] }[k]; if (!d) return;
      e.preventDefault(); cursor = cursor || { r: 0, c: 0 };
      cursor = { r: Math.max(0, Math.min(nr - 1, cursor.r + d[0])), c: Math.max(0, Math.min(nc - 1, cursor.c + d[1])) };
      const cx = labelW + cursor.c * cell, cy = labelH + cursor.r * cell;
      if (cx - scroll.scrollLeft < labelW) scroll.scrollLeft = cx - labelW; else if (cx + cell - scroll.scrollLeft > vw) scroll.scrollLeft = cx + cell - vw;
      if (cy - scroll.scrollTop < labelH) scroll.scrollTop = cy - labelH; else if (cy + cell - scroll.scrollTop > vh) scroll.scrollTop = cy + cell - vh;
      draw(); const b = cv.getBoundingClientRect();
      tipFor(cursor, b.left + labelW + cursor.c * cell - scroll.scrollLeft + cell, b.top + labelH + cursor.r * cell - scroll.scrollTop + cell);
    });
    new ResizeObserver(() => { if (scroll.offsetParent) resize(); }).observe(scroll);
    return {
      set(rows, cols) { nr = rows; nc = cols; cursor = null; hover = null; computeNorm(); measure(); layout(); scroll.scrollTop = 0; scroll.scrollLeft = 0; draw(); },
      redraw() { lut = null; computeNorm(); draw(); },
      max() { return gMax; },
    };
  }
  function sizeLegend(host, labels) {
    const NS = "http://www.w3.org/2000/svg", svg = document.createElementNS(NS, "svg"), cell = 16, maxR = cell * 0.46;
    let x = 0; const H = 20; svg.setAttribute("height", H); svg.setAttribute("aria-hidden", "true");
    labels.forEach(([s, text]) => {
      const r = maxR * Math.sqrt(s), c = document.createElementNS(NS, "circle");
      c.setAttribute("cx", x + maxR); c.setAttribute("cy", H / 2); c.setAttribute("r", r); c.setAttribute("fill", "var(--ink2)"); svg.append(c);
      const t = document.createElementNS(NS, "text"); t.setAttribute("x", x + 2 * maxR + 5); t.setAttribute("y", H / 2 + 4); t.setAttribute("font-size", "11.5"); t.setAttribute("fill", "var(--ink2)"); t.textContent = text; svg.append(t);
      x += 2 * maxR + 12 + text.length * 6.6 + 10;
    });
    svg.setAttribute("width", x); host.append(svg);
  }
  // Size legend: circles at area fractions 1, 1/2, 1/4, labelled in counts (one scale) or relative to the max.
  function sizeKey(host, scale, relativeTo, max) {
    const what = { log: "log(1+count)", sqrt: "√count", area: "count" }[scale];
    const box = el("div"); box.append(el("div", null, "Dot area ∝ " + what + (relativeTo ? ", relative to " + relativeTo : "")));
    const inv = f => scale === "log" ? Math.expm1(f * Math.log1p(max)) : scale === "sqrt" ? f * f * max : f * max;
    const lab = f => relativeTo ? (f === 1 ? "max" : scale === "log" ? "max^" + f : Math.round(100 * inv(f) / max) + "% of max") : fmtNum(inv(f));
    sizeLegend(box, [[1, lab(1)], [0.5, lab(0.5)], [0.25, lab(0.25)]]); host.append(box);
  }
  function colorLegend(host, L, caption) {
    const w = el("div"); w.append(el("div", null, caption));
    const bar = el("div", "cbar"), ticks = el("div", "cticks"); ticks.append(el("span", null, "≤ −" + L), el("span", null, "0"), el("span", null, "≥ +" + L));
    w.append(bar, ticks); host.append(w);
  }
  function download(name, text) {
    const a = el("a"); a.href = URL.createObjectURL(new Blob([text], { type: "text/tab-separated-values" })); a.download = name;
    document.body.append(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500);
  }
  function rowSwatch(f) { const v = compVar(compOf(f)); return cssVar(v.slice(4, -1)); }

  // ---------- gene dotplot ----------
  let gdot = null; const ds = { rows: [], cols: [], mark: null, sortCol: null };
  function initGeneDot() {
    if (!D.dot || !D.dot.genes.length) { $("#dhost").append(el("p", "note", "No factor model matrix was available, so there is no dotplot.")); return; }
    const G = D.dot.genes, nG = G.length, CNT = decode(D.dot.counts, Float32Array), LFC = decode(D.dot.lfc, Int16Array);
    const GI = new Map(G.map((g, i) => [g, i])), GL = new Map(G.map((g, i) => [g.toLowerCase(), i]));
    const GT = new Float64Array(nG); for (let g = 0; g < nG; g++) { let s = 0; for (let f = 0; f < nF; f++) s += CNT[g * nF + f]; GT[g] = s; }
    const c = { genes: $("#dgenes"), k: $("#dk"), gorder: $("#dgorder"), order: $("#dorder"), custom: $("#dcustom"), q: $("#dq"), comp: $("#dcomp"), hide: $("#dhideunres"), size: $("#dsize"), scale: $("#dscale"), clamp: $("#dclamp") };
    fillCompSelect(c.comp);
    const cnt = (r, col) => CNT[ds.cols[col] * nF + ds.rows[r]], lfc = (r, col) => LFC[ds.cols[col] * nF + ds.rows[r]] / 100;
    gdot = Dotplot($("#dhost"), {
      ariaLabel: "Dotplot of factors by genes", colItalic: true, colMono: true,
      rowLabel: r => { const f = F[ds.rows[r]]; return [f.id, shortOf(f) || "", rowSwatch(f), fcolor(f)]; },
      colLabel: col => G[ds.cols[col]], count: cnt, lfc: lfc,
      rowTotal: F.some(f => f.total != null) ? r => F[ds.rows[r]].total || 0 : null, colTotal: col => GT[ds.cols[col]],
      sizeMode: () => c.size.value, scale: () => c.scale.value, clamp: () => +c.clamp.value,
      markRow: () => ds.mark == null ? null : (ds.rows.indexOf(ds.mark) >= 0 ? ds.rows.indexOf(ds.mark) : null),
      cellTip: (t, r, col) => {
        const f = F[ds.rows[r]], g = ds.cols[col], v = cnt(r, col);
        tipTitle(t, G[g] + " in factor " + f.id); if (f.label) t.append(el("div", "note", f.label));
        tipRow(t, fmtNum(v), "count in factor"); tipRow(t, fmtLfc(lfc(r, col)), "log2 fold change vs other factors");
        tipRow(t, fmtPct(GT[g] > 0 ? v / GT[g] : 0), "of the gene's counts");
      },
      rowTip: (t, r) => { const f = F[ds.rows[r]]; tipTitle(t, "Factor " + f.id + " · " + nameOf(f)); if (f.label) tipRow(t, KIND[f.kind] || "", compOf(f)); if (f.total != null) tipRow(t, fmtNum(f.total), "total count"); if (f.share != null) tipRow(t, fmtPct(f.share), "of all counts"); t.append(el("div", "note", "click to open in the factor table")); },
      colTip: (t, col) => { const g = ds.cols[col]; tipTitle(t, G[g]); tipRow(t, fmtNum(GT[g]), "total count"); t.append(el("div", "note", ds.sortCol === g ? "click to restore factor order" : "click to order factors by this gene")); },
      onRow: r => goFactor(F[ds.rows[r]].id),
      onCol: col => { const g = ds.cols[col]; ds.sortCol = ds.sortCol === g ? null : g; update(); },
    });
    function rowsNow() {
      const q = c.q.value.trim().toLowerCase(), cp = c.comp.value;
      let rows = F.map((f, i) => i).filter(i => { const f = F[i]; return (!cp || compOf(f) === cp) && (!c.hide.checked || f.label) && (!q || searchText(f).indexOf(q) >= 0); });
      const o = c.order.value;
      if (o === "compartment") rows.sort((a, b) => (compRank(compOf(F[a])) - compRank(compOf(F[b]))) || (a - b));
      else if (o === "share") rows.sort((a, b) => (F[b].share || 0) - (F[a].share || 0));
      return rows;
    }
    function genesFor(rows) {
      const mode = c.genes.value, k = Math.max(1, +c.k.value || 2), out = [], seen = new Set(), unknown = [];
      const add = g => { const i = GI.has(g) ? GI.get(g) : GL.get(String(g).toLowerCase()); if (i == null) { unknown.push(g); return; } if (!seen.has(i)) { seen.add(i); out.push(i); } };
      if (mode === "custom") c.custom.value.split(/[\s,;]+/).filter(Boolean).forEach(add);
      else rows.forEach(fi => { const f = F[fi]; const src = mode === "key" && f.key_genes.length ? f.key_genes : f.markers; src.slice(0, k).forEach(add); });
      $("#dcustomnote").textContent = mode === "custom" && unknown.length ? "Not in the plot's gene set: " + unknown.slice(0, 20).join(", ") + (unknown.length > 20 ? " …" : "") : "";
      return out;
    }
    function update() {
      $("#dcustomwrap").hidden = c.genes.value !== "custom"; $("#dklab").hidden = c.genes.value === "custom";
      let rows = rowsNow(); let cols = genesFor(rows);
      const colMax = cols.map(g => { let m = 0; rows.forEach(fi => { const v = CNT[g * nF + fi]; if (v > m) m = v; }); return m || 1; });
      if (c.order.value === "match" && cols.length) {
        const score = fi => cols.reduce((a, g, j) => a + CNT[g * nF + fi] / colMax[j], 0);
        const sc = new Map(rows.map(fi => [fi, score(fi)])); rows.sort((a, b) => sc.get(b) - sc.get(a));
      }
      if (ds.sortCol != null) { const g = ds.sortCol; rows.sort((a, b) => CNT[g * nF + b] - CNT[g * nF + a]); }
      if (c.gorder.value === "peak" && cols.length) {
        const pos = new Map(rows.map((fi, i) => [fi, i]));
        const peak = cols.map(g => { let best = -1, bi = 0; rows.forEach(fi => { const v = CNT[g * nF + fi]; if (v > best) { best = v; bi = pos.get(fi); } }); return bi; });
        cols = cols.map((g, j) => [g, peak[j]]).sort((a, b) => a[1] - b[1]).map(x => x[0]);
      }
      ds.rows = rows; ds.cols = cols;
      $("#dcount").textContent = rows.length + " factors × " + cols.length + " genes" + (ds.sortCol != null ? " · ordered by " + G[ds.sortCol] : "");
      gdot.set(rows.length, cols.length); legend();
    }
    function legend() {
      const host = $("#dlegend"); host.replaceChildren();
      sizeKey(host, c.scale.value, c.size.value === "col" ? "the gene's max among shown factors" : null, gdot.max());
      colorLegend(host, c.clamp.value, "Colour: log2 fold change vs other factors");
    }
    ["genes", "k", "gorder", "order", "comp", "hide", "size", "scale", "clamp"].forEach(k => c[k].addEventListener("change", () => { if (k === "genes" || k === "order") ds.sortCol = null; update(); }));
    c.q.addEventListener("input", update);
    $("#dapply").addEventListener("click", update);
    $("#ddownload").addEventListener("click", () => {
      const lines = ["factor\tlabel\tgene\tcount\tlog2_fold_change"];
      ds.rows.forEach(fi => ds.cols.forEach(g => lines.push([F[fi].id, F[fi].label, G[g], CNT[g * nF + fi].toPrecision(6), (LFC[g * nF + fi] / 100).toFixed(2)].join("\t"))));
      download((M.stem || "dotplot") + ".dotplot.tsv", lines.join("\n") + "\n");
    });
    initGeneDot.update = update; initGeneDot.ctl = c;
    update();
  }
  function dotForFactor(f) {
    ensureGeneDot(); const c = initGeneDot.ctl; if (!c) return;
    c.genes.value = "custom"; c.custom.value = f.markers.join(" "); c.order.value = "match"; c.gorder.value = "listed";
    c.q.value = ""; c.comp.value = ""; ds.sortCol = null; ds.mark = FI.get(f.id);
    showTab("dotplot"); initGeneDot.update(); window.scrollTo({ top: 0 });
  }

  let dotReady = false;
  function ensureGeneDot() { if (!dotReady) { dotReady = true; initGeneDot(); } }

  // ---------- samples ----------
  let sdot = null;
  function initSamples() {
    if (!D.samples.length || !D.sdot) return;
    const S = D.samples, nS = S.length, CNT = decode(D.sdot.counts, Float32Array), LFC = decode(D.sdot.lfc, Int16Array);
    const ST = new Float64Array(nS), FT = new Float64Array(nF);
    for (let f = 0; f < nF; f++) for (let s = 0; s < nS; s++) { const v = CNT[f * nS + s]; ST[s] += v; FT[f] += v; }
    const host = $("#sample-bars");
    S.forEach((s, si) => {
      const lab = el("div", "sl"); lab.append(document.createTextNode(s.name)); lab.append(el("span", null, s.tissue));
      const segs = aggregate(fi => CNT[fi * nS + si]); const w = el("div"); stackedBar(w, segs, "of " + s.name);
      host.append(lab, w);
    });
    compLegend($("#sample-legend"), COMP_PRESENT.map(k => ({ key: k, value: null })));
    const c = { q: $("#sq"), comp: $("#scomp"), order: $("#sorder"), size: $("#ssize"), scale: $("#sscale"), clamp: $("#sclamp") };
    fillCompSelect(c.comp);
    const st = { rows: [] };
    const cnt = (r, s) => CNT[st.rows[r] * nS + s], lfc = (r, s) => LFC[st.rows[r] * nS + s] / 100;
    sdot = Dotplot($("#shost"), {
      ariaLabel: "Dotplot of factors by samples", cornerCols: "samples", cell: 26,
      rowLabel: r => { const f = F[st.rows[r]]; return [f.id, shortOf(f) || "", rowSwatch(f), fcolor(f)]; },
      colLabel: s => S[s].name, count: cnt, lfc: lfc, sizeMode: () => c.size.value, scale: () => c.scale.value, clamp: () => +c.clamp.value,
      rowTotal: r => FT[st.rows[r]], colTotal: s => ST[s],
      cellTip: (t, r, s) => {
        const f = F[st.rows[r]], v = cnt(r, s); tipTitle(t, "Factor " + f.id + " in " + S[s].name); if (f.label) t.append(el("div", "note", f.label));
        tipRow(t, fmtPct(ST[s] > 0 ? v / ST[s] : 0), "of the sample's counts"); tipRow(t, fmtPct(FT[st.rows[r]] > 0 ? v / FT[st.rows[r]] : 0), "of the factor's counts");
        tipRow(t, fmtLfc(lfc(r, s)), "log2 fold change vs other samples"); tipRow(t, fmtNum(v), "count");
      },
      rowTip: (t, r) => { const f = F[st.rows[r]]; tipTitle(t, "Factor " + f.id + " · " + nameOf(f)); tipRow(t, fmtNum(FT[st.rows[r]]), "decoded counts, all samples"); if (f.distribution) tipRow(t, f.distribution, "distribution"); t.append(el("div", "note", "click to open in the factor table")); },
      colTip: (t, s) => { tipTitle(t, S[s].name); tipRow(t, S[s].tissue, "tissue"); tipRow(t, fmtNum(ST[s]), "decoded counts"); },
      onRow: r => goFactor(F[st.rows[r]].id),
    });
    function update() {
      const q = c.q.value.trim().toLowerCase(), cp = c.comp.value;
      const rows = F.map((f, i) => i).filter(i => (!cp || compOf(F[i]) === cp) && (!q || searchText(F[i]).indexOf(q) >= 0));
      const o = c.order.value;
      if (o === "compartment") rows.sort((a, b) => (compRank(compOf(F[a])) - compRank(compOf(F[b]))) || (a - b));
      else if (o === "share") rows.sort((a, b) => (F[b].share || 0) - (F[a].share || 0));
      else if (o === "peak") { const pk = fi => { let b = 0, bs = 0; for (let s = 0; s < nS; s++) { const v = ST[s] > 0 ? CNT[fi * nS + s] / ST[s] : 0; if (v > b) { b = v; bs = s; } } return bs; }; const p = new Map(rows.map(fi => [fi, pk(fi)])); rows.sort((a, b) => (p.get(a) - p.get(b)) || (a - b)); }
      st.rows = rows; sdot.set(rows.length, nS);
      const hostL = $("#slegend"); hostL.replaceChildren(); const m = c.size.value;
      sizeKey(hostL, c.scale.value, m === "row" ? "the factor's max across samples" : m === "col" ? "the sample's max" : null, sdot.max());
      colorLegend(hostL, c.clamp.value, "Colour: log2 fold change vs other samples");
    }
    Object.entries(c).forEach(([k, e]) => e.addEventListener(k === "q" ? "input" : "change", update));
    update();
  }

  // ---------- run details ----------
  function renderRun() {
    const kv = $("#runkv");
    const add = (k, v) => { if (v == null || v === "") return; kv.append(el("dt", null, k), el("dd", null, String(v))); };
    add("Prefix", M.prefix); add("Factor model", M.model_file); add("DE table", M.de_file); add("Factor colours", M.rgb_file || "none (no rgb.tsv)"); add("Mode", M.mode);
    add("Organism", M.organism); add("Tissue context", M.tissues.join("; "));
    add("LLM", M.dry_run ? "none (dry run)" : M.api_type + " · " + M.model + " · effort " + M.effort); add("Prompt version", M.prompt_version);
    add("Annotation mode", M.annotation_mode); add("Candidate labels", M.n_candidates ? M.n_candidates + (M.closed ? " (closed list)" : " (suggestions)") : "none");
    add("Marker genes per factor", "top " + M.top_n + " by " + M.primary_rank + " and by " + M.secondary_rank + ", interleaved");
    add("LLM request/response files", M.llm_dir); add("Generated", M.date);
    const tb = $("#calls tbody");
    if (!D.calls.length) { const tr = el("tr"), td = el("td", "muted", "No LLM calls (dry run). Prompts were written to " + M.llm_dir + "."); td.colSpan = 6; tr.append(td); tb.append(tr); }
    D.calls.forEach(cl => { const u = cl.usage || {}; const tr = el("tr");
      tr.append(el("td", null, cl.name + (cl.reused ? " (reused)" : "")), el("td", "num", String(cl.n_factors)), el("td", null, cl.model || ""),
        el("td", "num", u.input_tokens != null ? u.input_tokens.toLocaleString() : ""), el("td", "num", u.output_tokens != null ? u.output_tokens.toLocaleString() : ""), el("td", null, cl.file)); tb.append(tr); });
    if (M.command) $("#cmd").textContent = M.command; else $("#cmdcard").remove();
  }

  // ---------- theme ----------
  const themeBtn = $("#theme"), MODES = ["auto", "light", "dark"];
  let theme = store.get("factor-report-theme") || document.documentElement.dataset.theme; if (MODES.indexOf(theme) < 0) theme = "auto";
  function redrawAll() { if (gdot) gdot.redraw(); if (sdot) sdot.redraw(); }
  function applyTheme() { if (theme === "auto") delete document.documentElement.dataset.theme; else document.documentElement.dataset.theme = theme; themeBtn.textContent = "Theme: " + theme; redrawAll(); }
  themeBtn.addEventListener("click", () => { theme = MODES[(MODES.indexOf(theme) + 1) % 3]; store.set("factor-report-theme", theme); applyTheme(); });
  try { matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redrawAll); } catch (e) { /* old browsers */ }

  renderOverview(); renderFactors(); renderRun();
  let samplesReady = false;
  onShow.dotplot = ensureGeneDot;
  onShow.samples = () => { if (!samplesReady) { samplesReady = true; initSamples(); } };
  applyTheme();
  showTab((location.hash || "").slice(1) || "overview");
})();
</script>
</body>
</html>
"""


# -----------------------------
# Main
# -----------------------------
def default_title(prefix, out):
    if not prefix:
        return os.path.basename(out).removesuffix(".tsv")
    parent = os.path.basename(os.path.dirname(os.path.abspath(prefix)))
    if parent == "cartl":
        parent = os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(prefix))))
    return f"{parent} · {os.path.basename(prefix)}"


def annotate_factors_with_llm(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Annotate the latent factors of a cartloader run with a generative LLM that labels all factors "
                    "together, and write an interactive HTML report.")
    io = parser.add_argument_group("Input/Output Parameters")
    io.add_argument("--prefix", type=str,
                    help="Factor prefix, e.g. cartl/<sample>/t18-f192 (reads <prefix>-bulk-de.tsv or -de.tsv, and "
                         "-model.tsv) or a multi-sample root cartl/t12-f192 (also reads each sample's pixel pseudobulk)")
    io.add_argument("--model", type=str, help="Factor model matrix TSV (gene x factor); default: <prefix>-model.tsv")
    io.add_argument("--de", type=str, help="Bulk DE TSV (gene, factor, Chi2, FoldChange, ...); default: "
                                           "<prefix>-bulk-de.tsv or <prefix>-de.tsv; computed from --model if absent")
    io.add_argument("--rgb", type=str, help=f"Factor colours (cartloader rgb.tsv, as shown in CartoScope); default: "
                                            f"<prefix>{RGB_SUFFIX} if present")
    io.add_argument("--out", type=str, help=f"Output alias TSV (index, alias); default: <prefix>{ALIAS_SUFFIX}. The "
                                            "companion files share its stem (.annotations.tsv, .html, .llm/)")
    io.add_argument("--organism", required=True, type=str, help="e.g. human, mouse")
    io.add_argument("--tissue", type=str, help="Tissue context of the dataset (all samples unless --sample-sheet gives "
                                               "a `tissue` per sample)")
    io.add_argument("--sample-sheet", "--samples", dest="sample_sheet", type=str,
                    help="Multi-sample: TSV with `id` and optional `tissue` (and `pseudobulk`) columns, e.g. the "
                         "run_together sample sheet. Samples are auto-discovered under --prefix without it")
    io.add_argument("--candidates", type=str, nargs="+",
                    help="Candidate label JSON file(s) from generate_candidate_labels_with_llm.py, or 'auto' to "
                         "generate (and cache) one list per tissue (optional)")
    io.add_argument("--api-type", required=True, choices=list(API_KEY_ENV), help="LLM provider")

    aux = parser.add_argument_group("Auxiliary Parameters")
    aux.add_argument("--single-sample", action="store_true", help="Treat --prefix as single-sample even if sample "
                                                                   "directories are found next to it")
    aux.add_argument("--multi-catalog", default="multi-catalog.yaml", help="Multi-sample catalog next to the prefix, "
                                                                            "used to find samples (default: multi-catalog.yaml)")
    aux.add_argument("--mode", default="either", choices=list(MODE_DESCRIPTIONS), help="What kind of label to assign")
    aux.add_argument("--closed", action="store_true", help="Restrict labels to the candidate labels")
    aux.add_argument("--candidate-cache-dir", default="~/.cache/cartloader/llm_candidates",
                     help="Cache of candidate lists for --candidates auto (default: ~/.cache/cartloader/llm_candidates)")
    aux.add_argument("--candidate-effort", default="low", help="Reasoning effort for --candidates auto (default: low)")
    aux.add_argument("--primary-rank", default="Chi2", help="Primary ranking column (default: Chi2)")
    aux.add_argument("--secondary-rank", default="FoldChange", help="Secondary ranking column (default: FoldChange)")
    aux.add_argument("--top-n", type=int, default=10, help="Top genes per ranking column (default: 10)")
    aux.add_argument("--min-marker-gene-count", type=int, default=5,
                     help="Factors with fewer marker genes are left Unresolved without an API call (default: 5)")
    aux.add_argument("--no-specificity", action="store_true", help="Omit per-gene specificity (needs the model matrix)")
    aux.add_argument("--no-sample-markers", action="store_true", help="Multi-sample: omit per-sample marker genes")
    aux.add_argument("--max-factors-per-call", type=int, default=0,
                     help="Factors per LLM call; 0 = all factors in one call (default). A call whose answer overflows "
                          "the output budget is split in halves automatically")
    aux.add_argument("--model-name", type=str, help=f"LLM model (defaults: {DEFAULT_MODEL})")
    aux.add_argument("--effort", default="high", help="Reasoning effort: low, medium, high, xhigh, max (default: high)")
    aux.add_argument("--max-output-tokens", type=int,
                     help=f"Output token budget per call (defaults: {DEFAULT_MAX_OUTPUT_TOKENS})")
    aux.add_argument("--no-fallbacks", action="store_true", help="claude: disable server-side refusal fallbacks")
    aux.add_argument("--api-base-url", type=str, help="Base URL for openai/umgpt (e.g. a proxy)")
    aux.add_argument("--request-timeout", type=int, default=1800, help="Request timeout in seconds (default: 1800)")
    aux.add_argument("--max-retries", type=int, default=3, help="Retries for failed requests (default: 3)")
    aux.add_argument("--threads", type=int, default=4, help="Parallel LLM calls when there are several (default: 4)")
    aux.add_argument("--report-only", action="store_true",
                     help="Rebuild the report and TSVs from the responses saved in <stem>.llm/ without calling the LLM "
                          "(each factor takes its newest saved answer, even if the prompt has changed since)")
    aux.add_argument("--dry-run", action="store_true",
                     help="Write the prompts to <stem>.llm/ and an evidence-only report, without calling the LLM")
    aux.add_argument("--html", type=str, help="Report path (default: <stem>.html)")
    aux.add_argument("--no-html", action="store_true", help="Do not write the report")
    aux.add_argument("--dotplot-genes", default="markers", choices=["markers", "all"],
                     help="Genes embedded in the report's dotplot: marker/key genes of all factors (default) or all "
                          "genes of the model (larger file)")
    aux.add_argument("--title", type=str, help="Report title (default: dataset and prefix name)")
    aux.add_argument("--no-reuse", action="store_true", help="Do not reuse responses saved in <stem>.llm/")
    aux.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args(_args)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="[%(asctime)s - %(levelname)s - %(message)s]", datefmt="%Y-%m-%d %H:%M:%S")
    if args.prefix:
        args.model = args.model or first_existing(args.prefix, MODEL_SUFFIXES)
        args.de = args.de or first_existing(args.prefix, DE_SUFFIXES)
        args.rgb = args.rgb or first_existing(args.prefix, (RGB_SUFFIX,))
        args.out = args.out or args.prefix + ALIAS_SUFFIX
    if not (args.model or args.de):
        parser.error(f"no factor model or DE table: pass --prefix (looked for {MODEL_SUFFIXES} and {DE_SUFFIXES}), "
                     "--model or --de")
    if not args.out:
        parser.error("--out is required without --prefix")
    args.model_name = args.model_name or DEFAULT_MODEL[args.api_type]
    args.max_output_tokens = args.max_output_tokens or DEFAULT_MAX_OUTPUT_TOKENS[args.api_type]
    if args.dry_run and args.report_only:
        parser.error("--dry-run and --report-only are mutually exclusive")
    if not (args.dry_run or args.report_only):
        env = API_KEY_ENV[args.api_type]
        if not os.environ.get(env) and not (args.api_type == "claude" and os.environ.get("ANTHROPIC_AUTH_TOKEN")):
            log.error("%s is not set (needed for --api-type %s)", env, args.api_type)
            sys.exit(1)

    samples = resolve_samples(args)
    if len(samples) == 1:
        log.info("Only one sample (%s); annotating as a single-sample dataset", samples[0]["name"])
        args.tissue = args.tissue or samples[0]["tissue"]
        samples = []
    if not samples and not args.tissue:
        parser.error("--tissue is required for a single-sample dataset")
    stem = args.out.removesuffix(".tsv")
    os.makedirs(os.path.dirname(os.path.abspath(stem)), exist_ok=True)

    log.info("Reading %s%s and computing marker genes and evidence", args.model or args.de,
             f" and {len(samples)} sample pseudobulks" if samples else "")
    bundle = build_profiles(args, samples)
    if args.candidates == ["auto"]:
        if args.report_only:  # no LLM calls: use the list saved by the earlier run, if any
            saved = stem + ".candidates.json"
            args.candidates = [saved] if os.path.exists(saved) else None
        else:
            args.candidates = [generate_candidates(args, bundle.tissues, stem + ".candidates.json")]
    candidates = load_candidates(args.candidates, bundle.genes) if args.candidates else []
    usable = [p for p in bundle.profiles if len(p["marker_genes"]) >= args.min_marker_gene_count]
    ctx = prompt_context(usable, samples, bundle.tissues, candidates, args.organism, args.mode, args.closed)
    reqs = [build_request(ctx, b, f"batch{i:02d}")
            for i, b in enumerate(make_batches(usable, args.max_factors_per_call), 1)]
    log.info("%d factors (%d with enough marker genes), %s, %d candidate labels, %s",
             len(bundle.profiles), len(usable),
             f"{len(samples)} samples from {len(bundle.tissues)} tissue context(s)" if samples else "single sample",
             len(candidates), "no LLM calls (report only)" if args.report_only else f"{len(reqs)} LLM call(s)")

    llm_dir = stem + ".llm"
    if not args.report_only:
        os.makedirs(llm_dir, exist_ok=True)
    raw, summaries, calls = {}, [], []
    if args.report_only:
        try:
            raw, summaries, calls, models, efforts = load_saved_responses(llm_dir, {p["factor"] for p in usable})
        except FileNotFoundError as e:
            log.error("%s", e)
            sys.exit(1)
        args.model_name = ", ".join(models) or args.model_name
        args.effort = ", ".join(efforts) or args.effort
        log.info("Report only: %d of %d annotated factors have a saved answer", len(raw), len(usable))
    elif args.dry_run:
        for name, system, user, _, _ in reqs:
            with open(os.path.join(llm_dir, f"{name}.prompt.txt"), "w", encoding="utf-8") as fh:
                fh.write(f"# SYSTEM\n{system}\n\n# USER\n{user}\n")
        log.info("Dry run: wrote %d prompt(s) to %s", len(reqs), llm_dir)
    else:
        with ThreadPoolExecutor(max_workers=max(1, min(args.threads, len(reqs) or 1))) as ex:
            for got in ex.map(lambda r: run_request(args, llm_dir, ctx, r), reqs):
                raw.update({str(a["factor"]): a for a in got["annotations"]})
                summaries += got["summaries"]
                calls += got["calls"]

    table = to_annotations(bundle.profiles, usable, raw, {c["label"] for c in candidates},
                           pending="dry run: no LLM call" if args.dry_run else None)
    if not args.dry_run:
        table.to_csv(stem + ".annotations.tsv", sep="\t", index=False)
        aliases = suffix_duplicates(list(zip(table["factor"], table["alias"])))
        pd.DataFrame(aliases, columns=["index", "alias"]).to_csv(args.out, sep="\t", index=False)
        log.info("Wrote %s and %s (%s)", args.out, stem + ".annotations.tsv", table["status"].value_counts().to_dict())
    if not args.no_html:
        html_path = args.html or stem + ".html"
        meta = {"title": args.title or default_title(args.prefix, args.out), "date": dt.date.today().isoformat(),  # noqa: DTZ011
                "organism": args.organism, "tissues": bundle.tissues,
                "mode": "multi-sample" if samples else "single-sample", "api_type": args.api_type,
                "model": args.model_name, "effort": args.effort, "prompt_version": PROMPT_VERSION,
                "annotation_mode": args.mode, "n_candidates": len(candidates), "closed": args.closed,
                "dry_run": args.dry_run, "prefix": args.prefix, "stem": os.path.basename(stem),
                "model_file": args.model, "de_file": args.de or "computed from the factor model", "rgb_file": args.rgb,
                "llm_dir": llm_dir, "top_n": args.top_n, "primary_rank": args.primary_rank,
                "secondary_rank": args.secondary_rank,
                "command": "cartloader annotate_factors_with_llm " + " ".join(shlex.quote(a) for a in _args)}
        colors = read_rgb(args.rgb) if args.rgb else {}
        missing = [p["factor"] for p in bundle.profiles if p["factor"] not in colors]
        if colors and missing:
            log.warning("%d factor(s) have no colour in %s", len(missing), args.rgb)
        payload = report_payload(bundle, table, raw, samples, summaries, calls, candidates, meta, args.dotplot_genes, colors)
        with open(html_path, "w", encoding="utf-8") as fh:
            fh.write(render_html(payload))
        log.info("Wrote %s", html_path)


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)
    # Call the function with command line arguments
    func(sys.argv[1:])
