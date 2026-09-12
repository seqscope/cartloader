import sys, os, argparse, inspect, json, copy, csv, datetime, glob

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import execute_makefile

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROFILE_DIR = os.path.join(repo_dir, "assets", "run_together_profiles")
SPATULA_BIN = os.path.join(repo_dir, "submodules", "spatula", "bin", "spatula")
IMAGE_TYPES_FILE = os.path.join(repo_dir, "assets", "run_together_image_types.json")

# Global fallbacks (a profile or --config may override; a CLI flag wins over both).
# The two feature-exclusion defaults split the work by what the feature is:
#  * technical artifacts (negative controls, blanks, unassigned/deprecated codewords) have
#    no biological meaning at all, so they are dropped at ingest and never enter the data;
#  * genes that are real but dominate or distort a factorization (predicted Gm* models,
#    mitochondrial, ribosomal) stay in the data — visible in the browser, packaged in the
#    tiles — and are only kept out of the model fitting.
# Both apply to every platform; a profile or --config may override either per run.
DEFAULT_INGEST_EXCLUDE_REGEX = "^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated|System|NCS-|NCP-)"
DEFAULT_EXCLUDE_REGEX = "^(Gm[0-9]|MT-|mt-|Rps|Rpl)"
DEFAULT_MIN_CT_PER_UNIT_HEXAGON = 50
# Sample id for a single-sample --out-dir run that names no id. Deliberately generic:
# the descriptive name lives in --out-dir, and the packaged directory / catalog id
# compose as <out_dir_basename>-<sample_id>, matching joint runs (see resolve_sample).
DEFAULT_SAMPLE_ID = "rep1"

# Publish (S3 upload) defaults for the CartoStore project.
DEFAULT_S3_PREFIX = "s3://cartostore/data"
DEFAULT_S3_PROFILE = "default"

# Sample-level input roles. A sample provides these either explicitly (sample
# sheet columns / JSON), by auto-detection inside its `in_dir` (a profile role's
# `file`), or by suffixing its `in_prefix` (a profile role's `suffix`).
ROLE_KEYS = ["transcript", "xy", "boundaries", "clusters", "mex", "cellxgene",
             "cell_tsv", "gef", "cellbin_gef"]
# Friendly sample-sheet column aliases (the generic platform advertises these
# names; each maps onto a canonical role key above).
ROLE_ALIASES = {
    "tsv": "transcript",
    "cell_xy": "xy",
    "cell_boundary": "boundaries",
    "cell_boundaries": "boundaries",
    "mex_dir": "mex",
    "cell_by_gene": "cellxgene",
}
# role -> the run_ficture2_multi_cells flag that consumes it
ROLE_LIST_FLAG = {
    "boundaries": "--list-boundaries",
    "xy": "--list-xy",
    "clusters": "--list-cluster",
    "mex": "--mex-list",
    "cell_tsv": "--tsv-list",
}


def explicit_lists(ca):
    """Ready-made --list-* files named for a cell analysis, as role -> path. run_together
    normally derives each list from the samples' resolved roles; naming one supplies it
    verbatim instead, either per analysis (`lists` in a cell_analyses entry) or run-wide
    (top-level `cell_lists` / the --list-* CLI flags, merged in by build_config).

    An explicitly named list also satisfies that role's requirement in the analysis's
    uses/any_uses gating, and is passed even for a role the analysis does not declare at
    all: the list was supplied directly, so no sample has to carry the role on disk. That
    is what makes externally assigned cluster labels reachable on a multi-sample run —
    without --list-cluster, run_ficture2_multi_cells computes Leiden clusters on demand,
    which is not always what the run wants. Same escape hatch for a curated xy/boundaries
    list that should replace the auto-detected per-sample files.
    """
    return ca.get("lists") or {}
# Recognized top-level keys in a --config JSON. Anything else aborts build_config
# (unless --allow-unknown-config-keys), since unknown keys are otherwise silently
# dropped. Nested/per-file overrides (e.g. csv_colnames) live under these, not here.
KNOWN_CONFIG_KEYS = frozenset({
    # sample selection / output
    "platform", "samples", "out_dir", "out_root", "saw",
    # scalar decode override consumed directly from cfg
    "colname_cell", "exclude_feature_regex", "include_feature_list", "exclude_feature_list",
    "ingest_include_feature_regex", "ingest_exclude_feature_regex",
    "ingest_include_feature_list", "ingest_exclude_feature_list",
    # ready-made --list-* files applied to every cell analysis (role -> path)
    "cell_lists",
    # dict blocks deep-merged into the profile (Layer 2)
    "ingest", "roles", "cartload", "ficture_defaults", "cell_defaults", "squares", "cell_import",
    "hne", "image_transform", "image_defaults", "publish", "resources",
    # list blocks merged by id
    "ficture", "cell_analyses", "images",
    # 10x MEX export of the hexagon files: true, or { "widths": "12,24" }
    "segment_10x",
})
# Recognized keys inside a `samples[]` entry. Anything else is a typo that would be
# silently dropped (a misspelled role looks provided but is never read), so build_config
# rejects it unless --allow-unknown-config-keys is set. Covers the canonical roles, their
# friendly aliases, the explicit mex triple, per-sample images, and the addressing keys.
KNOWN_SAMPLE_KEYS = frozenset(
    {"id", "in_dir", "in_prefix", "raw_transcript", "images", "dapi", "hne",
     "mex_bcd", "mex_ftr", "mex_mtx"}
    | set(ROLE_KEYS) | set(ROLE_ALIASES)
)
# Recognized sub-keys for the closed-schema nested config blocks. An unknown key here is
# deep_merged into the profile and then never read, so a typo (e.g. "bin_counts") looks
# applied but has no effect; build_config rejects it (same --allow-unknown-config-keys
# escape hatch). Only blocks with a fixed key schema are listed: open maps whose sub-keys
# are data, not schema (roles, csv_colnames, ficture_defaults per-entry fields, publish,
# hne, image_transform, cell_import, squares), are deliberately omitted.
KNOWN_SUBKEYS = {
    "ingest": frozenset({
        "method", "sge_platform", "units_per_um", "gef_suffix", "feature_file",
        "in_dir_flag", "raw_input_flag", "input_roles", "inputs", "produces",
        "autodetect", "assign_cell_id", "csv_colname_cell", "csv_colnames",
        "csv_colnames_others", "cellxgene_blank_prefix", "jitter_xy", "extra_flags",
    }),
    "cartload": frozenset({"use_pmpoint", "sge_scale", "bin_count"}),
    "resources": frozenset({"threads", "n_jobs"}),
    "image_defaults": frozenset({
        "um_per_pixel", "um_per_pixel_json", "um_per_pixel_key", "georeferenced",
        "georef_detect", "shrink_factor", "high_memory", "rescale", "rescale_range",
        "rescale_min", "rescale_max",
        "swapxy", "rotate", "flip_vertical", "flip_horizontal",
    }),
}
# transcript column-override key -> the sge_convert flag that names that input column
CSV_COLNAME_FLAGS = {
    "x": "--csv-colname-x",
    "y": "--csv-colname-y",
    "feature": "--csv-colname-feature-name",
    "count": "--csv-colname-count",
}
# Image modality registry: `type` -> its default colorize hex and kind. A single
# (colorized) type's `color` is overridable per image; an `rgb` type (e.g. hne) is a
# multi-channel passthrough (no colorize). Loaded from assets/ so it can be edited
# without touching code; keys beginning with '_' (e.g. "_comment") are ignored.
def load_image_types():
    if not os.path.exists(IMAGE_TYPES_FILE):
        sys.exit(f"ERROR: image type registry not found: {IMAGE_TYPES_FILE}")
    with open(IMAGE_TYPES_FILE) as f:
        return {k: v for k, v in json.load(f).items() if not k.startswith("_")}

IMAGE_TYPES = load_image_types()

# ---------------------------------------------------------------------------
# Merge helpers (the three-layer assembly: profile -> CLI -> JSON)
# ---------------------------------------------------------------------------

def deep_merge(base, override):
    """Recursively merge dicts; non-dict values (incl. lists) are replaced."""
    if not isinstance(base, dict) or not isinstance(override, dict):
        return copy.deepcopy(override)
    out = copy.deepcopy(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def merge_id_list(base, override):
    """Merge two id-keyed lists: same id -> deep-merge, new id -> append.

    ``override`` may instead be ``{"replace": [...]}`` to discard the base list.
    """
    if isinstance(override, dict):
        if "replace" in override:
            return copy.deepcopy(override["replace"])
        override = override.get("add", [])
    out = copy.deepcopy(base)
    index = {e["id"]: i for i, e in enumerate(out) if "id" in e}
    for e in override:
        eid = e.get("id")
        if eid is not None and eid in index:
            out[index[eid]] = deep_merge(out[index[eid]], e)
        else:
            out.append(copy.deepcopy(e))
            if eid is not None:
                index[eid] = len(out) - 1
    return out


def merge_images(base, override):
    """Images allow same-id fallback entries. An override id replaces every base
    entry sharing that id, then the overrides are appended."""
    if isinstance(override, dict):
        if "replace" in override:
            return copy.deepcopy(override["replace"])
        override = override.get("add", [])
    override_ids = {e["id"] for e in override if "id" in e}
    kept = [e for e in base if e.get("id") not in override_ids]
    return copy.deepcopy(kept) + copy.deepcopy(override)


def load_json(path):
    with open(path) as f:
        return json.load(f)


def load_builtin_profile(platform):
    path = os.path.join(PROFILE_DIR, f"{platform}.json")
    if not os.path.exists(path):
        avail = sorted(os.path.splitext(f)[0] for f in os.listdir(PROFILE_DIR) if f.endswith(".json"))
        sys.exit(f"ERROR: no built-in profile for platform '{platform}'. "
                 f"Available: {', '.join(avail)}. Provide one via --profile.")
    return load_json(path)


def first_existing(*paths):
    for p in paths:
        if p and os.path.exists(p):
            return p
    return None


SHEET_UNSET = {"", "-", ".", "NA"}


def read_sheet(path):
    """Read a wide sample sheet (TSV). Columns are input roles; the values
    '', '-', '.', 'NA' all mean unset."""
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            clean = {}
            for k, v in row.items():
                if k is None:
                    continue
                v = (v or "").strip()
                if v not in SHEET_UNSET:
                    clean[k.strip()] = v
            if clean:
                rows.append(clean)
    return rows


# ---------------------------------------------------------------------------
# FICTURE analysis resolution (the two tier-1 modes)
# ---------------------------------------------------------------------------

def build_projection_analyses(project_dirs, width):
    """Projection-only mode: read each existing FICTURE dir's ficture.params.json
    and reuse every trained model as a projection analysis."""
    analyses = []
    for d in [x for x in project_dirs.split(",") if x]:
        params_path = os.path.join(d, "ficture.params.json")
        if not os.path.exists(params_path):
            sys.exit(f"ERROR: --project-models expects a FICTURE directory containing "
                     f"ficture.params.json; not found: {params_path}")
        data = load_json(params_path)
        for tp in data.get("train_params", []):
            analyses.append({
                "id": tp["model_id"],
                "mode": "project",
                "model": tp["model_path"],
                "width": str(width) if width else str(tp.get("train_width")),
                "n_factor": tp.get("n_factor"),
            })
    if not analyses:
        sys.exit("ERROR: --project-models found no models in the referenced ficture.params.json file(s).")
    return analyses


def resolved_models(analyses):
    """Return (list of {model_id, n_factor}, default_model_id).

    A de-novo entry with widths W and factors F expands to model ids t{w}_f{f};
    a projection entry contributes its own id.
    """
    models = []
    for a in analyses:
        if a.get("mode") == "prepare":
            continue           # tiling only (--no-ficture): no models exist
        if a.get("mode") == "project":
            models.append({"model_id": a["id"], "n_factor": int(a.get("n_factor") or 0)})
        else:
            widths = str(a["width"]).split(",")
            factors = str(a["n_factor"]).split(",")
            for w in widths:
                for nf in factors:
                    models.append({"model_id": f"t{w}_f{nf}", "n_factor": int(nf)})
    default = max(models, key=lambda m: m["n_factor"])["model_id"] if models else None
    return models, default


# ---------------------------------------------------------------------------
# Configuration assembly
# ---------------------------------------------------------------------------

def _consumable_roles(ingest, cell_analyses):
    """Roles this platform can consume, given its ingest method and cell analyses. A
    sample role outside this set is read by no code path on this platform, so it is inert.

    Deliberately a question of platform *capability*, not of one run: `cell_analyses` here
    is the platform's (profile + config) declared set, NOT the run-mode-mutated one, so a
    --no-ficture / prepare-only run (which skips cells) does not make its cell roles look
    invalid. Verified against every role consumption site to be exact for a normal run
    (never under-reports), which is what lets the caller treat a miss as a hard error."""
    method = ingest.get("method")
    roles = {"transcript"}                                  # always the pixel transcript
    roles |= set(ingest.get("input_roles", {}).values())   # e.g. illumina mex -> --in-mex
    if method == "stereoseq":
        roles |= {"gef", "cellbin_gef", "cell_tsv"}         # SAW ingest + its cellbin cells
    if ingest.get("assign_cell_id"):
        roles.add("boundaries")                             # tsv-add-cell-id
    if method not in ("stereoseq", "reformat_cosmx"):
        roles.add("cellxgene")                              # cellxgene -> mex on the generic path
    for ca in cell_analyses:
        roles |= set(ca.get("uses", [])) | set(ca.get("any_uses", [])) | set(ca.get("optional_uses", []))
        if ca.get("generic_cell"):
            roles.add("mex")
    return roles & set(ROLE_KEYS)


def build_config(args):
    cfg = load_json(args.config) if args.config else {}

    # Reject unrecognized top-level config keys: they are otherwise silently dropped,
    # so a misplaced key (e.g. a top-level "csv_colnames" that belongs under "ingest")
    # looks applied but has no effect. --allow-unknown-config-keys downgrades this to
    # the old ignore-and-continue behavior.
    unknown = [k for k in cfg if k not in KNOWN_CONFIG_KEYS]
    if unknown and not args.allow_unknown_config_keys:
        sys.exit(f"ERROR: unrecognized top-level key(s) in --config: {', '.join(sorted(unknown))}. "
                 f"Recognized keys: {', '.join(sorted(KNOWN_CONFIG_KEYS))}. "
                 f"(Per-file column overrides like 'csv_colnames' go under 'ingest'.) "
                 f"Pass --allow-unknown-config-keys to ignore unknown keys instead of failing.")

    # Reject unknown keys one level down, in the same spirit: a misspelled key inside a
    # sample entry or a closed-schema block is deep_merged in and then never read, so it
    # looks applied but silently does nothing.
    if not args.allow_unknown_config_keys:
        for block, allowed in KNOWN_SUBKEYS.items():
            bad = [k for k in (cfg.get(block) or {}) if k not in allowed]
            if bad:
                sys.exit(f"ERROR: unrecognized key(s) in --config '{block}': {', '.join(sorted(bad))}. "
                         f"Recognized keys: {', '.join(sorted(allowed))}. "
                         f"Pass --allow-unknown-config-keys to ignore unknown keys instead of failing.")
        for i, samp in enumerate(cfg.get("samples") or []):
            if not isinstance(samp, dict):
                sys.exit(f"ERROR: --config 'samples[{i}]' must be an object, got {type(samp).__name__}.")
            bad = [k for k in samp if k not in KNOWN_SAMPLE_KEYS]
            if bad:
                sid = samp.get("id", f"index {i}")
                sys.exit(f"ERROR: unrecognized key(s) in --config sample '{sid}': {', '.join(sorted(bad))}. "
                         f"Recognized keys: {', '.join(sorted(KNOWN_SAMPLE_KEYS))}. "
                         f"Pass --allow-unknown-config-keys to ignore unknown keys instead of failing.")

    platform = args.platform or cfg.get("platform")
    if not platform:
        sys.exit("ERROR: --platform (or a 'platform' field in --config) is required.")

    # Layer 0: profile
    prof = load_builtin_profile(platform)
    if args.platform_json:
        prof = deep_merge(prof, load_json(args.platform_json))
    prof["platform"] = platform
    prof.setdefault("ficture_defaults", {})
    prof.setdefault("cell_defaults", {})
    prof.setdefault("ficture", [])
    prof.setdefault("cell_analyses", [])
    prof.setdefault("images", [])
    prof.setdefault("roles", {})
    # The platform's structural cell analyses, captured before Layer 1 can wipe them
    # (--no-ficture). Used only for the inert-role capability check, so a prepare-only run
    # does not make otherwise-valid cell roles look invalid. The actual run still uses the
    # (possibly cleared) prof["cell_analyses"].
    _profile_cell_analyses = list(prof["cell_analyses"])

    # Layer 1: tier-1 CLI selects the base FICTURE mode
    if args.no_ficture:
        if args.project_models or args.n_factor:
            sys.exit("ERROR: --no-ficture runs no factor analysis, so it cannot be combined "
                     "with --project-models or --n-factor.")
        # Only the tiling step runs. A width is still needed: prepare builds the hexagon
        # files alongside the tiles, so a later full FICTURE run in the same directory
        # resumes from them instead of re-tiling.
        width = args.width or next((str(a["width"]) for a in prof["ficture"] if a.get("width")), "12")
        prof["ficture"] = [{"id": "prepare", "mode": "prepare", "width": width}]
        prof["cell_analyses"] = []
    elif args.project_models:
        prof["ficture"] = build_projection_analyses(args.project_models, args.width)
    elif args.n_factor or args.width:
        # override the (single) de-novo entry, or create one
        train = next((a for a in prof["ficture"] if a.get("mode", "train") == "train"), None)
        if train is None:
            train = {"id": "denovo", "mode": "train"}
            prof["ficture"].append(train)
        if args.width:
            train["width"] = args.width
        if args.n_factor:
            train["n_factor"] = args.n_factor

    # Layer 2: tier-2 JSON augments (lists append/override-by-id by default)
    for k in ("exclude_feature_regex", "include_feature_list", "exclude_feature_list",
              "ingest_include_feature_regex", "ingest_exclude_feature_regex",
              "ingest_include_feature_list", "ingest_exclude_feature_list",
              "ingest", "roles", "cartload", "ficture_defaults", "cell_defaults",
              "squares", "cell_import", "hne", "image_transform", "image_defaults",
              "publish", "resources", "segment_10x"):
        if k in cfg:
            prof[k] = deep_merge(prof.get(k), cfg[k]) if isinstance(cfg[k], dict) else cfg[k]
    if "ficture" in cfg:
        prof["ficture"] = merge_id_list(prof["ficture"], cfg["ficture"])
    if "cell_analyses" in cfg:
        prof["cell_analyses"] = merge_id_list(prof["cell_analyses"], cfg["cell_analyses"])
    if "images" in cfg:
        prof["images"] = merge_images(prof["images"], cfg["images"])

    # Reject explicit sample roles the platform can never consume (e.g. cell_tsv on
    # 10x_xenium): the value resolves but no ingest step or cell analysis on this platform
    # reads it, so it silently has no effect — almost always a wrong --platform, a copy-paste
    # leftover, or a typo. Judged by platform *capability* (profile + config cell analyses,
    # not the run-mode-mutated set) so a --no-ficture / prepare-only run does not fail on
    # cell roles it merely skips. Checks explicit JSON `samples` keys only (auto-detected
    # roles are not second-guessed). Same --allow-unknown-config-keys escape hatch.
    if not args.allow_unknown_config_keys:
        cap_cas = (merge_id_list(_profile_cell_analyses, cfg["cell_analyses"])
                   if "cell_analyses" in cfg else _profile_cell_analyses)
        consumable = _consumable_roles(prof.get("ingest", {}), cap_cas)
        role_of = {**{r: r for r in ROLE_KEYS}, **ROLE_ALIASES}
        for samp in (cfg.get("samples") or []):
            inert = sorted({role_of[k] for k in samp
                            if k in role_of and role_of[k] not in consumable})
            if inert:
                sys.exit(f"ERROR: sample '{samp.get('id', '?')}' provides role(s) "
                         f"{', '.join(inert)} that platform '{platform}' cannot consume; they "
                         f"would have no effect (wrong --platform, or a leftover/misplaced role?). "
                         f"Recognized roles for this platform: {', '.join(sorted(consumable))}. "
                         f"Pass --allow-unknown-config-keys to ignore instead of failing.")

    # --- resolve common decode defaults (CLI > config/profile > hardcoded) ---
    # Feature filters. These apply to the factor analyses only (pixel FICTURE and the
    # cell-based analysis): ingest keeps every gene, so the transcript TSV, the feature
    # list and the packaged tiles stay complete regardless of what is filtered here.
    # Passing an empty string turns the default off without substituting anything.
    if args.exclude_feature_regex is not None:
        prof["exclude_feature_regex"] = args.exclude_feature_regex
    elif "exclude_feature_regex" not in prof:
        prof["exclude_feature_regex"] = DEFAULT_EXCLUDE_REGEX
    # Ingest-level filters are a separate, independent knob: they drop features from the
    # transcript TSV itself (the "very first process"), so a feature removed here is gone
    # from every later stage, packaging included.
    # Which of them this run actually asked for (rather than inherited from the default),
    # so that "these had no effect" is reported only for a filter someone chose. Recorded
    # here because `cfg` is merged into `prof` and the two become indistinguishable after.
    prof["_ingest_filters_explicit"] = [
        f"--{k.replace('_', '-')}"
        for k in ("ingest_include_feature_list", "ingest_exclude_feature_list",
                  "ingest_include_feature_regex", "ingest_exclude_feature_regex")
        if getattr(args, k) is not None or k in cfg]
    if args.ingest_exclude_feature_regex is not None:
        prof["ingest_exclude_feature_regex"] = args.ingest_exclude_feature_regex
    elif "ingest_exclude_feature_regex" not in prof:
        prof["ingest_exclude_feature_regex"] = DEFAULT_INGEST_EXCLUDE_REGEX
    if args.ingest_include_feature_regex is not None:
        prof["ingest_include_feature_regex"] = args.ingest_include_feature_regex
    for key, val in (("include_feature_list", args.include_feature_list),
                     ("exclude_feature_list", args.exclude_feature_list),
                     ("ingest_include_feature_list", args.ingest_include_feature_list),
                     ("ingest_exclude_feature_list", args.ingest_exclude_feature_list)):
        if val is not None:
            prof[key] = val
        if prof.get(key):
            path = os.path.abspath(os.path.expanduser(prof[key]))
            if not os.path.exists(path):
                sys.exit(f"ERROR: file not found: {prof[key]} (--{key.replace('_', '-')})")
            prof[key] = path
    # Ingest-time coordinate jitter, in um: each transcript's X and Y independently get a
    # uniform random offset in [-jitter, +jitter] as the transcript TSV is written. It is for
    # coarse-resolution platforms whose molecules all land on the same lattice points (e.g.
    # Visium HD's 2um bins): without it every transcript in a bin shares one coordinate, and
    # hexagon/pixel-level analysis and rendering see a grid rather than a distribution. Like
    # the ingest filters this changes the data itself, so every later stage sees the jittered
    # coordinates. A CLI flag wins over `ingest.jitter_xy` from the profile/config.
    if args.jitter_xy is not None:
        prof.setdefault("ingest", {})["jitter_xy"] = args.jitter_xy
    jitter = prof.get("ingest", {}).get("jitter_xy")
    if jitter is not None:
        try:
            jitter = float(jitter)
        except (TypeError, ValueError):
            sys.exit(f"ERROR: 'ingest.jitter_xy' (--jitter-xy) must be a number, got {jitter!r}.")
        if jitter < 0:
            sys.exit(f"ERROR: --jitter-xy is a half-width in um and must be >= 0 (got {jitter}).")
        prof["ingest"]["jitter_xy"] = jitter
    # Only the sge_convert ingest path can apply jitter. The CosMx reformat writes its
    # transcript TSV itself and has no jitter option, so asking for it there is an error
    # rather than a flag that is silently dropped.
    ingest_method = prof.get("ingest", {}).get("method", "sge_convert")
    if jitter and ingest_method == "reformat_cosmx":
        sys.exit(f"ERROR: --jitter-xy is not supported on platform '{platform}': its ingest "
                 f"('{ingest_method}') writes the transcript TSV without going through "
                 f"sge_convert. Drop the flag, or jitter the transcript yourself and supply it "
                 f"with --in-transcript.")
    # Stereo-seq ingests bin1 through sge_convert (so jitter applies) but writes its cell-bin
    # TSV with convert_stereoseq_cellbin, which has no jitter option.
    if jitter and ingest_method == "stereoseq":
        print("WARNING: --jitter-xy applies to the Stereo-seq bin1 pixel transcript only. The "
              "cell-bin TSV (the cell_tsv role) is written by convert_stereoseq_cellbin, which "
              "has no jitter option, so its coordinates stay on the original grid.",
              file=sys.stderr)
    # An explicit jitter, like an explicit ingest filter, silently does nothing for a sample
    # that supplies an already-ingested transcript; report it there alongside the filters.
    if jitter:
        prof["_ingest_filters_explicit"].append("--jitter-xy")
    # Coordinate units per um in the raw input, i.e. how many input units make one micron
    # (1000 for nanometer coordinates, 1 for microns, 2 for Stereo-seq's 500nm bins). Left
    # unset, each platform's sge_convert preset decides -- including Illumina, which reads
    # its barcode file to tell nanometer barcodes ("SBC:433503:2393851") from micron ones
    # ("SBC:686.951:4668.15"). A CLI flag wins over `ingest.units_per_um`.
    if args.units_per_um is not None:
        prof.setdefault("ingest", {})["units_per_um"] = args.units_per_um
    units_per_um = prof.get("ingest", {}).get("units_per_um")
    if units_per_um is not None:
        try:
            units_per_um = float(units_per_um)
        except (TypeError, ValueError):
            sys.exit(f"ERROR: 'ingest.units_per_um' (--units-per-um) must be a number, got "
                     f"{units_per_um!r}.")
        if units_per_um <= 0:
            sys.exit(f"ERROR: --units-per-um is input units per micron and must be > 0 (got "
                     f"{units_per_um}).")
        prof["ingest"]["units_per_um"] = units_per_um
        # Same story as jitter: only the sge_convert ingest path takes a coordinate scale.
        if ingest_method == "reformat_cosmx":
            sys.exit(f"ERROR: --units-per-um is not supported on platform '{platform}': its "
                     f"ingest ('{ingest_method}') writes the transcript TSV without going "
                     f"through sge_convert. Drop the flag, or rescale the transcript yourself "
                     f"and supply it with --in-transcript.")
        # Only a value this run chose is worth reporting as "had no effect" on a
        # pre-ingested sample; a profile's own default (e.g. Stereo-seq's 2) is not.
        if args.units_per_um is not None or "units_per_um" in cfg.get("ingest", {}):
            prof["_ingest_filters_explicit"].append("--units-per-um")
    # min count per unit hexagon / per unit trained (both apply to the pixel FICTURE analyses)
    fd = prof.setdefault("ficture_defaults", {})
    if args.min_ct_per_unit_hexagon is not None:
        fd["min_ct_per_unit_hexagon"] = args.min_ct_per_unit_hexagon
    elif "min_ct_per_unit_hexagon" not in fd:
        fd["min_ct_per_unit_hexagon"] = DEFAULT_MIN_CT_PER_UNIT_HEXAGON
    if args.min_ct_per_unit_train is not None:
        fd["min_ct_per_unit_train"] = args.min_ct_per_unit_train
    # per-cell count thresholds (applied by the cell analyses after feature filtering)
    cd = prof.setdefault("cell_defaults", {})
    if args.cell_min_cell_count is not None:
        cd["min_cell_count"] = args.cell_min_cell_count
    if args.cell_min_feature_count is not None:
        cd["min_feature_count"] = args.cell_min_feature_count
    # single-molecule: default ON for pixel FICTURE, OFF for cell decode;
    # --always/--never force the same value for both.
    if args.always_single_molecule and args.never_single_molecule:
        sys.exit("ERROR: --always-single-molecule and --never-single-molecule are mutually exclusive.")
    if args.always_single_molecule:
        prof["_sm_pixel"], prof["_sm_cells"] = True, True
    elif args.never_single_molecule:
        prof["_sm_pixel"], prof["_sm_cells"] = False, False
    else:
        prof["_sm_pixel"], prof["_sm_cells"] = True, False
    fd.pop("single_molecule", None)   # now controlled by _sm_pixel/_sm_cells

    # Packaging knobs: a CLI flag wins over the profile's `cartload` block.
    if args.bin_count is not None:
        prof.setdefault("cartload", {})["bin_count"] = args.bin_count

    # Packaging without any factor analysis (tiling only; see cmd_ficture_analysis).
    prof["_no_ficture"] = args.no_ficture

    # 10x MEX export of the hexagon files (run_ficture2_multi --segment-10x). A top-level
    # key rather than a per-analysis one so it survives --no-ficture, which rebuilds the
    # ficture list; it applies to every analysis (targets are per width, and idempotent).
    # Normalized to None (off) or {"widths": "12,24" | None} (on; None = each analysis'
    # own --width). Config: `"segment_10x": true` or `{"widths": "..."}`; CLI wins.
    s10 = prof.get("segment_10x")
    if isinstance(s10, bool) or s10 is None:
        s10 = {} if s10 else None
    elif isinstance(s10, dict):
        bad = [k for k in s10 if k != "widths"]
        if bad and not args.allow_unknown_config_keys:
            sys.exit(f"ERROR: unrecognized key(s) in 'segment_10x': {', '.join(sorted(bad))}. "
                     f"Recognized keys: widths. Pass --allow-unknown-config-keys to ignore "
                     f"unknown keys instead of failing.")
        s10 = {"widths": str(s10["widths"])} if s10.get("widths") is not None else {}
    else:
        sys.exit("ERROR: 'segment_10x' must be true/false or an object like {\"widths\": \"12,24\"}.")
    if args.segment_10x:
        s10 = s10 if s10 is not None else {}
    if args.segment_width_10x is not None:
        if s10 is None:
            sys.exit("ERROR: --segment-width-10x requires --segment-10x (or 'segment_10x' in the config).")
        s10["widths"] = args.segment_width_10x
    prof["_segment_10x"] = s10

    # Tolerate corrupt histology images in the images stage (opt-in; see plan_images).
    prof["_skip_image_errors"] = args.skip_image_errors

    # Explicit cell-analysis selection for the generic platform (see plan_cell_analyses):
    # --colname-cell names the transcript's cell-id column and turns on cell_id-based
    # clustering; mex columns turn on mex-based clustering; neither -> cells skipped.
    prof["colname_cell"] = args.colname_cell or cfg.get("colname_cell")

    # Per-file column-name overrides (CLI). The inputs may be arbitrarily-named
    # individual files, so their column names are not fixed by the profile: the
    # transcript columns flow to sge_convert, the xy/boundary columns into their
    # roles (consumed by run_ficture2_multi_cells / spatula tsv-add-cell-id).
    tx_cols = {"x": args.colname_transcript_x, "y": args.colname_transcript_y,
               "feature": args.colname_transcript_feature, "count": args.colname_transcript_count}
    tx_cols = {k: v for k, v in tx_cols.items() if v}
    if tx_cols:
        prof.setdefault("ingest", {}).setdefault("csv_colnames", {}).update(tx_cols)
    # A transcript CSV that already carries a per-molecule cell-id column (rare for
    # MERFISH): naming it here carries the column through ingest to transcript column 5
    # (X, Y, gene, count, cell_id) and turns on cell analysis from that column, so
    # spatula tsv-add-cell-id is not needed. A cellxgene MEX, if also present, still
    # drives the clustering (mex2sptsv). Applies to every sge_convert sample in the run.
    if args.colname_transcript_cell:
        prof.setdefault("ingest", {})["csv_colname_cell"] = args.colname_transcript_cell
    # cell-id colnames use `is not None` because "" is meaningful (an unnamed
    # pandas-index first column, e.g. MERFISH cell metadata).
    xy_over = {}
    if args.colname_xy_cell is not None: xy_over["colname_cell"] = args.colname_xy_cell
    if args.colname_xy_x: xy_over["colname_x"] = args.colname_xy_x
    if args.colname_xy_y: xy_over["colname_y"] = args.colname_xy_y
    if xy_over:
        prof.setdefault("roles", {}).setdefault("xy", {}).update(xy_over)
    bnd_over = {}
    if args.colname_boundary_cell is not None: bnd_over["colname_cell"] = args.colname_boundary_cell
    if args.colname_boundary_x: bnd_over["colname_x"] = args.colname_boundary_x
    if args.colname_boundary_y: bnd_over["colname_y"] = args.colname_boundary_y
    if bnd_over:
        prof.setdefault("roles", {}).setdefault("boundaries", {}).update(bnd_over)

    # Apply ficture_defaults / cell_defaults to every analysis (per-entry keys win).
    prof["ficture"] = [deep_merge(prof["ficture_defaults"], a) for a in prof["ficture"]]
    prof["cell_analyses"] = [deep_merge(prof["cell_defaults"], ca) for ca in prof["cell_analyses"]]

    # Ready-made --list-* files (see explicit_lists). A run-wide default comes from the
    # top-level `cell_lists` block or the --list-* CLI flags (CLI wins) and applies to
    # every cell analysis; a per-analysis `lists` entry overrides it role by role.
    run_lists = dict(cfg.get("cell_lists") or {})
    for role, val in (("clusters", args.list_cluster), ("xy", args.list_xy),
                      ("boundaries", args.list_boundaries), ("mex", args.list_mex),
                      ("cell_tsv", args.list_cell_tsv)):
        if val:
            run_lists[role] = val
    # Validate the merged overrides in one place, before any planning: an unknown role
    # would otherwise be a --list-* flag that never appears, and a missing file would
    # surface as a make failure deep inside the cells stage. Paths are made absolute here
    # so they do not depend on where make is invoked from.
    for ca in prof["cell_analyses"]:
        # A display `name` is recorded into the catalogs by a generated make recipe, so
        # reject the characters that would break out of its quoting (or be eaten by make)
        # rather than emitting a recipe that fails or, worse, runs something else.
        name = ca.get("name")
        if name is not None:
            bad = [c for c in "'\"$`\\\n" if c in str(name)]
            if bad:
                sys.exit(f"ERROR: cell analysis '{ca.get('id', '?')}' has a name containing "
                         f"unsupported character(s) {' '.join(repr(c) for c in bad)}: {name!r}. "
                         f"Plain text (letters, digits, spaces, -_.,()/:) only.")
        if run_lists:
            ca["lists"] = {**run_lists, **explicit_lists(ca)}
        lists = explicit_lists(ca)
        if not isinstance(lists, dict):
            sys.exit(f"ERROR: 'lists' in cell analysis '{ca.get('id', '?')}' must be a map "
                     f"of role -> list file path (got {type(lists).__name__}).")
        for role, path in lists.items():
            if role not in ROLE_LIST_FLAG:
                sys.exit(f"ERROR: cell analysis '{ca.get('id', '?')}' names a list for role "
                         f"'{role}', which feeds no --list-* flag. Roles that take a list: "
                         f"{', '.join(sorted(ROLE_LIST_FLAG))}.")
            full = os.path.abspath(os.path.expanduser(path))
            if not os.path.exists(full):
                sys.exit(f"ERROR: file not found: {path} (cell analysis "
                         f"'{ca.get('id', '?')}', lists.{role})")
            lists[role] = full

    prof["out_dir"] = args.out_dir or cfg.get("out_dir")
    prof["out_root"] = args.out_root or cfg.get("out_root")
    prof.setdefault("resources", {})
    prof["resources"].setdefault("n_jobs", args.n_jobs)
    prof["resources"].setdefault("threads", args.threads)

    # Samples: --in-dir / explicit per-file flags (one), --samples sheet (many),
    # or JSON 'samples'. The single-sample CLI accepts either a directory (--in-dir,
    # roles auto-detected inside it) or individual files by explicit path
    # (--in-transcript raw CSV to ingest, --in-cell-xy, --in-cell-boundary) which
    # need not share a directory or use any standard filename.
    # The SAW binary (Stereo-seq): the GEF inputs are binary and only SAW can read
    # them. It is required only when a GEF is actually being converted, so the check
    # lives in cmds_stereoseq_ingest (a sample that supplies a pre-converted transcript
    # + cell_tsv skips SAW entirely, and should not need the binary).
    prof["saw"] = args.saw or cfg.get("saw")

    raw_samples = list(cfg.get("samples", []))
    if args.samples:
        raw_samples += read_sheet(args.samples)
    if args.in_transcript and args.raw_transcript:
        sys.exit("ERROR: --in-transcript (a pre-converted transcript that skips ingest) and "
                 "--raw-transcript (a raw CSV that goes through sge_convert) are mutually "
                 "exclusive; provide only one.")
    if (args.in_dir or args.in_prefix or args.in_transcript or args.raw_transcript
            or args.in_cell_xy or args.in_cell_boundary or args.in_cellxgene):
        s = {}
        if args.in_dir: s["in_dir"] = args.in_dir
        if args.in_prefix: s["in_prefix"] = args.in_prefix
        if args.id: s["id"] = args.id
        # --in-transcript is an already-converted pixel transcript (the transcripts.unsorted.tsv.gz
        # layout: X, Y, gene, count[, cell_id]); it skips ingest entirely via the `transcript`
        # role, on every platform. A raw platform CSV that still needs sge_convert is passed with
        # --raw-transcript instead (raw_transcript). Stereo-seq has no raw-CSV path (SAW expands
        # the binary GEF), so only --in-transcript applies there.
        if args.in_transcript: s["transcript"] = args.in_transcript
        if args.raw_transcript: s["raw_transcript"] = args.raw_transcript
        if args.in_cell_xy: s["xy"] = args.in_cell_xy
        if args.in_cell_boundary: s["boundaries"] = args.in_cell_boundary
        if args.in_cellxgene: s["cellxgene"] = args.in_cellxgene
        raw_samples.append(s)
    if not raw_samples:
        sys.exit("ERROR: no samples. Use --in-dir, --in-prefix, the "
                 "--in-transcript/--raw-transcript/--in-cell-xy/--in-cell-boundary file flags, "
                 "--samples <sheet>, or a 'samples' config block.")
    prof["_raw_samples"] = raw_samples

    # Per-image specs: the --images TSV (one row per image) plus any repeatable
    # --image CLI entries (each a comma-separated key=value row, for the single-sample
    # case). Both attach to samples by the `sample` column during resolve_sample.
    prof["_image_rows"] = read_sheet(args.images) if args.images else []
    prof["_image_rows"] += [parse_image_arg(s) for s in (args.image or [])]

    if not prof["out_dir"] and not prof["out_root"]:
        sys.exit("ERROR: --out-dir (or --out-root for independent per-sample models) is required.")
    return prof


def _abs_in_dir(path, in_dir):
    """Resolve a possibly-relative sheet path against the sample's in_dir."""
    if in_dir and not os.path.isabs(path) and not os.path.exists(path):
        cand = os.path.join(in_dir, path)
        return cand if os.path.exists(cand) else path
    return path


def sheet_image_specs(raw, in_dir, cfg):
    """Turn friendly single-file image columns (dapi/hne) into image specs the
    image stage understands. `dapi` is a single-channel colorized layer; `hne` is
    an RGB layer (unless the profile already provides a dedicated hne handler,
    e.g. Visium HD, which needs scale metadata and is left to that path)."""
    specs = []
    if raw.get("dapi"):
        src = _abs_in_dir(raw["dapi"], in_dir)
        is_ome = src.lower().endswith((".ome.tif", ".ome.tiff"))
        specs.append({"id": "dapi", "source": src, "kind": "single",
                      "color": "0F73E6", "convert": "ome2png" if is_ome else "none"})
    if raw.get("hne") and not cfg.get("hne"):
        specs.append({"id": "hne", "source": _abs_in_dir(raw["hne"], in_dir), "kind": "rgb"})
    return specs


def image_row_to_spec(row, cfg, in_dir):
    """Turn one --images TSV row into an image spec. The `type` sets the default
    color and kind (from IMAGE_TYPES); `id` defaults to `type`; `source`/`src` is
    the image path (.tif or .png). A platform's transform column (declared by the
    profile's `image_transform`, e.g. MERFISH `merfish_csv` -> --micron2pixel-csv)
    is carried as transform_flag/transform_path. Blank cells were already dropped
    by read_sheet, so any present key is meaningful."""
    itype = row.get("type")
    if not itype:
        sys.exit(f"ERROR: --images row is missing the required 'type' column: {row}")
    reg = IMAGE_TYPES.get(itype, {})
    kind = row.get("kind") or reg.get("kind", "single")
    src = row.get("source") or row.get("src") or row.get("tif")
    if not src:
        sys.exit(f"ERROR: --images row for type '{itype}' is missing an image path "
                 f"('source'/'src' column): {row}")
    spec = {"id": row.get("id") or row.get("img_id") or itype, "type": itype,
            "kind": kind, "source": _abs_in_dir(src, in_dir)}
    if kind != "rgb":
        # colorized single-channel: explicit column wins, else the type's default.
        color = row.get("color") or reg.get("color")
        if not color:
            sys.exit(f"ERROR: image type '{itype}' has no default color and the --images row "
                     f"gives none. Add a 'color' column (hex) or register the type in IMAGE_TYPES.")
        spec["color"] = color
    # Platform-specific geometric transform: the profile declares its column name
    # and import_image flag; a generic `transform` key is also accepted so the same
    # spelling works across platforms (the platform-named column wins).
    tr = cfg.get("image_transform")
    if tr:
        val = row.get(tr["column"]) or row.get("transform")
        if val:
            spec["transform_flag"] = tr["flag"]
            spec["transform_path"] = _abs_in_dir(val, in_dir)
    for k in ("shrink_factor", "high_memory", "use_middle_page", "convert", "um_per_pixel",
              "georef_plain", "georeferenced", "georef_detect",
              "swapxy", "rotate", "flip_vertical", "flip_horizontal",
              "rescale", "rescale_range", "rescale_min", "rescale_max"):
        if row.get(k):
            spec[k] = row[k]
    # A plain .png needs no OME->PNG conversion; a .tif/.ome.tif does (default).
    spec.setdefault("convert", "none" if spec["source"].lower().endswith(".png") else "ome2png")
    return spec


# OME-TIFFs known to be multi-page z-stacks: Xenium's morphology.ome.tif is the 3D DAPI
# stack (the 2D projections are morphology_focus*.ome.tif). import_image refuses a
# multi-page file without a page, so these default to --use-middle-page wherever the
# image comes from (profile match, sample sheet `dapi` column, --images row, --image,
# config). Matched by suffix so renamed copies (e.g. GEO's GSM123_rep1_morphology.ome.tif)
# qualify too, while morphology_focus.ome.tif does not (its stem ends in "_focus"). A
# `use_middle_page` key/column on the image (true/false) overrides the default.
ZSTACK_OME_SUFFIXES = ("morphology.ome.tif", "morphology.ome.tiff")


def default_use_middle_page(src):
    return os.path.basename(src).lower().endswith(ZSTACK_OME_SUFFIXES)


def _truthy(v):
    """Coerce a sheet string ('true'/'1'/…) or a JSON bool to a boolean."""
    return str(v).strip().lower() in ("1", "true", "yes", "t", "y") if v is not None else False


def _image_orient_flags(settings, idef):
    """Rotate/flip flags for image_png2pmtiles, from a per-image spec (else image_defaults).

    `swapxy` is the friendly spelling of an x/y transpose — a reflection about the main
    diagonal, which is what a histology image needs when its axes are stored swapped relative
    to the transcript coordinates. The rotate/flip vocabulary already expresses it as a 90
    degree rotation followed by a vertical flip (orient_helper's orient2axisorder maps that
    pair to the gdalwarp axisswap order "2,1"), so it is expanded here rather than carried
    any further: downstream there stays exactly one orientation mechanism.
    """
    def get(k):
        return settings.get(k, idef.get(k))

    rotate = get("rotate")
    vflip = _truthy(get("flip_vertical"))
    hflip = _truthy(get("flip_horizontal"))
    if _truthy(get("swapxy")):
        if rotate or vflip or hflip:
            sys.exit("ERROR: image 'swapxy' is shorthand for rotate=90 plus flip_vertical, so it "
                     "cannot be combined with an explicit rotate/flip_vertical/flip_horizontal. "
                     "Drop swapxy and spell the orientation out instead.")
        rotate, vflip = "90", True
    parts = []
    if rotate:
        parts.append(f"--rotate {rotate}")
    if vflip:
        parts.append("--flip-vertical")
    if hflip:
        parts.append("--flip-horizontal")
    return (" " + " ".join(parts)) if parts else ""


def parse_image_arg(s):
    """Parse one --image value (comma-separated key=value pairs) into an image row,
    the CLI equivalent of a single --images TSV row (same field vocabulary)."""
    row = {}
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        k, sep, v = tok.partition("=")
        if not sep:
            sys.exit(f"ERROR: --image expects comma-separated key=value pairs (e.g. "
                     f"type=dapi,source=/p/dapi.tif,transform=/p/x.csv); got '{tok}' in '{s}'.")
        row[k.strip()] = v.strip()
    return row


def resolve_sample(raw, cfg):
    """Resolve one sample's id, out_dir, and input roles."""
    n = len(cfg["_raw_samples"])
    if cfg.get("out_dir") and n == 1 and "id" not in raw:
        # A lone --out-dir sample takes the generic id rather than the out-dir basename,
        # which is typically a long collection-style name and uninformative as a sample
        # id. The out-dir name is not lost: the packaged directory and catalog id become
        # <out_dir_basename>-rep1. An --out-root run keeps naming each sample after its
        # input instead, since those ids distinguish samples from one another.
        sid = DEFAULT_SAMPLE_ID
    else:
        sid = raw.get("id")
        if not sid and raw.get("in_dir"):
            sid = os.path.basename(os.path.normpath(raw["in_dir"]))
        # A prefix-addressed sample (Stereo-seq) has no directory of its own; its
        # trailing path component is the chip/sample name (e.g. .../C04687E314).
        if not sid and raw.get("in_prefix"):
            sid = os.path.basename(raw["in_prefix"])
    if not sid:
        sys.exit(f"ERROR: cannot determine id for sample {raw}")

    out_dir = cfg["out_dir"] if cfg.get("out_dir") else os.path.join(cfg["out_root"], sid)
    in_dir = raw.get("in_dir")
    in_prefix = raw.get("in_prefix")

    def sheet_value(role):
        # explicit role column wins, else a friendly alias column (e.g. tsv -> transcript)
        if raw.get(role):
            return raw[role]
        for alias, target in ROLE_ALIASES.items():
            if target == role and raw.get(alias):
                return raw[alias]
        return None

    # Resolve roles: explicit column/alias/JSON value wins; else auto-detect, either as
    # a named `file` inside in_dir or as a `suffix` appended to in_prefix.
    roles = {}
    role_specs = cfg.get("roles", {})
    for role in ROLE_KEYS:
        val = sheet_value(role)
        if val:
            roles[role] = _abs_in_dir(val, in_dir)
            continue
        spec = role_specs.get(role, {})
        cand = None
        if spec.get("file") and in_dir:
            cand = os.path.join(in_dir, spec["file"])
        elif spec.get("suffix") and in_prefix:
            cand = in_prefix + spec["suffix"]
        if cand and os.path.exists(cand):
            roles[role] = cand

    # The mex role may instead be given as an explicit bcd/ftr/mtx triple (a dict
    # value, written as a 4-column --mex-list line) when a single directory does
    # not apply. An explicit `mex`/`mex_dir` column (handled above) takes precedence.
    if "mex" not in roles and all(raw.get(k) for k in ("mex_bcd", "mex_ftr", "mex_mtx")):
        roles["mex"] = {k: _abs_in_dir(raw[f"mex_{k}"], in_dir) for k in ("bcd", "ftr", "mtx")}

    # Images from the --images TSV whose `sample` targets this sample (an exact id,
    # or '*'/blank meaning every sample).
    tsv_images = [image_row_to_spec(row, cfg, in_dir)
                  for row in cfg.get("_image_rows", [])
                  if (row.get("sample") or row.get("sample_id") or "*") in ("*", sid)]

    return {
        "id": sid, "in_dir": in_dir, "in_prefix": in_prefix, "out_dir": out_dir,
        "roles": roles,
        # The transcript already carries a cell-id column (declared via
        # ingest.csv_colname_cell / --colname-transcript-cell): use it directly for cell
        # analysis and skip tsv-add-cell-id. See build_config / add_targets.
        "has_tx_cell_id": bool(cfg.get("ingest", {}).get("csv_colname_cell")),
        # raw transcript CSV to ingest (sge_convert --in-csv); distinct from the
        # `transcript` role, which is an already-ingested TSV that skips ingest.
        "raw_transcript": raw.get("raw_transcript"),
        "hne": raw.get("hne"),
        "images": list(raw.get("images", [])) + tsv_images + sheet_image_specs(raw, in_dir, cfg),
    }


# ---------------------------------------------------------------------------
# Command builders
# ---------------------------------------------------------------------------

def ingest_feature_filter_flags(cfg):
    """Feature filters for the ingest step (sge_convert / reformat_cosmx / the Stereo-seq
    cell-bin conversion), i.e. the step that writes the transcript TSV.

    A feature dropped here never reaches any later stage — it is absent from the transcript
    TSV, the feature list, the packaged tiles and the browser's gene list. That is the whole
    difference from the factor-analysis filters (see feature_filter_flags), which leave the
    data complete and narrow only the model fitting; the two are independent.

    With no ingest filter set, an empty exclude regex is passed explicitly so that
    sge_convert's per-platform default (e.g. Xenium's negative-probe pattern) does not
    silently reintroduce filtering.
    """
    keys = {"--include-feature-list": "ingest_include_feature_list",
            "--exclude-feature-list": "ingest_exclude_feature_list",
            "--include-feature-regex": "ingest_include_feature_regex",
            "--exclude-feature-regex": "ingest_exclude_feature_regex"}
    parts = [f'{flag} "{cfg[key]}"' if "regex" in flag else f"{flag} {cfg[key]}"
             for flag, key in keys.items() if cfg.get(key)]
    if not cfg.get("ingest_exclude_feature_regex"):
        parts.append('--exclude-feature-regex ""')
    return parts


# Sample inputs whose consumers read CSV/TSV only: the raw transcript goes to sge_convert
# --in-csv, and the cell roles to spatula / the per-sample importers' --csv-* readers. A
# `.parquet` given for one of them (Xenium ships transcripts.parquet, cells.parquet and
# cell_boundaries.parquet beside the .csv.gz files) is converted to .csv.gz first.
PARQUET_CONVERTIBLE_ROLES = ("xy", "boundaries", "clusters")


def _is_parquet(path):
    return isinstance(path, str) and path.lower().endswith(".parquet")


def plan_parquet_conversions(s, sge_dir, cfg):
    """Convert this sample's .parquet inputs (raw_transcript and the CSV-only roles) to
    .csv.gz under <sge_dir>/parquet2csv/ with `cartloader parquet_to_csv_rapid`, rewriting
    the sample's paths to the converted files. Returns the shell commands (empty when
    nothing needs converting); the caller makes them an ingest-stage target that every
    consumer depends on. A failed or empty conversion aborts the run with an explicit
    error rather than letting a downstream reader choke on the parquet."""
    res = cfg["resources"]
    out_dir = os.path.join(sge_dir, "parquet2csv")
    cmds = []
    items = [("raw_transcript", s.get("raw_transcript"))] + \
            [(r, s["roles"].get(r)) for r in PARQUET_CONVERTIBLE_ROLES]
    for role, path in items:
        if not _is_parquet(path):
            continue
        out = os.path.join(out_dir, os.path.basename(path)[:-len(".parquet")] + ".csv.gz")
        print(f"NOTE: sample {s['id']}: '{role}' is a parquet file ({path}); it will be converted "
              f"to {out} with 'cartloader parquet_to_csv_rapid' before it is used.", file=sys.stderr)
        cmds.append(f"echo 'Converting parquet to csv.gz ({role}): {path} -> {out}'")
        cmds.append(f"cartloader parquet_to_csv_rapid --in-parquet {path} --out-csv-gz {out} "
                    f"--threads {res['threads']} || {{ rm -f {out}; "
                    f"echo 'ERROR: parquet -> csv.gz conversion failed for {path} ({role} of sample {s['id']})' >&2; exit 1; }}")
        cmds.append(f"[ -s {out} ] || {{ echo 'ERROR: parquet -> csv.gz conversion of {path} produced no output ({out})' >&2; exit 1; }}")
        if role == "raw_transcript":
            s["raw_transcript"] = out
        else:
            s["roles"][role] = out
    return cmds


def cmd_sge_convert(cfg, sge_dir, s):
    ing = cfg.get("ingest", {})
    res = cfg["resources"]
    in_dir = s.get("in_dir")
    parts = ["cartloader", "sge_convert",
             f"--platform {ing.get('sge_platform', cfg['platform'])}",
             f"--out-dir {sge_dir}", f"--n-jobs {res['n_jobs']}",
             f"--pigz-threads {res['threads']}", "--gzip pigz"]
    parts.extend(ingest_feature_filter_flags(cfg))
    # An explicit raw transcript CSV (single-sample --in-transcript) is ingested
    # directly and takes precedence over in_dir autodetection.
    raw_tx = s.get("raw_transcript")
    if raw_tx:
        parts.append(f"{ing.get('raw_input_flag', '--in-csv')} {raw_tx}")
    # inputs taken from resolved sample roles, e.g. illumina's --in-mex from the
    # `mex` role (a mex_dir sample-sheet column) since its layout is not standardized.
    for flag, role in ing.get("input_roles", {}).items():
        val = s["roles"].get(role)
        if not val:
            sys.exit(f"ERROR: {cfg['platform']} ingest needs the '{role}' role for {flag}; "
                     f"provide it as a sample-sheet column (e.g. {role}_dir/{role}) or in_dir.")
        if isinstance(val, dict):
            sys.exit(f"ERROR: {cfg['platform']} ingest expects a single path for the '{role}' role "
                     f"for {flag}, not a bcd/ftr/mtx triple.")
        parts.append(f"{flag} {val}")
    # A platform whose whole input IS the directory (SeqScope: --in-dir is the MEX
    # directory holding barcodes/features/matrix, with no standard parent layout).
    # Skipped when the sample names a raw transcript instead: that is a second, complete
    # input route (SeqScope's raw per-molecule TSV via --raw-transcript), not a supplement
    # to the directory, so requiring --in-dir as well would be wrong.
    if ing.get("in_dir_flag") and not raw_tx:
        if not in_dir:
            sys.exit(f"ERROR: {cfg['platform']} ingest reads its input directory directly; "
                     f"provide it with --in-dir <dir> (or an 'in_dir' sample-sheet column).")
        parts.append(f"{ing['in_dir_flag']} {in_dir}")
    if not raw_tx and "autodetect" in ing:
        if not in_dir:
            names = ", ".join(c["file"] for c in ing["autodetect"])
            sys.exit(f"ERROR: {cfg['platform']} ingest needs a transcript input. Provide it explicitly "
                     f"with --in-transcript <raw CSV>, or point --in-dir at a directory containing one "
                     f"of: {names}.")
        chosen = None
        for cand in ing["autodetect"]:
            p = os.path.join(in_dir, cand["file"])
            if os.path.exists(p):
                chosen = (cand["flag"], p); break
        if chosen is None:
            cand = ing["autodetect"][0]
            chosen = (cand["flag"], os.path.join(in_dir, cand["file"]))
        parts.append(f"{chosen[0]} {chosen[1]}")
    for flag, rel in ing.get("inputs", {}).items():
        if in_dir:
            parts.append(f"{flag} {os.path.join(in_dir, rel)}")
    # Per-input-column name overrides (the input files may use non-standard column names).
    for key, val in ing.get("csv_colnames", {}).items():
        flag = CSV_COLNAME_FLAGS.get(key)
        if flag and val:
            parts.append(f"{flag} {val}")
    # Columns to carry through beyond X/Y/gene/count. A declared transcript cell-id
    # column is kept first so it lands at column 5 (X, Y, gene, count, cell_id) — the
    # position the cell decode reads (run_ficture2_multi_cells --colidx-cell-id 5).
    others = list(ing.get("csv_colnames_others", []))
    if ing.get("csv_colname_cell"):
        others = [ing["csv_colname_cell"]] + [c for c in others if c != ing["csv_colname_cell"]]
    if others:
        parts.append("--csv-colnames-others " + " ".join(others))
    # Random jitter (um) added to X and Y as the transcript TSV is written. Both
    # sge_convert routes accept it -- spatula convert-sge for MEX input, sge_format_generic
    # for CSV input -- so every platform whose ingest method is sge_convert supports it
    # (build_config rejects the platforms whose ingest is something else).
    if ing.get("jitter_xy"):
        parts.append(f"--jitter-xy {ing['jitter_xy']}")
    # Coordinate units per um in the raw input. Only passed when the profile/config/CLI
    # names one: left out, sge_convert applies its own per-platform default (e.g. the
    # Illumina barcode-unit detection), which is the better answer than a blanket 1.
    if ing.get("units_per_um") is not None:
        parts.append(f"--units-per-um {ing['units_per_um']}")
    parts.extend(ing.get("extra_flags", []))
    return " ".join(p for p in parts if p)


def cmd_reformat_cosmx(cfg, sge_dir, sid, in_dir):
    """Custom CosMx ingest: reformat raw *_tx_file/*_metadata_file/*-polygons CSVs into
    the transcript TSV + cell metadata (xy) + polygon (boundaries) files. Input patterns
    are globbed against the sample's in_dir at plan time so a missing file fails early."""
    ing = cfg.get("ingest", {})
    if not in_dir:
        sys.exit("ERROR: cosmx ingest requires a sample 'in_dir'.")
    parts = ["cartloader", "reformat_cosmx"]
    for flag, pattern in ing.get("inputs", {}).items():
        if not pattern:      # a config may null out an optional input to drop it
            continue
        matches = sorted(glob.glob(os.path.join(in_dir, pattern)))
        if not matches:
            avail = ", ".join(sorted(os.listdir(in_dir))) if os.path.isdir(in_dir) else "(not a directory)"
            sys.exit(f"ERROR: cosmx ingest: no file matching '{pattern}' for {flag} in {in_dir}.\n"
                     f"       Files present: {avail}\n"
                     f"       Override the glob via the config 'ingest.inputs' block (e.g. "
                     f"{{\"ingest\": {{\"inputs\": {{\"{flag}\": \"*your_suffix.csv.gz\"}}}}}}).")
        if len(matches) > 1:
            names = ", ".join(os.path.basename(m) for m in matches)
            sys.exit(f"ERROR: cosmx ingest: pattern '{pattern}' for {flag} matched multiple files in "
                     f"{in_dir}: {names}.\n       Make the pattern more specific in 'ingest.inputs'.")
        parts.append(f"{flag} {matches[0]}")
    parts.append(f"--out {os.path.join(sge_dir, sid)}")
    # Same ingest-level feature filters as sge_convert; this step writes the transcript TSV
    # for CosMx, so a filtered feature never enters the run. (No platform default to
    # suppress here, so an empty exclude regex is simply a no-op.)
    parts.extend(ingest_feature_filter_flags(cfg))
    parts.extend(ing.get("extra_flags", []))
    return " ".join(parts)


def cmds_stereoseq_ingest(cfg, sge_dir, s):
    """Stereo-seq ingest: SAW expands the binary GEFs into text GEMs, which are then
    read by the ordinary pixel path. Returns (commands, cell_tsv or None).

    Two GEMs are produced from the same run:

      * `{prefix}.tissue.gef` -> a bin1 GEM (one row per 0.5um position per gene). This
        is the pixel-level transcript; sge_convert ingests it with --units-per-um 2 so
        the 0.5um grid lands in um.
      * `{prefix}.cellbin.gef` -> a cell-bin GEM, the subset of MIDs that segmentation
        placed inside a cell. Its coordinates cannot be mapped back onto the bin1 rows,
        so it is NOT merged into the transcript: it becomes a standalone pixel TSV that
        the cells stage feeds to run_ficture2_multi_cells as a --tsv-list entry.

    The GEMs are tens of GB of text, so each is deleted as soon as it has been read.
    A re-run therefore re-invokes SAW; that is the tradeoff for not parking the
    intermediates on disk for the life of the output tree.
    """
    ing = cfg.get("ingest", {})
    saw = cfg["saw"]
    gef = s["roles"].get("gef")
    if not gef:
        sys.exit(f"ERROR: stereo-seq ingest for sample '{s['id']}' found no "
                 f"{ing.get('gef_suffix', '.tissue.gef')} input. Point --in-prefix at the "
                 f"sample prefix (e.g. --in-prefix /data/C04687E314), or give the path in a "
                 f"'gef' sample-sheet column. To skip SAW entirely, supply a pre-converted "
                 f"'transcript' TSV (and a 'cell_tsv' for cell analysis) instead.")
    # SAW is needed only here, where a GEF is actually converted; a sample that
    # supplies a pre-converted transcript never reaches this function.
    if not saw:
        sys.exit(f"ERROR: converting the .gef inputs for sample '{s['id']}' needs the SAW "
                 f"binary; provide it via --saw <path to saw>. To skip SAW, supply a "
                 f"pre-converted 'transcript' TSV (and 'cell_tsv') instead of a GEF.")

    bin1_gem = os.path.join(sge_dir, "bin1.gem")
    cmds = [f"{saw} convert gef2gem --bin-size 1 --gef {gef} --gem {bin1_gem}",
            cmd_sge_convert(cfg, sge_dir, {**s, "raw_transcript": bin1_gem}),
            f"rm -f {bin1_gem}"]

    cellbin_gef = s["roles"].get("cellbin_gef")
    if not cellbin_gef:
        return cmds, None

    # The cell clusters are projected onto the model trained from the bin1 features, so
    # --check-features fails the run on a feature-naming mismatch between the two GEMs
    # rather than letting it surface as a decode over near-empty cells.
    cellbin_gem = os.path.join(sge_dir, "cellbin.gem")
    cell_tsv = os.path.join(sge_dir, "cellbin.tsv")
    # Must match sge_convert's --out-feature default (singular "feature.clean.tsv.gz").
    feature_f = os.path.join(sge_dir, ing.get("feature_file", "feature.clean.tsv.gz"))
    conv = ["cartloader", "convert_stereoseq_cellbin",
            f"--in-gem {cellbin_gem}", f"--out {cell_tsv}",
            f"--units-per-um {ing.get('units_per_um', 2)}",
            f"--check-features {feature_f}"]
    # The cell-bin GEM is a second ingest product of the same run, so it takes the same
    # ingest filters as bin1: without them the cells would carry features the pixel-level
    # transcript (and the model trained on it) no longer has.
    conv.extend(ingest_feature_filter_flags(cfg))
    # Carry the same feature/count column overrides sge_convert used for bin1 so the
    # cellbin GEM is read with matching columns (e.g. protein GEFs use MIDCount, not
    # the ExonCount default).
    feature_col = ing.get("csv_colnames", {}).get("feature")
    if feature_col:
        conv.append(f"--colname-feature {feature_col}")
    count_col = ing.get("csv_colnames", {}).get("count")
    if count_col:
        conv.append(f"--colname-count {count_col}")
    cmds += [f"{saw} convert gef2gem --cellbin-gef {cellbin_gef} --gef {gef} --cellbin-gem {cellbin_gem}",
             " ".join(conv),
             f"rm -f {cellbin_gem}"]
    return cmds, cell_tsv


def cmd_tsv_add_cell_id(cfg, in_tsv, out_tsv, boundaries):
    """Append a cell_id column to the ingested transcript by assigning each transcript
    to the cell-boundary polygon that contains it (point-in-polygon). MERFISH
    transcripts carry no per-transcript cell assignment (unlike CosMx, whose reformat
    writes one), so this reproduces the CosMx transcript contract (X, Y, gene, count,
    cell_id) that the standard cell decode (run_ficture2_multi_cells) consumes. The
    boundary CSV's cell-id/vertex columns come from the profile's roles.boundaries."""
    ing = cfg.get("ingest", {})
    ac = ing.get("assign_cell_id")
    ac = ac if isinstance(ac, dict) else {}
    bnd = cfg.get("roles", {}).get("boundaries", {})
    res = cfg["resources"]
    parts = [SPATULA_BIN, "tsv-add-cell-id",
             f"--tsv {in_tsv}", f"--out {out_tsv}",
             f"--boundaries-csv {boundaries}",
             f"--colname-x {ac.get('colname_x', 'X')}",
             f"--colname-y {ac.get('colname_y', 'Y')}",
             f"--csv-colname-cell-id {bnd.get('colname_cell', 'cell_id')}",
             f"--csv-colname-x {bnd.get('colname_x', 'vertex_x')}",
             f"--csv-colname-y {bnd.get('colname_y', 'vertex_y')}",
             f"--threads {res['threads']}"]
    if ac.get("expand_um"):
        # assign a transcript outside every polygon to the nearest cell within this
        # distance (µm); off by default (transcripts outside all cells stay UNASSIGNED).
        parts.append(f"--expand-um {ac['expand_um']}")
    return " ".join(parts)


def cmd_convert_cellxgene(cfg, csv_path, mex_dir):
    """Convert a cell-by-gene matrix CSV (e.g. MERSCOPE cell_by_gene.csv) into a MEX
    directory (barcodes/features/matrix). The produced dir becomes the sample's `mex`
    role, which drives the cell decode's clustering (mex2sptsv) — independent of, and
    taking precedence over, the transcript/boundary cell-count path."""
    res = cfg["resources"]
    ing = cfg.get("ingest", {})
    parts = ["cartloader", "convert_cellxgene",
             f"--csv {csv_path}", f"--out-dir {mex_dir}",
             "--gzip pigz", f"--threads {res['threads']}"]
    if ing.get("cellxgene_blank_prefix"):
        parts.append(f"--blank-prefix {ing['cellxgene_blank_prefix']}")
    return " ".join(parts)


def feature_filter_flags(cfg):
    """FICTURE feature restriction for the pixel factor analysis (run_ficture2_multi).

    run_ficture2_multi resolves these into a feature list and passes it to `lda4hex
    --features`, so they restrict the *model* only: the tiled transcript, the packaged
    feature list, and the hexagon counts keep every gene. A list and a regex combine.

    Deliberately NOT forwarded to the cell path. There the same FICTURE restriction reaches
    clustering through projection onto this (already --features-restricted) pixel model,
    while the cell *counts* get the ingest filter instead (see cells_count_filter_flags), so
    FICTURE-excluded and non-HVG genes stay in the cell pseudobulk and DE.
    """
    parts = []
    if cfg.get("include_feature_list"):
        parts.append(f"--include-feature-list {cfg['include_feature_list']}")
    if cfg.get("exclude_feature_list"):
        parts.append(f"--exclude-feature-list {cfg['exclude_feature_list']}")
    if cfg.get("exclude_feature_regex"):
        parts.append(f"--exclude-feature-regex \"{cfg['exclude_feature_regex']}\"")
    return parts


def cells_count_filter_flags(cfg):
    """Ingest-level feature filters for the cell count matrices (run_ficture2_multi_cells).

    These flags reach `mex2sptsv` / `pixel2sptsv`, i.e. the per-cell COUNTS, so they must be
    the *ingest* exclusions (technical artifacts), never the FICTURE ones. Two reasons:
      * a MEX-derived cell (cell_by_gene, or an external mex role) never passes through
        sge_convert, so this is the only place its Neg/BLANK/etc. are removed — closing the
        gap where the transcript was ingest-filtered but the MEX was not;
      * omitting the FICTURE exclude/include keeps those genes in the cell counts, so they
        reappear in the pseudobulk (sptsv2model, unrestricted) and DE, per design. The
        FICTURE restriction still governs clustering, via the projected pixel model.
    """
    keys = {"--include-feature-list": "ingest_include_feature_list",
            "--exclude-feature-list": "ingest_exclude_feature_list",
            "--include-feature-regex": "ingest_include_feature_regex",
            "--exclude-feature-regex": "ingest_exclude_feature_regex"}
    parts = []
    for flag, key in keys.items():
        if cfg.get(key):
            parts.append(f'{flag} "{cfg[key]}"' if "regex" in flag else f"{flag} {cfg[key]}")
    return parts


def cmd_ficture_analysis(a, in_list, fic_dir, cfg):
    res = cfg["resources"]
    parts = ["cartloader", "run_ficture2_multi",
             f"--in-list {in_list}", f"--out-dir {fic_dir}",
             f"--width {a['width']}", f"--threads {res['threads']}",
             f"--n-jobs {res['n_jobs']}", "--gzip pigz"]
    if a.get("mode") == "prepare":
        # Tiling only: no model is trained or projected, and the manifest it writes
        # carries just the tiled transcript for packaging (--no-ficture).
        parts.append("--prepare-only")
    elif a.get("mode") == "project":
        parts += [f"--pretrained-model {a['model']}", f"--model-id {a['id']}"]
    else:
        parts.append(f"--n-factor {a['n_factor']}")
    if a.get("min_ct_per_unit_hexagon") is not None:
        parts.append(f"--min-ct-per-unit-hexagon {a['min_ct_per_unit_hexagon']}")
    if a.get("min_ct_per_unit_train") is not None:
        parts.append(f"--min-ct-per-unit-train {a['min_ct_per_unit_train']}")
    if a.get("decode_scale"):
        parts.append(f"--decode-scale {a['decode_scale']}")
    if cfg.get("_sm_pixel"):   # single-molecule for pixel FICTURE (default ON)
        parts.append("--single-molecule")
    s10 = cfg.get("_segment_10x")
    if s10 is not None:        # 10x MEX export of the hexagon files (see build_config)
        parts.append("--segment-10x")
        if s10.get("widths"):
            parts.append(f"--segment-width-10x {s10['widths']}")
    parts.extend(feature_filter_flags(cfg))
    return " ".join(parts)


def cmd_cells(ca, list_files, fic_dir, model_path, cfg):
    res = cfg["resources"]
    parts = ["cartloader", "run_ficture2_multi_cells", "--all",
             f"--out-prefix {ca['id']}", f"--out-dir {fic_dir}",
             f"--threads {res['threads']}", f"--n-jobs {res['n_jobs']}",
             f"--pretrained-model {model_path}", "--gzip pigz"]
    xy_cfg = cfg.get("roles", {}).get("xy", {})
    for role, path in list_files.items():
        parts.append(f"{ROLE_LIST_FLAG[role]} {path}")
        if role == "xy":
            parts.append(f"--xy-colname-x {xy_cfg.get('colname_x', 'X')}")
            parts.append(f"--xy-colname-y {xy_cfg.get('colname_y', 'Y')}")
            # Cell-id column of the xy metadata file: a role-level colname_cell wins
            # (may be "" for an unnamed pandas-index first column, e.g. MERFISH cell
            # metadata), else the global --colname-cell, else the tool default (cell_id).
            cell_col = xy_cfg.get("colname_cell", cfg.get("colname_cell"))
            if cell_col is not None:
                parts.append(f"--xy-colname-cell-id '{cell_col}'")
        if role == "boundaries":
            # A boundaries `format` opts the decode into deriving per-cell centroids from
            # the polygons (for samples that have boundaries but no cell XY, e.g. default
            # Visium HD segmentation). Without it, boundaries stay pass-through.
            bnd = cfg.get("roles", {}).get("boundaries", {})
            if bnd.get("format"):
                parts.append(f"--boundaries-format {bnd['format']}")
                if bnd.get("cell_id_format"):
                    parts.append(f"--boundaries-cell-id-format '{bnd['cell_id_format']}'")
                if bnd.get("cell_id_prop"):
                    parts.append(f"--boundaries-cell-id-prop {bnd['cell_id_prop']}")
                if bnd.get("units_key"):
                    parts.append(f"--boundaries-units-key {bnd['units_key']}")
    if cfg.get("_sm_cells"):   # single-molecule for cell decode (default OFF)
        parts.append("--single-molecule")
    # Count thresholds are applied after feature filtering, so a run restricted to a small
    # panel usually wants them lowered (see cell_defaults / --cell-min-*-count).
    if ca.get("min_cell_count") is not None:
        parts.append(f"--min-cell-count {ca['min_cell_count']}")
    if ca.get("min_feature_count") is not None:
        parts.append(f"--min-feature-count {ca['min_feature_count']}")
    # Cell counts get the ingest filter, NOT the FICTURE filter: FICTURE-excluded / non-HVG
    # genes must survive into the cell pseudobulk and DE (the FICTURE restriction reaches
    # clustering via the projected pixel model instead). See cells_count_filter_flags.
    parts.extend(cells_count_filter_flags(cfg))
    # Escape hatch for run_ficture2_multi_cells flags this builder does not model. Needed
    # in practice by `lists`-supplied clusters: cluster ids are assumed 1-based and are
    # decremented, so externally assigned 0-based labels need --zero-based-clust-id here.
    parts.extend(ca.get("extra_flags", []))
    return " ".join(parts)


def cmd_cartload(fic_sample_dir, cart_dir, sid, cfg, cell_params):
    cart = cfg.get("cartload", {})
    res = cfg["resources"]
    parts = ["cartloader", "run_cartload2",
             f"--fic-dir {fic_sample_dir}", f"--out-dir {cart_dir}",
             f"--id {sid}", f"--n-jobs {res['n_jobs']}", f"--threads {res['threads']}",
             "--colname-count count", "--gzip pigz"]
    if cell_params:
        parts.append("--in-cell-params " + " ".join(cell_params))
    if cart.get("use_pmpoint"):
        parts.append("--use-pmpoint")
    if cart.get("sge_scale"):
        parts.append(f"--sge-scale {cart['sge_scale']}")
    if cart.get("bin_count") is not None:
        parts.append(f"--bin-count {cart['bin_count']}")
    return " ".join(parts)


def cmd_record_alias(cart_dir, catalog_path, oid, alias_path):
    """Deploy a companion alias (manual factor labels) beside a projection model:
    copy it into `cart_dir` as `<oid>-alias.tsv` and record it under the factor's
    `alias` key in `catalog_path`. Distinct from the AI-generated `alias_ai` that
    anno_cartload_folder writes; supplied via an `alias` field on a ficture entry."""
    dst = f"{oid}-alias.tsv"
    return [
        f"cp {alias_path} {os.path.join(cart_dir, dst)}",
        f"python3 -c \"from cartloader.utils.cartload_helper import record_catalog_alias; "
        f"record_catalog_alias('{catalog_path}', '{oid}', '{dst}')\"",
    ]


def cmd_record_catalog_name(catalog_path, fid, name):
    """Set a factor's human-readable `name` in a catalog. The factor **id** is untouched —
    every file name, asset key and cross-reference keeps using it; this is only the label
    shown for the layer. Without it a cell analysis displays its bare id (factor_id_to_name
    passes a non-FICTURE id through), and the multi-catalog's cell factors carry no name at
    all. Unrelated to `alias`, which points at a companion factor-label TSV."""
    return (f"python3 -c \"from cartloader.utils.cartload_helper import record_catalog_name; "
            f"record_catalog_name('{catalog_path}', '{fid}', '{name}')\"")


def cell_display_names(active_cells):
    """(analysis_id, name) for the active cell analyses that declare a display name."""
    return [(c["id"], c["name"]) for c in active_cells if c.get("name")]


def cmd_cartload_multi(fic_dir, cart_root, multi_id, cfg):
    """Package all samples of a joint run via run_cartload2_multi (one call,
    parallel internally, writes per-sample dirs + multi-catalog.yaml)."""
    cart = cfg.get("cartload", {})
    res = cfg["resources"]
    parts = ["cartloader", "run_cartload2_multi",
             f"--fic-dir {fic_dir}", f"--out-dir {cart_root}", f"--id {multi_id}",
             f"-j {res['n_jobs']}", f"--threads {res['threads']}",
             "--colname-count count", "--gzip pigz"]
    if cart.get("use_pmpoint"):
        parts.append("--use-pmpoint")
    if cart.get("sge_scale"):
        parts.append(f"--sge-scale {cart['sge_scale']}")
    if cart.get("bin_count") is not None:
        parts.append(f"--bin-count {cart['bin_count']}")
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Image / import planning
# ---------------------------------------------------------------------------

def resolve_image_ops(cfg, s):
    """Return concrete image operations (source resolved, deduped by id)."""
    ops, seen = [], set()
    # Per-sample images take precedence over the profile's for the same id.
    entries = s.get("images", []) + cfg.get("images", [])
    for spec in entries:
        iid = spec.get("id")
        if iid in seen:
            continue
        if "source" in spec:
            src = spec["source"]
            if s.get("in_dir") and not os.path.isabs(src) and not os.path.exists(src):
                cand = os.path.join(s["in_dir"], src)
                src = cand if os.path.exists(cand) else src
            if not os.path.exists(src):
                continue
        elif "suffix" in spec:   # appended to the sample's in_prefix
            src = first_existing(s["in_prefix"] + spec["suffix"]) if s.get("in_prefix") else None
            if not src:
                continue
        else:  # match relative to in_dir
            src = first_existing(os.path.join(s["in_dir"], spec["match"])) if s.get("in_dir") else None
            if not src:
                continue
        seen.add(iid)
        ops.append({**spec, "source": src})
    return ops


def catalog_image_line(cfg, catalog, iid, cart_dir):
    """Append an image entry to the catalog. When --skip-image-errors is on, only
    add it if the PMTiles actually exists, so a skipped corrupt image is left out
    of the catalog instead of leaving a dangling reference."""
    line = f"echo '    {iid}: {iid}.pmtiles' >> {catalog}"
    if cfg.get("_skip_image_errors"):
        pmt = os.path.join(cart_dir, f"{iid}.pmtiles")
        return f"[ -f {pmt} ] && {line} || echo 'WARNING: image {iid} skipped; omitted from catalog' >&2"
    return line


def plan_images(cfg, s, cart_dir, multi, transcript=None):
    cmds = []
    catalog = os.path.join(cart_dir, "catalog.yaml")
    skip_img = "--skip-image-errors " if cfg.get("_skip_image_errors") else ""

    for op in resolve_image_ops(cfg, s):
        iid, kind, src = op["id"], op.get("kind", "single"), op["source"]
        if kind == "prebuilt":
            cmds.append(f"cp {src} {os.path.join(cart_dir, iid + '.pmtiles')}")
        elif kind == "rgb":
            cmds += _cmd_rgb_image(cfg, s, iid, src, cart_dir, op)
            continue  # _cmd_rgb_image appends its own catalog line
        else:  # single-channel colorized
            conv = "--ome2png " if op.get("convert", "ome2png") == "ome2png" else ""
            idef = cfg.get("image_defaults", {})
            # platform geometric transform (e.g. MERFISH --micron2pixel-csv), and
            # shrink/high-memory (per-image, else profile image_defaults for big mosaics).
            transform = f"{op['transform_flag']} {op['transform_path']} " if op.get("transform_path") else ""
            sfval = op.get("shrink_factor", idef.get("shrink_factor"))
            sf = f"--shrink-factor {sfval} " if sfval is not None else ""
            hm = "--high-memory " if _truthy(op.get("high_memory", idef.get("high_memory"))) else ""
            # A plain (non-OME) TIFF carries no pixel size, so state the scale directly:
            # a Stereo-seq *_regist.tif at 0.5 um/pixel georeferences as px_per_um = 2.
            upp = op.get("um_per_pixel", idef.get("um_per_pixel"))
            if upp:
                ppu = 1.0 / float(upp)
                transform += f"--px-per-um-x {ppu:g} --px-per-um-y {ppu:g} "
            extra_flags = list(op.get("extra_flags", []))
            # multi-page z-stack: pick the middle page (by file name unless overridden)
            ump = op.get("use_middle_page")
            ump = default_use_middle_page(src) if ump is None else _truthy(ump)
            if ump and "--use-middle-page" not in extra_flags:
                extra_flags.append("--use-middle-page")
            extra = " ".join(extra_flags)
            cmds.append(
                f"cartloader import_image {conv}{skip_img}--png2pmtiles --georeference "
                f"--in-img {src} --out-dir {cart_dir} --img-id {iid} "
                f"--upper-thres-quantile 0.95 --level 0 --colorize {op['color']} "
                f"--transparent-below 5 {transform}{sf}{hm}{extra}".strip())
        cmds.append(catalog_image_line(cfg, catalog, iid, cart_dir))

    # Visium HD H&E (rgb via a per-sample path + profile hne settings)
    if s.get("hne") and cfg.get("hne"):
        cmds += _cmd_rgb_image(cfg, s, cfg["hne"].get("out_id", "hne"), s["hne"], cart_dir, cfg["hne"])

    # Visium HD square-bin imports
    for sq in cfg.get("squares", []):
        sub = os.path.join(s["in_dir"], sq["subdir"]) if s.get("in_dir") else None
        if sub and os.path.exists(sub):
            cmds.append(f"cartloader import_visiumhd_square --in-dir {sub} "
                        f"--outprefix {os.path.join(cart_dir, sq['suffix'])} "
                        f"--bin-size {sq['bin_size']} --update-catalog")
    ci = cfg.get("cell_import")
    if ci and s.get("in_dir") and os.path.exists(os.path.join(s["in_dir"], ci["detect_dir"])):
        cmds.append(f"cartloader import_visiumhd_cell --in-dir {s['in_dir']} "
                    f"--outprefix {os.path.join(cart_dir, ci['suffix'])} --all --update-catalog")

    # Per-sample cluster imports for joint runs (e.g. Xenium Ranger clusters):
    # a cell analysis flagged with `multi_import` is imported per sample here
    # instead of being jointly decoded by run_ficture2_multi_cells. Sheet-provided
    # role paths are forwarded as --csv-* overrides so GEO-style layouts (with
    # non-standard filenames or scattered paths) also work on joint runs.
    IMPORT_ROLE_FLAG = {"xy": "--csv-cells", "boundaries": "--csv-boundaries",
                        "clusters": "--csv-clust"}
    if multi:
        for ca in cfg.get("cell_analyses", []):
            imp = ca.get("multi_import")
            if not imp:
                continue
            # A sample supplies the inputs via sheet role columns, or a Ranger-style
            # --in-dir whose standard files were auto-detected into the same roles (a role
            # is only resolved from in_dir when its file exists). The import needs every
            # role in `uses` -- e.g. import_xenium_cell cannot run without the clusters
            # CSV -- so it is skipped, with a note, when any is missing.
            has_indir = bool(s.get("in_dir"))
            uses = ca.get("uses", [])
            role_paths = {r: s["roles"][r] for r in uses if s["roles"].get(r)}
            missing = [r for r in uses if r not in role_paths]
            if missing:
                if role_paths or has_indir:
                    print(f"NOTE: sample {s['id']}: skipping the per-sample import '{ca['id']}' "
                          f"({imp}); missing input role(s): {', '.join(missing)}.", file=sys.stderr)
                continue
            parts = [f"cartloader {imp}",
                     f"--in-dir {s['in_dir'] if has_indir else '.'}",
                     f"--outprefix {os.path.join(cart_dir, ca['id'])}",
                     "--all --update-catalog"]
            for role, path in role_paths.items():
                flag = IMPORT_ROLE_FLAG.get(role)
                if flag:
                    parts.append(f"{flag} {path}")
            # Regenerate pseudobulk/DE from the run_together transcript TSV (carries
            # cell_id + gene + count, matching the importer's --pixel defaults) so a
            # separate cell-feature MEX directory is not required.
            if imp == "import_xenium_cell" and transcript:
                parts.append(f"--pixel {transcript}")
            cmds.append(" ".join(parts))
    return cmds


def _cmd_rgb_image(cfg, s, iid, src, cart_dir, settings):
    catalog = os.path.join(cart_dir, "catalog.yaml")
    prefix = os.path.join(cart_dir, iid)
    cmds, upp = [], ""
    idef = cfg.get("image_defaults", {})
    # Rescale controls for 16-bit imagery (e.g. some Stereo-seq H&E TIFs), forwarded to
    # image_png2pmtiles -> geotiff2pmtiles, which rejects 16-bit input without a range.
    # The range is either rescale_range ("min,max", for JSON/TSV) or rescale_min +
    # rescale_max (so it survives the comma split in a --image CLI value). Applies to
    # every branch below (all run geotiff2pmtiles).
    rs = settings.get("rescale", idef.get("rescale"))
    rrange = settings.get("rescale_range", idef.get("rescale_range"))
    if not rrange:
        rmin = settings.get("rescale_min", idef.get("rescale_min"))
        rmax = settings.get("rescale_max", idef.get("rescale_max"))
        if rmin is not None and rmax is not None:
            rrange = f"{rmin},{rmax}"
    rescale = ""
    if rs:
        rescale += f" --rescale {rs}"
    if rrange:
        rescale += f" --rescale-range {rrange}"
    # Rotate/flip (incl. swapxy, the x/y transpose). Orthogonal to how the bounds are
    # obtained, so it applies to every branch below.
    orient = _image_orient_flags(settings, idef)
    # An image that already carries the geotransform we want (e.g. a Seq-Scope H&E TIF
    # registered upstream, whose corner coordinates are already the transcript um extent):
    # its bounds are taken as-is, so no bounds have to be synthesized and --georef-plain/
    # --um-per-pixel do not apply. It still goes through the georeference step, because
    # such a TIF typically carries no CRS, and the tilers derive their zoom range from one
    # — given none, geotiff2pmtiles settles on max zoom 0 and writes an empty PMTiles
    # without failing. --georef-detect gtiff assigns the CRS and touches nothing else.
    if _truthy(settings.get("georeferenced", idef.get("georeferenced"))):
        cmds.append(f"cartloader image_png2pmtiles --in-img {src} --out-prefix {prefix} "
                    f"--geotif2mbtiles --mbtiles2pmtiles --georeference "
                    f"--georef-detect gtiff{orient}{rescale}")
        cmds.append(catalog_image_line(cfg, catalog, iid, cart_dir))
        return cmds
    # Bounds source. An OME-TIFF carries its pixel size in embedded metadata, so the
    # bounds are read from it (--georef-detect ome) — this is the default for an
    # .ome.tif(f) RGB image (e.g. an Illumina registered H&E) unless the profile/sheet
    # gives an explicit scale instead. A plain image has no such metadata and needs
    # --georef-plain with a stated um/pixel below.
    detect = settings.get("georef_detect", idef.get("georef_detect"))
    if (not detect and not settings.get("georef_plain")
            and not settings.get("um_per_pixel") and not settings.get("um_per_pixel_json")
            and src.lower().endswith((".ome.tif", ".ome.tiff"))):
        detect = "ome"
    if detect:
        cmds.append(f"cartloader image_png2pmtiles --in-img {src} --out-prefix {prefix} "
                    f"--geotif2mbtiles --mbtiles2pmtiles --georeference --georef-detect {detect}"
                    f"{orient}{rescale}")
        cmds.append(catalog_image_line(cfg, catalog, iid, cart_dir))
        return cmds
    jrel = settings.get("um_per_pixel_json")
    if jrel and s.get("in_dir"):
        jpath = os.path.join(s["in_dir"], jrel)
        key = settings.get("um_per_pixel_key", "microns_per_pixel")
        # Inline the substitution: make runs each recipe line in its own shell,
        # so a UPP=... on a separate line would not survive to this command.
        upp = f"--um-per-pixel \"$(jq -r '.{key}' {jpath})\""
    elif settings.get("um_per_pixel"):
        # A fixed scale, for a registered image whose um/pixel is a property of the
        # platform rather than of the run (e.g. Stereo-seq *_regist.tif at 0.5).
        upp = f"--um-per-pixel {settings['um_per_pixel']}"
    plain = "--georef-plain" if settings.get("georef_plain") else ""
    cmds.append(f"cartloader image_png2pmtiles --in-img {src} --out-prefix {prefix} "
                f"--geotif2mbtiles --mbtiles2pmtiles --georeference {plain} {upp}{orient}{rescale}".strip())
    cmds.append(catalog_image_line(cfg, catalog, iid, cart_dir))
    return cmds


def cmd_anno(cart_dir, args, multi=False):
    """AI-annotate a packaged directory. For a joint run (multi=True) point at the
    cartl/ root: shared factors are annotated once and reused into every sample."""
    return (f"cartloader anno_cartload_folder --cartl-dir {cart_dir} "
            + ("--multi-sample " if multi else "")
            + f"--tissue \"{args.tissue}\" --organism {args.organism} "
            f"--api-type {args.anno_api_type} --model {args.anno_model} --threads {args.anno_threads} --profile {args.aws_profile}")


def cmd_upload(cart_dir, args, batch):
    """Upload one sample's packaged directory to S3.

    Destination: <s3_prefix>/batch=<YYYY_MM>/<collection>/<dir>, where <dir> is the
    self-contained sample directory name (<sample_id> or <multi_id>-<sample_id>).
    For a joint --out-dir run <collection> is the out-dir basename; for a
    single-sample --out-dir run it is the out-dir's PARENT basename, so the
    sample dir does not repeat as <collection>/<collection>. See run_together().
    """
    dest_id = os.path.basename(os.path.normpath(cart_dir))
    dest = f"{args.s3_prefix.rstrip('/')}/batch={batch}/{args.collection}/{dest_id}"
    catalog = os.path.join(cart_dir, "catalog.yaml")
    return ("(grep -E \"\\.\" " + catalog + " | perl -lane 'print $F[$#F]' | sort | uniq; "
            "echo catalog.yaml;) | xargs -I {} -P " + str(args.s3_jobs) + " " + args.aws + " s3 cp " + cart_dir + "/{} "
            + dest + "/{} --profile " + args.aws_profile)


def cmd_upload_multi_catalog(cart_root, args, batch):
    """Upload multi-catalog.yaml and the shared factor files at the cartl/ root to
    <s3_prefix>/batch=<YYYY_MM>/<collection>/ — the parent of the per-sample dirs,
    so the catalog's relative pointers resolve on S3."""
    mc = os.path.join(cart_root, "multi-catalog.yaml")
    dest = f"{args.s3_prefix.rstrip('/')}/batch={batch}/{args.collection}"
    # Shared files are the flat (no-slash) basenames in multi-catalog.yaml; the
    # per-sample 'samples:' entries contain '/' and are already uploaded elsewhere.
    return ("(grep -E \"\\.\" " + mc + " | perl -lane 'print $F[$#F]' | grep -v / | sort | uniq; "
            "echo multi-catalog.yaml;) | xargs -I {} -P " + str(args.s3_jobs) + " " + args.aws + " s3 cp " + cart_root + "/{} "
            + dest + "/{} --profile " + args.aws_profile)


# ---------------------------------------------------------------------------
# Makefile assembly
# ---------------------------------------------------------------------------

def add_targets(mm, samples, cfg, args):
    stages = set(args.only.split(",")) if args.only else None
    skip = set(args.skip.split(",")) if args.skip else set()

    def on(stage):
        return stage not in skip and (stages is None or stage in stages)

    _, default_model_id = resolved_models(cfg["ficture"])
    # Every sample id in the run (not just this group's): a user-supplied list file may
    # legitimately cover the whole run while each --out-root group sees one sample.
    cfg["_all_sids"] = [s["id"] for s in samples]

    groups = {}
    for s in samples:
        groups.setdefault(s["out_dir"], []).append(s)

    for out_dir, grp in groups.items():
        sge_root = os.path.join(out_dir, "tsv")
        fic_dir = os.path.join(out_dir, "fic")
        cart_root = os.path.join(out_dir, "cartl")
        mkdir = os.path.join(out_dir, "mk")
        os.makedirs(mkdir, exist_ok=True)
        os.makedirs(sge_root, exist_ok=True)

        # --- ingest (skipped for samples that already provide a transcript) ---
        # An ingest 'method' of "reformat_cosmx" replaces sge_convert with a custom
        # step that emits multiple role files (transcript + xy + boundaries); their
        # produced paths are injected as this sample's roles so the cells stage picks
        # them up (unless the sample already supplies that role explicitly).
        ing = cfg.get("ingest", {})
        method = ing.get("method", "sge_convert")
        produces = ing.get("produces", {})
        sge_flags, transcript = [], {}
        preingested = []
        for s in grp:
            sge_dir = os.path.join(sge_root, s["id"])
            # .parquet inputs are converted to .csv.gz first (see plan_parquet_conversions);
            # the conversion is its own ingest-stage target that the sample's ingest and,
            # through sge_flags, every later stage depend on.
            conv_flags = []
            conv_cmds = plan_parquet_conversions(s, sge_dir, cfg)
            if conv_cmds:
                cflag = os.path.join(mkdir, f"parquet2csv.{s['id']}.done")
                if on("ingest"):
                    mm.add_target(cflag, [], [f"mkdir -p {os.path.join(sge_dir, 'parquet2csv')}"]
                                  + conv_cmds + [f"touch {cflag}"])
                conv_flags.append(cflag)
            if s["roles"].get("transcript"):
                transcript[s["id"]] = s["roles"]["transcript"]
                preingested.append(s["id"])
                sge_flags.extend(conv_flags)
                continue
            if method == "reformat_cosmx":
                prefix = os.path.join(sge_dir, s["id"])
                transcript[s["id"]] = prefix + produces["transcript"]
                for role, suffix in produces.items():
                    if role != "transcript" and not s["roles"].get(role):
                        s["roles"][role] = prefix + suffix
                ingest_cmds = [cmd_reformat_cosmx(cfg, sge_dir, s["id"], s["in_dir"])]
            elif method == "stereoseq":
                # SAW-driven ingest. The bin1 GEM becomes the ordinary pixel transcript;
                # the cell-bin GEM, when present, becomes a standalone `cell_tsv` role
                # that the cells stage passes as --tsv-list (the cell assignment cannot
                # be carried on the transcript itself). An explicit role still wins.
                transcript[s["id"]] = os.path.join(sge_dir, "transcripts.unsorted.tsv.gz")
                ingest_cmds, cell_tsv = cmds_stereoseq_ingest(cfg, sge_dir, s)
                if cell_tsv and not s["roles"].get("cell_tsv"):
                    s["roles"]["cell_tsv"] = cell_tsv
            else:
                base_tx = os.path.join(sge_dir, "transcripts.unsorted.tsv.gz")
                ingest_cmds = [cmd_sge_convert(cfg, sge_dir, s)]
                # The transcript already carries a cell_id column (carried to column 5 by
                # sge_convert): use it directly and skip boundary assignment.
                if s.get("has_tx_cell_id"):
                    transcript[s["id"]] = base_tx
                # Otherwise optionally assign each transcript to its overlapping cell
                # boundary so the transcript gains a cell_id column (MERFISH ships no
                # per-transcript cell assignment). Only when the sample supplies a
                # boundaries role; without it the run stays pixel-level (no cell decode).
                elif ing.get("assign_cell_id") and s["roles"].get("boundaries"):
                    tx_cellid = os.path.join(sge_dir, "transcripts.cell_id.tsv.gz")
                    ingest_cmds.append(cmd_tsv_add_cell_id(cfg, base_tx, tx_cellid, s["roles"]["boundaries"]))
                    transcript[s["id"]] = tx_cellid
                else:
                    transcript[s["id"]] = base_tx
                # A cell-by-gene matrix CSV (e.g. MERSCOPE cell_by_gene.csv) is converted
                # to a MEX directory that becomes this sample's `mex` role. It drives the
                # cell decode's clustering (mex2sptsv), so cell-based analysis works with
                # no cell boundaries; when boundaries are also present they still supply
                # the transcript cell_id column (above) and polygon rendering, while the
                # MEX overrides the clustering source. An explicit mex role wins.
                if s["roles"].get("cellxgene") and not s["roles"].get("mex"):
                    mex_dir = os.path.join(sge_dir, "cellxgene_mex")
                    ingest_cmds.append(cmd_convert_cellxgene(cfg, s["roles"]["cellxgene"], mex_dir))
                    s["roles"]["mex"] = mex_dir
            flag = os.path.join(mkdir, f"sge.{s['id']}.done")
            sge_flags.append(flag)
            if on("ingest"):
                mm.add_target(flag, conv_flags, [f"mkdir -p {sge_dir}"] + ingest_cmds + [f"touch {flag}"])

        # The ingest filters act while the transcript TSV is written, so a sample that
        # supplies a ready-made transcript never sees them: its TSV is taken as given.
        # Only filters this run explicitly asked for are worth reporting; the built-in
        # default not applying to a pre-ingested file is unremarkable.
        ingest_filters = cfg.get("_ingest_filters_explicit", [])
        if ingest_filters and preingested:
            print(f"WARNING: {', '.join(ingest_filters)} had no effect on sample(s) "
                  f"{', '.join(preingested)}: they provide an already-ingested transcript, which is "
                  f"used as-is. Apply the same beforehand to that file, or supply the raw input "
                  f"instead so it goes through ingest.", file=sys.stderr)

        in_list = os.path.join(sge_root, "in_list.tsv")
        with open(in_list, "w") as f:
            for s in grp:
                f.write(f"{s['id']}\t{transcript[s['id']]}\n")

        # --- ficture: one target, each analysis run sequentially (params merge) ---
        fic_flag = os.path.join(mkdir, "ficture.done")
        if on("ficture"):
            cmds = [cmd_ficture_analysis(a, in_list, fic_dir, cfg) for a in cfg["ficture"]]
            mm.add_target(fic_flag, list(sge_flags), cmds + [f"touch {fic_flag}"])

        # A joint run (>1 sample sharing out_dir) packages every sample with a
        # single run_cartload2_multi call, laid out as <multi_id>-<sample_id>/ plus
        # a multi-catalog.yaml. A single sample uses run_cartload2 directly.
        multi = len(grp) > 1
        multi_id = os.path.basename(os.path.normpath(out_dir))

        # Companion alias files supplied on ficture entries (e.g. a curated label set
        # shipped with a projection model): each is deployed beside its factor as
        # <oid>-alias.tsv and recorded under the catalog's `alias` key. oid is the
        # hyphenated model id, matching the factor ids in the deployed catalog(s).
        alias_specs = [(str(a["id"]).replace("_", "-"), a["alias"])
                       for a in cfg["ficture"] if a.get("alias")]

        # --- cell analyses (platform default; run those whose roles are present) ---
        cells_flag = os.path.join(mkdir, "cells.done")
        # Cell analyses decode against a trained model, so they cannot run without one.
        active_cells = [] if cfg.get("_no_ficture") else \
            plan_cell_analyses(grp, sge_root, cfg, fic_dir, default_model_id, multi)
        if active_cells and on("cells"):
            cmds = [c["cmd"] for c in active_cells]
            mm.add_target(cells_flag, [fic_flag], cmds + [f"touch {cells_flag}"])

        cart_prereq = cells_flag if (active_cells and on("cells")) else fic_flag
        multi_cart_flag = os.path.join(mkdir, "cartload.done")
        if multi and on("cartload"):
            cmds = [f"mkdir -p {cart_root}",
                    cmd_cartload_multi(fic_dir, cart_root, multi_id, cfg)]
            # Record each alias in the shared multi-catalog and in every per-sample
            # catalog (run_cartload2_multi has written all of them by this point).
            for oid, alias_path in alias_specs:
                cmds += cmd_record_alias(cart_root, os.path.join(cart_root, "multi-catalog.yaml"), oid, alias_path)
                for s in grp:
                    sample_cart = os.path.join(cart_root, f"{multi_id}-{s['id']}")
                    cmds += cmd_record_alias(sample_cart, os.path.join(sample_cart, "catalog.yaml"), oid, alias_path)
            # Display names for cell analyses that declare one, in the shared multi-catalog
            # (whose factor keys are hyphenated) and in each contributing sample's catalog
            # (whose factor ids are the analysis id verbatim). run_cartload2_multi has
            # written all of them by this point.
            for cid, cname in cell_display_names(active_cells):
                cmds.append(cmd_record_catalog_name(
                    os.path.join(cart_root, "multi-catalog.yaml"), cid.replace("_", "-"), cname))
                for c in active_cells:
                    if c["id"] != cid:
                        continue
                    for sid in c["sids"]:
                        sample_cart = os.path.join(cart_root, f"{multi_id}-{sid}")
                        cmds.append(cmd_record_catalog_name(
                            os.path.join(sample_cart, "catalog.yaml"), cid, cname))
            mm.add_target(multi_cart_flag, [cart_prereq], cmds + [f"touch {multi_cart_flag}"])

        # --- cartload + images (per sample); collect each sample's post-images flag ---
        sample_ctx = []   # (sample, cart_dir, base_prereq)
        for s in grp:
            fic_sample_dir = os.path.join(fic_dir, "samples", s["id"])
            if multi:
                cart_dir = os.path.join(cart_root, f"{multi_id}-{s['id']}")
                cart_flag = multi_cart_flag
            else:
                # Compose the sample directory name and catalog id/title as
                # <out_dir>-<sample>, matching the <multi_id>-<sample_id> layout and ids
                # used in multi-sample runs, so a single-sample output carries the
                # collection context rather than just the bare sample name (e.g.
                # my-collection-rep1). When the out_dir basename already equals the
                # sample id — an --out-root run, or an --out-dir whose basename happens
                # to be the id — keep the bare id to avoid a redundant "rep1-rep1".
                catalog_id = s["id"] if multi_id == s["id"] else f"{multi_id}-{s['id']}"
                cart_dir = os.path.join(cart_root, catalog_id)
                cart_flag = os.path.join(mkdir, f"cartload.{s['id']}.done")
                cell_params = [os.path.join(fic_sample_dir, f"ficture.{c['id']}.params.json")
                               for c in active_cells if s["id"] in c["sids"]]
                if on("cartload"):
                    cmds = [f"mkdir -p {cart_dir}",
                            cmd_cartload(fic_sample_dir, cart_dir, catalog_id, cfg, cell_params)]
                    for oid, alias_path in alias_specs:
                        cmds += cmd_record_alias(cart_dir, os.path.join(cart_dir, "catalog.yaml"), oid, alias_path)
                    for cid, cname in cell_display_names(active_cells):
                        if s["id"] in next(c["sids"] for c in active_cells if c["id"] == cid):
                            cmds.append(cmd_record_catalog_name(
                                os.path.join(cart_dir, "catalog.yaml"), cid, cname))
                    mm.add_target(cart_flag, [cart_prereq], cmds + [f"touch {cart_flag}"])

            img_prereq = cart_flag if on("cartload") else cart_prereq
            img_flag = os.path.join(mkdir, f"images.{s['id']}.done")
            img_cmds = plan_images(cfg, s, cart_dir, multi, transcript.get(s["id"]))
            if img_cmds and on("images"):
                mm.add_target(img_flag, [img_prereq], img_cmds + [f"touch {img_flag}"])
            base_prereq = img_flag if (img_cmds and on("images")) else img_prereq
            sample_ctx.append((s, cart_dir, base_prereq))

        # --- annotate (opt-in) ---
        # A joint run annotates the shared factors once at the cartl/ root
        # (--multi-sample) and reuses them into every sample; a single run
        # annotates its one directory. `anno_prereq[sid]` is what upload waits on.
        anno_prereq = {}
        if args.anno and on("anno"):
            if multi:
                anno_flag = os.path.join(mkdir, "anno.done")
                mm.add_target(anno_flag, [bp for (_, _, bp) in sample_ctx],
                              [cmd_anno(cart_root, args, multi=True), f"touch {anno_flag}"])
                anno_prereq = {s["id"]: anno_flag for (s, _, _) in sample_ctx}
            else:
                for (s, cart_dir, bp) in sample_ctx:
                    anno_flag = os.path.join(mkdir, f"anno.{s['id']}.done")
                    mm.add_target(anno_flag, [bp], [cmd_anno(cart_dir, args), f"touch {anno_flag}"])
                    anno_prereq[s["id"]] = anno_flag

        # --- S3 upload (opt-in); upload waits on annotation when both run ---
        if args.s3_upload and on("upload"):
            upload_flags = []
            for (s, cart_dir, bp) in sample_ctx:
                up_flag = os.path.join(mkdir, f"upload.{s['id']}.done")
                mm.add_target(up_flag, [anno_prereq.get(s["id"], bp)],
                              [cmd_upload(cart_dir, args, args.batch), f"touch {up_flag}"])
                upload_flags.append(up_flag)
            # joint run: also upload multi-catalog.yaml + shared files to the parent dir
            if multi:
                mc_flag = os.path.join(mkdir, "upload.multi-catalog.done")
                mm.add_target(mc_flag, upload_flags,
                              [cmd_upload_multi_catalog(cart_root, args, args.batch), f"touch {mc_flag}"])


def _role_list_line(sid, role, val):
    """One --list-* line: a 4-column bcd/ftr/mtx triple for an explicit mex dict,
    else the usual `id<TAB>path`."""
    if role == "mex" and isinstance(val, dict):
        return f"{sid}\t{val['bcd']}\t{val['ftr']}\t{val['mtx']}\n"
    return f"{sid}\t{val}\n"


def _check_list_ids(path, role, ca_id, grp_sids, all_sids):
    """Sanity-check a user-supplied --list-* file against the run's sample ids.

    The first column is the join key that ties each listed file to a sample's ingested
    transcript, FICTURE subdir and packaged output, so an id matching nothing in this
    analysis's group is a hard error (a typo, or a list written for a different sample
    set). Everything else is only reported: incomplete coverage is legitimate (a run may
    carry labels for some samples only), and an id outside the group but inside the run is
    normal for an --out-root run, where one list can serve every per-sample group."""
    ids = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                ids.append(line.split("\t")[0])
    if not ids:
        sys.exit(f"ERROR: empty list file: {path} (cell analysis '{ca_id}', lists.{role})")
    grp, seen = set(grp_sids), set(ids)
    if not (seen & grp):
        sys.exit(f"ERROR: none of the sample ids in {path} (cell analysis '{ca_id}', "
                 f"lists.{role}) match the samples it applies to. Listed: "
                 f"{', '.join(sorted(seen))}. Expected one or more of: {', '.join(sorted(grp))}.")
    unknown = sorted(seen - set(all_sids))
    if unknown:
        print(f"WARNING: {path} (cell analysis '{ca_id}', lists.{role}) lists sample id(s) "
              f"not in this run, which will be ignored: {', '.join(unknown)}", file=sys.stderr)
    missing = sorted(grp - seen)
    if missing:
        print(f"WARNING: {path} (cell analysis '{ca_id}', lists.{role}) has no entry for "
              f"sample(s): {', '.join(missing)}", file=sys.stderr)


def _resolve_cell_inputs(ca, grp, cfg):
    """Decide whether a cell analysis runs for this group and which file roles feed
    it. Returns (contributing_samples, required_list_roles) or (None, None) if it
    should not run.

    A `generic_cell` analysis is selected explicitly (generic platform):
      (a) --colname-cell given      -> cell_id-based (cell_id read positionally from
                                        the tiled transcript; no file role required),
      (b) mex present in every sample -> mex-based clustering,
      (c) neither                   -> skipped.
    Otherwise the analysis is role-driven: a sample contributes when it provides all of
    `uses` and (if `any_uses` is set) at least one of `any_uses` — so an analysis can
    accept alternative cell-count sources (e.g. MERSCOPE: boundaries OR a cellxgene MEX).
    `require_all` demands every sample in the group contribute. The pseudo-role
    ``tx_cell_id`` matches a sample whose transcript already carries a cell_id column
    (has_tx_cell_id); it feeds no --list file (cell_id is read positionally from the
    tiled transcript, like the generic platform's --colname-cell path).

    A role named in the analysis's `lists` override counts as present for every sample:
    the list was supplied directly, so no sample has to carry that role on disk for the
    analysis to run (see explicit_lists).
    """
    explicit = explicit_lists(ca)

    def _has(s, r):
        if r in explicit:
            return True
        return bool(s.get("has_tx_cell_id")) if r == "tx_cell_id" else bool(s["roles"].get(r))
    if ca.get("generic_cell"):
        mex_ok = all(s["roles"].get("mex") for s in grp)
        cellid_ok = bool(cfg.get("colname_cell"))
        if mex_ok and cellid_ok:
            sys.exit("ERROR: ambiguous cell analysis for the generic platform: both "
                     "--colname-cell and mex inputs are present. Provide only one.")
        if mex_ok:
            required = ["mex"]
        elif cellid_ok:
            required = []          # cell_id comes from the tiled transcript, no list file
        else:
            return None, None
        contributing = list(grp)   # applies to every sample in the group
        candidate_roles = required + list(ca.get("optional_uses", []))
    else:
        uses = list(ca.get("uses", []))
        any_uses = list(ca.get("any_uses", []))
        def _contributes(s):
            if not all(_has(s, r) for r in uses):
                return False
            if any_uses and not any(_has(s, r) for r in any_uses):
                return False
            return True
        contributing = [s for s in grp if _contributes(s)]
        if not contributing:
            return None, None
        # Analyses that decode all samples jointly need every sample to contribute.
        if ca.get("require_all") and len(contributing) < len(grp):
            return None, None
        candidate_roles = uses + any_uses + list(ca.get("optional_uses", []))

    # A role feeds a --list file when at least one contributing sample supplies it; the
    # list then contains only those samples. run_ficture2_multi_cells resolves each
    # sample's cell-count source per-sample (MEX vs tiled transcript) and tolerates
    # partial xy/boundaries/cluster coverage, so a mixed joint run (some samples with
    # boundaries, some with a cellxgene MEX) is packaged from one call.
    # An explicitly named list is always passed, including for a role the analysis does
    # not declare in uses/any_uses/optional_uses (e.g. adding cluster labels to an
    # analysis that would otherwise cluster on demand).
    list_roles = []
    for r in list(explicit) + candidate_roles:
        if r not in ROLE_LIST_FLAG or r in list_roles:
            continue
        if r in explicit or any(s["roles"].get(r) for s in contributing):
            list_roles.append(r)
    return contributing, list_roles


def plan_cell_analyses(grp, sge_root, cfg, fic_dir, default_model_id, multi):
    """Write per-cell-analysis role lists; return the analyses that have inputs, each
    tagged with the sample ids that contribute (`sids`).

    An analysis with a ``multi_import`` command is a per-sample import for joint
    runs (e.g. Xenium Ranger clusters, which are sample-specific and cannot be
    jointly decoded) — it is skipped here for multi and handled in the image stage. A
    `lists` override does not change that routing: such an analysis is per-sample by
    nature, so a joint list belongs on one of the jointly decoded analyses instead.
    """
    active = []
    for ca in cfg.get("cell_analyses", []):
        if multi and ca.get("multi_import"):
            continue
        contributing, list_roles = _resolve_cell_inputs(ca, grp, cfg)
        if contributing is None:
            continue
        bnd = cfg.get("roles", {}).get("boundaries", {})
        explicit = explicit_lists(ca)
        list_files = {}
        for role in list_roles:
            # A list named in `lists` is passed through verbatim — the sample ids in it
            # are the join key for every later stage, so they are checked against the
            # run's own ids rather than trusted blindly.
            if role in explicit:
                list_files[role] = explicit[role]
                _check_list_ids(explicit[role], role, ca["id"],
                                [s["id"] for s in contributing], cfg.get("_all_sids", []))
                continue
            # Only the contributing samples that actually supply this role (a mixed run
            # lists, e.g., boundaries for some samples and mex for others).
            samples_with = [s for s in contributing if s["roles"].get(role)]
            path = os.path.join(sge_root, f"in_{role}.{ca['id']}.tsv")
            with open(path, "w") as f:
                for s in samples_with:
                    # For geojson boundaries, append the sample's scale JSON as a third
                    # column so the decode can rescale polygon coords into microns when
                    # deriving centroids (units_json is relative to the sample's in_dir).
                    if role == "boundaries" and bnd.get("units_json") and s.get("in_dir"):
                        scale_json = _abs_in_dir(bnd["units_json"], s["in_dir"])
                        f.write(f"{s['id']}\t{s['roles'][role]}\t{scale_json}\n")
                    else:
                        f.write(_role_list_line(s["id"], role, s["roles"][role]))
            list_files[role] = path
        model_id = ca.get("model_id", default_model_id)
        model_path = os.path.join(fic_dir, f"{model_id}.model.tsv")
        active.append({"id": ca["id"], "sids": [s["id"] for s in contributing],
                       "name": ca.get("name"),
                       "cmd": cmd_cells(ca, list_files, fic_dir, model_path, cfg)})
    return active


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_arguments(_args):
    p = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="End-to-end multi-platform pipeline (ingest -> FICTURE -> cells -> "
                    "cartload -> images -> publish) assembled as a single resumable Makefile.")
    r = p.add_argument_group("Run Options")
    r.add_argument("--dry-run", action="store_true", help="Write the Makefile and print commands without executing")
    r.add_argument("--restart", action="store_true", help="Ignore existing outputs and rebuild (make -B)")
    r.add_argument("-j", "--n-jobs", type=int, default=1, help="Parallel jobs (default: 1)")
    r.add_argument("--threads", type=int, default=4, help="Threads per job (default: 4)")
    r.add_argument("--makefn", type=str, default="run_together.mk", help="Master Makefile name")
    r.add_argument("--only", type=str, help="Comma-separated stages to run exclusively (ingest,ficture,cells,cartload,images,anno,upload)")
    r.add_argument("--skip", type=str, help="Comma-separated stages to skip")
    r.add_argument("--skip-image-errors", action="store_true",
                   help="Tolerate an unreadable/corrupt OME-TIFF in the images stage: the image is "
                        "warned-and-skipped (and omitted from the catalog) instead of failing the run. "
                        "Re-run with this flag to regenerate the Makefile and resume past a corrupt image.")
    r.add_argument("--allow-unknown-config-keys", action="store_true",
                   help="Do not fail on unrecognized keys in --config (they are ignored). By default an "
                        "unknown key aborts the run, since a misplaced key (e.g. a top-level 'csv_colnames' "
                        "that belongs under 'ingest') is silently dropped otherwise. Checks top-level keys, "
                        "sample entries, the closed-schema blocks (ingest, cartload, resources, "
                        "image_defaults), and sample roles the platform cannot consume.")

    io = p.add_argument_group("Input/Output")
    io.add_argument("--platform", type=str, help="Platform preset, e.g. 10x_xenium, 10x_visium_hd, cosmx_smi, merfish, generic")
    io.add_argument("--in-dir", type=str, help="Input directory for a single sample (roles auto-detected inside it)")
    io.add_argument("--in-prefix", type=str,
                    help="Input path prefix for a single sample, for platforms whose files are named by "
                         "suffix rather than laid out in a directory (Stereo-seq). E.g. --in-prefix "
                         "/data/C04687E314 picks up C04687E314.tissue.gef, C04687E314.cellbin.gef and "
                         "C04687E314_{HE,ssDNA,DAPI}_regist.tif. Missing files are skipped; --image and "
                         "the role flags override any of them. The sample id defaults to the basename.")
    io.add_argument("--in-transcript", type=str, help="Single sample: a pre-converted (already ingested) transcript TSV in the "
                                                      "transcripts.unsorted.tsv.gz layout (X, Y, gene, count[, cell_id]). It skips ingest entirely "
                                                      "(the `transcript` role) on every platform — use it to re-run analysis on a filtered/edited "
                                                      "transcript without re-ingesting. For a RAW platform CSV that still needs sge_convert, use "
                                                      "--raw-transcript instead.")
    io.add_argument("--raw-transcript", type=str, help="Single sample: a raw transcript CSV to ingest through sge_convert. Use instead of --in-dir "
                                                      "when the raw file is arbitrarily named / not in a standard directory. Pair with the "
                                                      "--colname-transcript-* flags if its columns are non-standard. Not applicable to "
                                                      "--platform stereoseq (which ingests binary GEFs via SAW, not a raw CSV).")
    io.add_argument("--in-cell-xy", type=str, help="Single sample: cell metadata (centroids) file for the xy role (optional)")
    io.add_argument("--in-cell-boundary", type=str, help="Single sample: cell boundary polygon file for the boundaries role (optional)")
    io.add_argument("--in-cellxgene", type=str, help="Single sample: cell-by-gene matrix CSV (e.g. MERSCOPE cell_by_gene.csv). "
                                                     "Converted to a MEX directory that drives cell clustering, so cell-based "
                                                     "analysis works without cell boundaries (optional).")
    io.add_argument("--samples", type=str, help="TSV sample sheet (wide table of input roles) for a joint multi-sample run")
    io.add_argument("--images", type=str, help="TSV of per-image specs, one row per image (columns: sample, type, id, "
                                              "source/src, color, and the platform transform column e.g. merfish_csv). "
                                              "Images attach to samples by the `sample` column ('*'/blank = all samples).")
    io.add_argument("--image", action="append", metavar="type=..,source=..,..",
                    help="Single-sample image as comma-separated key=value pairs (keys: type, source/src, id, color, "
                         "transform (or the platform column e.g. merfish_csv), shrink_factor, high_memory, and the "
                         "orientation keys swapxy / rotate / flip_vertical / flip_horizontal — swapxy=true transposes "
                         "the image, for histology stored with its axes swapped relative to the transcript "
                         "coordinates). Same fields as a --images TSV row; repeat --image per image.")
    io.add_argument("--out-dir", type=str, help="Output directory (single sample / shared joint model)")
    io.add_argument("--out-root", type=str, help="Output root; each sample gets its own dir and an independent model")
    io.add_argument("--id", type=str, help=f"Sample id for a single-sample run (default: {DEFAULT_SAMPLE_ID}). The packaged directory and catalog id become <out-dir basename>-<id>, so the descriptive name lives in --out-dir. With --out-root instead, each sample is named after its input directory/prefix.")
    io.add_argument("--config", type=str, help="JSON config that augments the profile/CLI (full spec for complex runs)")
    io.add_argument("--platform-json", type=str, help="External JSON profile that overrides the built-in platform profile")
    io.add_argument("--saw", type=str, help="Path to the SAW binary (required for --platform stereoseq: the "
                                            ".gef inputs are binary and only SAW can expand them into text GEMs)")

    l = p.add_argument_group(
        "Ready-made cell input lists (override what run_together derives from the samples)")
    l.add_argument("--list-cluster", type=str, default=None,
                   help="TSV of externally assigned cell clusters, '[SAMPLE_ID] [CLUSTER_FILE]' per "
                        "line, passed to run_ficture2_multi_cells as --list-cluster. Without it the "
                        "cells stage computes Leiden clusters on demand. Each cluster file is 2 "
                        "columns (cell_id, cluster id; 1-based unless the analysis passes "
                        "--zero-based-clust-id via extra_flags).")
    l.add_argument("--list-xy", type=str, default=None,
                   help="TSV of cell metadata (centroid) files, '[SAMPLE_ID] [XY_FILE]' per line, "
                        "replacing the per-sample xy role files")
    l.add_argument("--list-boundaries", type=str, default=None,
                   help="TSV of cell boundary files, '[SAMPLE_ID] [BOUNDARY_FILE]' per line, "
                        "replacing the per-sample boundaries role files")
    l.add_argument("--list-mex", type=str, default=None,
                   help="TSV of cell-feature MEX dirs, '[SAMPLE_ID] [MEX_DIR]' (or a 4-column "
                        "bcd/ftr/mtx triple) per line, replacing the per-sample mex role files")
    l.add_argument("--list-cell-tsv", type=str, default=None,
                   help="TSV of standalone cell-level TSVs, '[SAMPLE_ID] [CELL_TSV]' per line, "
                        "replacing the per-sample cell_tsv role files")

    c = p.add_argument_group("Single-sample input column overrides (else profile / platform defaults)")
    c.add_argument("--colname-transcript-x", type=str, default=None, help="X column name in --raw-transcript")
    c.add_argument("--colname-transcript-y", type=str, default=None, help="Y column name in --raw-transcript")
    c.add_argument("--colname-transcript-feature", type=str, default=None, help="Gene/feature column name in --raw-transcript")
    c.add_argument("--colname-transcript-count", type=str, default=None, help="Count column name in --raw-transcript (default: none, count of 1 per row)")
    c.add_argument("--colname-transcript-cell", type=str, default=None,
                   help="Cell-id column name already present in the raw transcript CSV (rare). Carries the column "
                        "through ingest to transcript column 5 and runs cell analysis from it, skipping "
                        "spatula tsv-add-cell-id. A cellxgene MEX, if present, still drives clustering.")
    c.add_argument("--colname-xy-cell", type=str, default=None, help="Cell-id column name in --in-cell-xy (use '' for an unnamed pandas-index first column)")
    c.add_argument("--colname-xy-x", type=str, default=None, help="X (centroid) column name in --in-cell-xy")
    c.add_argument("--colname-xy-y", type=str, default=None, help="Y (centroid) column name in --in-cell-xy")
    c.add_argument("--colname-boundary-cell", type=str, default=None, help="Cell-id column name in --in-cell-boundary")
    c.add_argument("--colname-boundary-x", type=str, default=None, help="Vertex X column name in --in-cell-boundary")
    c.add_argument("--colname-boundary-y", type=str, default=None, help="Vertex Y column name in --in-cell-boundary")

    f = p.add_argument_group("FICTURE mode (choose one; default de-novo from the profile)")
    f.add_argument("--width", type=str, help="De-novo: hexagon width(s) in um (comma-separated)")
    f.add_argument("--n-factor", type=str, help="De-novo: factor count(s) (comma-separated)")
    f.add_argument("--project-models", type=str, help="Projection-only: existing FICTURE dir(s) (comma-separated); "
                                                      "reads each ficture.params.json and reuses its models. No LDA training.")
    f.add_argument("--no-ficture", action="store_true",
                   help="No factor analysis at all: run only FICTURE's tiling step and package the tiled "
                        "transcripts (points + raster + images, no factor layers). Cell analyses are skipped too.")
    f.add_argument("--segment-10x", action="store_true",
                   help="Also export the hexagon files as 10x MEX directories, one per sample and width "
                        "(<fic_dir>/samples/<id>/<id>.hex_<width>.mex/ with barcodes/features/matrix; "
                        "barcodes are hexagon centers 'x:y'). Uses the very same hexagons as the factor "
                        "analysis (after --min-ct-per-unit-hexagon). Works in every mode, including "
                        "--no-ficture. Config equivalent: \"segment_10x\": true.")
    f.add_argument("--segment-width-10x", type=str, default=None,
                   help="Hexagon width(s) in um to export with --segment-10x, comma-separated (default: "
                        "each analysis' own --width). Extra widths get their hexagon files built too. "
                        "Config equivalent: \"segment_10x\": {\"widths\": \"12,24\"}.")

    d = p.add_argument_group("Common decode overrides (else profile / built-in defaults)")
    d.add_argument("--exclude-feature-regex", type=str, default=None,
                   help=f"Regex of features to exclude from the factor analyses (default: "
                        f"'{DEFAULT_EXCLUDE_REGEX}' — predicted, mitochondrial and ribosomal genes, "
                        f"which stay in the data but out of the models). Ingest and packaging are "
                        f"unaffected. Pass '' to disable.")
    d.add_argument("--include-feature-list", type=str, default=None,
                   help="File listing the feature names (one per line) the factor analyses are restricted "
                        "to — pixel FICTURE (LDA training and pixel decoding) and the cell-based analysis. "
                        "Ingest, the tiled transcripts, the feature list and the packaged tiles keep every "
                        "gene. Combines with --exclude-feature-regex.")
    d.add_argument("--exclude-feature-list", type=str, default=None,
                   help="File listing the feature names (one per line) to drop from the factor analyses. "
                        "Same scope as --include-feature-list.")
    g = p.add_argument_group("Ingest feature filters (drop features from the data itself)")
    g.add_argument("--ingest-include-feature-regex", type=str, default=None,
                   help="Regex of features to keep when the transcript TSV is written. Unlike the "
                        "options above, these act on the very first step, so a dropped feature is "
                        "absent from the transcript TSV, the feature list, the packaged tiles and "
                        "every analysis. Independent of the factor-analysis filters. Ignored (with a "
                        "warning) for samples that supply an already-ingested transcript.")
    g.add_argument("--ingest-exclude-feature-regex", type=str, default=None,
                   help=f"Regex of features to drop when the transcript TSV is written (default: "
                        f"'{DEFAULT_INGEST_EXCLUDE_REGEX}' — negative controls, blanks and unassigned/"
                        f"deprecated codewords, which are technical artifacts rather than genes). "
                        f"Replaces sge_convert's per-platform default. Pass '' to disable.")
    g.add_argument("--ingest-include-feature-list", type=str, default=None,
                   help="File listing the feature names (one per line) to keep when the transcript TSV "
                        "is written. Same scope as --ingest-include-feature-regex.")
    g.add_argument("--ingest-exclude-feature-list", type=str, default=None,
                   help="File listing the feature names (one per line) to drop when the transcript TSV "
                        "is written. Same scope as --ingest-exclude-feature-regex.")

    j = p.add_argument_group("Ingest coordinate options (change the transcript coordinates themselves)")
    j.add_argument("--jitter-xy", type=float, default=None,
                   help="Add uniform random jitter of +/- this many microns to each transcript's X and Y "
                        "(independently) as the transcript TSV is written; 0 disables it (the default). "
                        "Use it on coarse-resolution platforms whose molecules all sit on the same lattice "
                        "points -- e.g. 10x Visium HD 2um bins, where every transcript in a bin otherwise "
                        "shares one coordinate -- so hexagon/pixel-level analysis and rendering see a "
                        "distribution instead of a grid. Pick roughly half the bin pitch (0.8 for 2um bins). "
                        "Like the ingest feature filters this rewrites the data itself, so every later stage "
                        "sees the jittered coordinates. Supported on every platform whose ingest runs "
                        "sge_convert (all but cosmx_smi); ignored (with a warning) for samples that supply an "
                        "already-ingested transcript. Equivalent config key: {\"ingest\": {\"jitter_xy\": 0.8}}.")
    j.add_argument("--units-per-um", type=float, default=None,
                   help="How many coordinate units of the raw input make one micron, i.e. the "
                        "factor sge_convert divides the input X/Y by (1 = the input is already "
                        "in microns, 1000 = nanometers, 2 = Stereo-seq's 500nm bins). Only needed "
                        "when a dataset departs from its platform's convention -- notably Illumina "
                        "StrataMap, whose older barcodes encode nanometers ('SBC:433503:2393851', "
                        "--units-per-um 1000) and whose current ones encode microns "
                        "('SBC:686.951:4668.15', --units-per-um 1); left unset, ingest detects "
                        "which of the two a barcode file uses. Applies to every platform whose "
                        "ingest runs sge_convert (all but cosmx_smi); ignored (with a warning) for "
                        "samples that supply an already-ingested transcript. Equivalent config key: "
                        "{\"ingest\": {\"units_per_um\": 1}}.")

    d.add_argument("--min-ct-per-unit-hexagon", type=int, default=None,
                   help=f"Minimum count per hexagon for FICTURE (default: {DEFAULT_MIN_CT_PER_UNIT_HEXAGON})")
    d.add_argument("--min-ct-per-unit-train", type=int, default=None,
                   help="Minimum count per hexagon during LDA training, counted over the features the "
                        "analysis uses (default: FICTURE2's own). Worth lowering when a feature list "
                        "restricts the analysis to a small panel, since --min-ct-per-unit-hexagon is "
                        "applied earlier, over all genes.")
    d.add_argument("--cell-min-cell-count", type=int, default=None,
                   help="Minimum count per cell in the cell-based analysis, counted over the features it "
                        "uses (default: run_ficture2_multi_cells' own). Same rationale as "
                        "--min-ct-per-unit-train.")
    d.add_argument("--cell-min-feature-count", type=int, default=None,
                   help="Minimum total count per feature in the cell-based analysis "
                        "(default: run_ficture2_multi_cells' own)")
    d.add_argument("--always-single-molecule", action="store_true",
                   help="Force single-molecule ON for both pixel FICTURE and cell decode")
    d.add_argument("--never-single-molecule", action="store_true",
                   help="Force single-molecule OFF for both pixel FICTURE and cell decode "
                        "(default: ON for pixel FICTURE, OFF for cell decode)")
    d.add_argument("--colname-cell", type=str, default=None,
                   help="Name of the cell-id column in the transcript TSV (generic platform). "
                        "Providing it turns on cell_id-based cell clustering; the column must be "
                        "the 5th TSV column (X, Y, gene, count, cell_id). Mutually exclusive with "
                        "mex inputs. Omit it (and provide no mex inputs) to skip cell analysis.")

    k = p.add_argument_group("Packaging overrides (else profile / run_cartload2 defaults)")
    k.add_argument("--bin-count", type=int, default=None,
                   help="Number of gene bins for the point PMTiles layers (profile default: "
                        "500 on most platforms)")

    pub = p.add_argument_group("Publish (opt-in; enable with --anno and/or --s3-upload)")
    pub.add_argument("--anno", action="store_true", help="AI-annotate each sample (requires --tissue and --organism)")
    pub.add_argument("--s3-upload", action="store_true", help="Upload each sample to S3 (requires --collection)")
    # annotation (tissue/organism required; the rest default)
    pub.add_argument("--tissue", type=str, default=None, help="Tissue for --anno (required for --anno)")
    pub.add_argument("--organism", type=str, default=None, help="Organism/species for --anno (required for --anno)")
    pub.add_argument("--anno-api-type", type=str, default="umgpt", help="AI annotation API type (default: umgpt)")
    pub.add_argument("--anno-model", type=str, default="claude-opus-4-7", help="AI annotation model (default: claude-opus-4-7)")
    pub.add_argument("--anno-threads", type=int, default=10, help="AI annotation threads (default: 10)")
    # S3 upload
    pub.add_argument("--collection", type=str, default=None, help="Collection name (required for --s3-upload)")
    pub.add_argument("--batch", type=str, default=None, help="Batch segment (default: current YYYY_MM)")
    pub.add_argument("--s3-prefix", type=str, default=DEFAULT_S3_PREFIX, help=f"S3 destination prefix (default: {DEFAULT_S3_PREFIX})")
    pub.add_argument("--aws-profile", type=str, default=DEFAULT_S3_PROFILE, help=f"AWS CLI profile (default: {DEFAULT_S3_PROFILE})")
    pub.add_argument("--aws", type=str, default="aws", help="Path to the aws CLI binary (default: aws)")
    pub.add_argument("--s3-jobs", type=int, default=4, help="Parallel S3 copies (xargs -P) per upload (default: 4)")
    return p.parse_args(_args)


def run_together(_args):
    args = parse_arguments(_args)
    cfg = build_config(args)

    # Publish is opt-in per action. Each action requires its mandatory inputs.
    if args.anno and not (args.tissue and args.organism):
        sys.exit("ERROR: --anno requires --tissue and --organism (no defaults).")
    # Collection defaults to the run id, matching the <multi_id>-<sample_id>
    # per-sample naming; override with --collection. This only affects the S3
    # destination (<s3_prefix>/batch=/<collection>/<dir>), not the local layout.
    #
    # For a single-sample --out-dir run the out-dir basename names the RUN, and the
    # per-sample dir is <out-dir basename>-<sample id>; using the out-dir basename as
    # the collection too would produce a redundant <name>/<name>-rep1, so its PARENT
    # dir is the collection. A joint --out-dir run (>1 sample) or an --out-root batch
    # keeps the out-dir/out-root basename as the collection.
    if args.s3_upload and not args.collection:
        single_outdir = bool(cfg.get("out_dir")) and len(cfg["_raw_samples"]) == 1
        if single_outdir:
            args.collection = os.path.basename(os.path.dirname(os.path.abspath(cfg["out_dir"])))
        else:
            args.collection = os.path.basename(os.path.normpath(cfg.get("out_dir") or cfg.get("out_root") or "."))
    if not args.batch:
        args.batch = datetime.date.today().strftime("%Y_%m")

    samples = [resolve_sample(raw, cfg) for raw in cfg["_raw_samples"]]

    out_dirs = sorted({s["out_dir"] for s in samples})
    anchor = cfg.get("out_root") or out_dirs[0]
    os.makedirs(anchor, exist_ok=True)

    mm = minimake()
    add_targets(mm, samples, cfg, args)
    if not mm.targets:
        sys.exit("ERROR: no targets generated. Check --only/--skip and your configuration.")

    # Prune prereqs to emitted flags so excluded upstream stages don't break make.
    valid = set(mm.targets.keys())
    for tgt, (srcs, cmds) in mm.targets.items():
        mm.targets[tgt] = ([s for s in srcs if s in valid], cmds)

    resolved = {k: v for k, v in cfg.items() if not k.startswith("_")}
    resolved["samples"] = samples
    with open(os.path.join(anchor, "run_together.resolved.json"), "w") as f:
        json.dump(resolved, f, indent=2, default=str)

    make_f = os.path.join(anchor, args.makefn)
    mm.write_makefile(make_f)
    print(f"Wrote master Makefile: {make_f}", flush=True)
    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)


if __name__ == "__main__":
    run_together(sys.argv[1:])
