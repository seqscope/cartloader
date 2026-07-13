import sys, os, argparse, inspect, json, copy, csv, datetime

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import execute_makefile

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROFILE_DIR = os.path.join(repo_dir, "assets", "run_together_profiles")

# Global fallbacks (a profile or --config may override; a CLI flag wins over both).
DEFAULT_EXCLUDE_REGEX = "^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated|System|Gm[0-9]|MT-|mt-|Rps|Rpl|NCS-|NCP-)"
DEFAULT_MIN_CT_PER_UNIT_HEXAGON = 50

# Publish (S3 upload) defaults for the CartoStore project.
DEFAULT_S3_PREFIX = "s3://cartostore/data"
DEFAULT_S3_PROFILE = "default"

# Sample-level input roles. A sample provides these either explicitly (sample
# sheet columns / JSON) or by auto-detection inside its `in_dir` (profile.roles).
ROLE_KEYS = ["transcript", "xy", "boundaries", "clusters", "mex"]
# role -> the run_ficture2_multi_cells flag that consumes it
ROLE_LIST_FLAG = {
    "boundaries": "--list-boundaries",
    "xy": "--list-xy",
    "clusters": "--list-cluster",
    "mex": "--mex-list",
}

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


def infer_meta_from_outdir(out_dir):
    parts = os.path.normpath(out_dir).split(os.sep)
    return parts[-1] if parts else out_dir


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

def build_config(args):
    cfg = load_json(args.config) if args.config else {}

    platform = args.platform or cfg.get("platform")
    if not platform:
        sys.exit("ERROR: --platform (or a 'platform' field in --config) is required.")

    # Layer 0: profile
    prof = load_builtin_profile(platform)
    if args.platform_json:
        prof = deep_merge(prof, load_json(args.platform_json))
    prof["platform"] = platform
    prof.setdefault("ficture_defaults", {})
    prof.setdefault("ficture", [])
    prof.setdefault("cell_analyses", [])
    prof.setdefault("images", [])
    prof.setdefault("roles", {})

    # Layer 1: tier-1 CLI selects the base FICTURE mode
    if args.project_models:
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
    for k in ("exclude_feature_regex", "ingest", "roles", "cartload", "ficture_defaults",
              "squares", "cell_import", "hne", "publish", "resources"):
        if k in cfg:
            prof[k] = deep_merge(prof.get(k), cfg[k]) if isinstance(cfg[k], dict) else cfg[k]
    if "ficture" in cfg:
        prof["ficture"] = merge_id_list(prof["ficture"], cfg["ficture"])
    if "cell_analyses" in cfg:
        prof["cell_analyses"] = merge_id_list(prof["cell_analyses"], cfg["cell_analyses"])
    if "images" in cfg:
        prof["images"] = merge_images(prof["images"], cfg["images"])

    # --- resolve common decode defaults (CLI > config/profile > hardcoded) ---
    # exclude-feature regex
    if args.exclude_feature_regex is not None:
        prof["exclude_feature_regex"] = args.exclude_feature_regex
    elif not prof.get("exclude_feature_regex"):
        prof["exclude_feature_regex"] = DEFAULT_EXCLUDE_REGEX
    # min count per unit hexagon (applies to the pixel FICTURE analyses)
    fd = prof.setdefault("ficture_defaults", {})
    if args.min_ct_per_unit_hexagon is not None:
        fd["min_ct_per_unit_hexagon"] = args.min_ct_per_unit_hexagon
    elif "min_ct_per_unit_hexagon" not in fd:
        fd["min_ct_per_unit_hexagon"] = DEFAULT_MIN_CT_PER_UNIT_HEXAGON
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

    # Apply ficture_defaults to every analysis (per-entry keys win).
    prof["ficture"] = [deep_merge(prof["ficture_defaults"], a) for a in prof["ficture"]]

    prof["out_dir"] = args.out_dir or cfg.get("out_dir")
    prof["out_root"] = args.out_root or cfg.get("out_root")
    prof.setdefault("resources", {})
    prof["resources"].setdefault("n_jobs", args.n_jobs)
    prof["resources"].setdefault("threads", args.threads)

    # Samples: --in-dir (one), --samples sheet (many), or JSON 'samples'.
    raw_samples = list(cfg.get("samples", []))
    if args.samples:
        raw_samples += read_sheet(args.samples)
    if args.in_dir:
        s = {"in_dir": args.in_dir}
        if args.id:
            s["id"] = args.id
        raw_samples.append(s)
    if not raw_samples:
        sys.exit("ERROR: no samples. Use --in-dir, --samples <sheet>, or a 'samples' config block.")
    prof["_raw_samples"] = raw_samples

    if not prof["out_dir"] and not prof["out_root"]:
        sys.exit("ERROR: --out-dir (or --out-root for independent per-sample models) is required.")
    return prof


def resolve_sample(raw, cfg):
    """Resolve one sample's id, out_dir, and input roles."""
    n = len(cfg["_raw_samples"])
    if cfg.get("out_dir") and n == 1 and "id" not in raw:
        sid = infer_meta_from_outdir(cfg["out_dir"])
    else:
        sid = raw.get("id") or (os.path.basename(os.path.normpath(raw["in_dir"])) if raw.get("in_dir") else None)
    if not sid:
        sys.exit(f"ERROR: cannot determine id for sample {raw}")

    out_dir = cfg["out_dir"] if cfg.get("out_dir") else os.path.join(cfg["out_root"], sid)
    in_dir = raw.get("in_dir")

    # Resolve roles: explicit column/JSON value wins; else auto-detect in in_dir.
    roles = {}
    role_specs = cfg.get("roles", {})
    for role in ROLE_KEYS:
        if raw.get(role):
            path = raw[role]
            if in_dir and not os.path.isabs(path) and not os.path.exists(path):
                path = os.path.join(in_dir, path)
            roles[role] = path
        elif role in role_specs and in_dir:
            cand = os.path.join(in_dir, role_specs[role]["file"])
            if os.path.exists(cand):
                roles[role] = cand
    return {
        "id": sid, "in_dir": in_dir, "out_dir": out_dir,
        "roles": roles,
        "hne": raw.get("hne"),
        "images": raw.get("images", []),
    }


# ---------------------------------------------------------------------------
# Command builders
# ---------------------------------------------------------------------------

def cmd_sge_convert(cfg, sge_dir, in_dir):
    ing = cfg.get("ingest", {})
    res = cfg["resources"]
    parts = ["cartloader", "sge_convert",
             f"--platform {ing.get('sge_platform', cfg['platform'])}",
             f"--out-dir {sge_dir}", f"--n-jobs {res['n_jobs']}",
             f"--pigz-threads {res['threads']}", "--gzip pigz"]
    if cfg.get("exclude_feature_regex"):
        parts.append(f"--exclude-feature-regex \"{cfg['exclude_feature_regex']}\"")
    if "autodetect" in ing:
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
        parts.append(f"{flag} {os.path.join(in_dir, rel)}")
    if ing.get("csv_colnames_others"):
        parts.append("--csv-colnames-others " + " ".join(ing["csv_colnames_others"]))
    parts.extend(ing.get("extra_flags", []))
    return " ".join(p for p in parts if p)


def cmd_ficture_analysis(a, in_list, fic_dir, cfg):
    res = cfg["resources"]
    parts = ["cartloader", "run_ficture2_multi",
             f"--in-list {in_list}", f"--out-dir {fic_dir}",
             f"--width {a['width']}", f"--threads {res['threads']}",
             f"--n-jobs {res['n_jobs']}", "--gzip pigz"]
    if a.get("mode") == "project":
        parts += [f"--pretrained-model {a['model']}", f"--model-id {a['id']}"]
    else:
        parts.append(f"--n-factor {a['n_factor']}")
    if a.get("min_ct_per_unit_hexagon") is not None:
        parts.append(f"--min-ct-per-unit-hexagon {a['min_ct_per_unit_hexagon']}")
    if a.get("decode_scale"):
        parts.append(f"--decode-scale {a['decode_scale']}")
    if cfg.get("_sm_pixel"):   # single-molecule for pixel FICTURE (default ON)
        parts.append("--single-molecule")
    if cfg.get("exclude_feature_regex"):
        parts.append(f"--exclude-feature-regex \"{cfg['exclude_feature_regex']}\"")
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
    if cfg.get("_sm_cells"):   # single-molecule for cell decode (default OFF)
        parts.append("--single-molecule")
    if cfg.get("exclude_feature_regex"):
        parts.append(f"--exclude-feature-regex \"{cfg['exclude_feature_regex']}\"")
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
        else:  # match relative to in_dir
            src = first_existing(os.path.join(s["in_dir"], spec["match"])) if s.get("in_dir") else None
            if not src:
                continue
        seen.add(iid)
        ops.append({**spec, "source": src})
    return ops


def plan_images(cfg, s, cart_dir, multi, transcript=None):
    cmds = []
    catalog = os.path.join(cart_dir, "catalog.yaml")

    for op in resolve_image_ops(cfg, s):
        iid, kind, src = op["id"], op.get("kind", "single"), op["source"]
        if kind == "prebuilt":
            cmds.append(f"cp {src} {os.path.join(cart_dir, iid + '.pmtiles')}")
        elif kind == "rgb":
            cmds += _cmd_rgb_image(cfg, s, iid, src, cart_dir, op)
            continue  # _cmd_rgb_image appends its own catalog line
        else:  # single-channel colorized
            conv = "--ome2png " if op.get("convert", "ome2png") == "ome2png" else ""
            extra = " ".join(op.get("extra_flags", []))
            cmds.append(
                f"cartloader import_image {conv}--png2pmtiles --georeference "
                f"--in-img {src} --out-dir {cart_dir} --img-id {iid} "
                f"--upper-thres-quantile 0.95 --level 0 --colorize {op['color']} "
                f"--transparent-below 5 {extra}".strip())
        cmds.append(f"echo '    {iid}: {iid}.pmtiles' >> {catalog}")

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
            # A sample can supply inputs via a Ranger-style --in-dir or via sheet
            # role columns (or both). Emit only when at least one is present.
            has_indir = bool(s.get("in_dir"))
            role_paths = {r: s["roles"][r] for r in ca.get("uses", []) if s["roles"].get(r)}
            if not (has_indir or role_paths):
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
    jrel = settings.get("um_per_pixel_json")
    if jrel and s.get("in_dir"):
        jpath = os.path.join(s["in_dir"], jrel)
        key = settings.get("um_per_pixel_key", "microns_per_pixel")
        # Inline the substitution: make runs each recipe line in its own shell,
        # so a UPP=... on a separate line would not survive to this command.
        upp = f"--um-per-pixel \"$(jq -r '.{key}' {jpath})\""
    plain = "--georef-plain" if settings.get("georef_plain") else ""
    cmds.append(f"cartloader image_png2pmtiles --in-img {src} --out-prefix {prefix} "
                f"--geotif2mbtiles --mbtiles2pmtiles --georeference {plain} {upp}".strip())
    cmds.append(f"echo '    {iid}: {iid}.pmtiles' >> {catalog}")
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
        sge_flags, transcript = [], {}
        for s in grp:
            if s["roles"].get("transcript"):
                transcript[s["id"]] = s["roles"]["transcript"]
                continue
            sge_dir = os.path.join(sge_root, s["id"])
            transcript[s["id"]] = os.path.join(sge_dir, "transcripts.unsorted.tsv.gz")
            flag = os.path.join(mkdir, f"sge.{s['id']}.done")
            sge_flags.append(flag)
            if on("ingest"):
                mm.add_target(flag, [], [f"mkdir -p {sge_dir}",
                                         cmd_sge_convert(cfg, sge_dir, s["in_dir"]),
                                         f"touch {flag}"])

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

        # --- cell analyses (platform default; run those whose roles are present) ---
        cells_flag = os.path.join(mkdir, "cells.done")
        active_cells = plan_cell_analyses(grp, sge_root, cfg, fic_dir, default_model_id, multi)
        if active_cells and on("cells"):
            cmds = [c["cmd"] for c in active_cells]
            mm.add_target(cells_flag, [fic_flag], cmds + [f"touch {cells_flag}"])

        cart_prereq = cells_flag if (active_cells and on("cells")) else fic_flag
        active_ids = [c["id"] for c in active_cells]
        multi_cart_flag = os.path.join(mkdir, "cartload.done")
        if multi and on("cartload"):
            mm.add_target(multi_cart_flag, [cart_prereq], [
                f"mkdir -p {cart_root}",
                cmd_cartload_multi(fic_dir, cart_root, multi_id, cfg),
                f"touch {multi_cart_flag}"])

        # --- cartload + images (per sample); collect each sample's post-images flag ---
        sample_ctx = []   # (sample, cart_dir, base_prereq)
        for s in grp:
            fic_sample_dir = os.path.join(fic_dir, "samples", s["id"])
            if multi:
                cart_dir = os.path.join(cart_root, f"{multi_id}-{s['id']}")
                cart_flag = multi_cart_flag
            else:
                cart_dir = os.path.join(cart_root, s["id"])
                cart_flag = os.path.join(mkdir, f"cartload.{s['id']}.done")
                cell_params = [os.path.join(fic_sample_dir, f"ficture.{cid}.params.json")
                               for cid in active_ids if _sample_in_cell(s, cfg, cid, active_cells)]
                if on("cartload"):
                    mm.add_target(cart_flag, [cart_prereq], [
                        f"mkdir -p {cart_dir}",
                        cmd_cartload(fic_sample_dir, cart_dir, s["id"], cfg, cell_params),
                        f"touch {cart_flag}"])

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


def _sample_in_cell(s, cfg, cid, active_cells):
    ca = next((a for a in cfg["cell_analyses"] if a["id"] == cid), None)
    if not ca:
        return False
    return all(s["roles"].get(r) for r in ca.get("uses", []))


def plan_cell_analyses(grp, sge_root, cfg, fic_dir, default_model_id, multi):
    """Write per-cell-analysis role lists; return the analyses that have inputs.

    An analysis with a ``multi_import`` command is a per-sample import for joint
    runs (e.g. Xenium Ranger clusters, which are sample-specific and cannot be
    jointly decoded) — it is skipped here for multi and handled in the image stage.
    """
    active = []
    for ca in cfg.get("cell_analyses", []):
        if multi and ca.get("multi_import"):
            continue
        uses = ca.get("uses", [])
        # samples that have every required role for this analysis
        contributing = [s for s in grp if all(s["roles"].get(r) for r in uses)]
        if not contributing:
            continue
        list_files = {}
        for role in uses:
            path = os.path.join(sge_root, f"in_{role}.{ca['id']}.tsv")
            with open(path, "w") as f:
                for s in contributing:
                    f.write(f"{s['id']}\t{s['roles'][role]}\n")
            list_files[role] = path
        model_id = ca.get("model_id", default_model_id)
        model_path = os.path.join(fic_dir, f"{model_id}.model.tsv")
        active.append({"id": ca["id"],
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

    io = p.add_argument_group("Input/Output")
    io.add_argument("--platform", type=str, help="Platform preset, e.g. 10x_xenium, 10x_visium_hd")
    io.add_argument("--in-dir", type=str, help="Input directory for a single sample")
    io.add_argument("--samples", type=str, help="TSV sample sheet (wide table of input roles) for a joint multi-sample run")
    io.add_argument("--out-dir", type=str, help="Output directory (single sample / shared joint model)")
    io.add_argument("--out-root", type=str, help="Output root; each sample gets its own dir and an independent model")
    io.add_argument("--id", type=str, help="Sample id for a single --in-dir run")
    io.add_argument("--config", type=str, help="JSON config that augments the profile/CLI (full spec for complex runs)")
    io.add_argument("--platform-json", type=str, help="External JSON profile that overrides the built-in platform profile")

    f = p.add_argument_group("FICTURE mode (choose one; default de-novo from the profile)")
    f.add_argument("--width", type=str, help="De-novo: hexagon width(s) in um (comma-separated)")
    f.add_argument("--n-factor", type=str, help="De-novo: factor count(s) (comma-separated)")
    f.add_argument("--project-models", type=str, help="Projection-only: existing FICTURE dir(s) (comma-separated); "
                                                      "reads each ficture.params.json and reuses its models. No LDA training.")

    d = p.add_argument_group("Common decode overrides (else profile / built-in defaults)")
    d.add_argument("--exclude-feature-regex", type=str, default=None,
                   help=f"Regex of features to exclude (default: profile's, else '{DEFAULT_EXCLUDE_REGEX}')")
    d.add_argument("--min-ct-per-unit-hexagon", type=int, default=None,
                   help=f"Minimum count per hexagon for FICTURE (default: {DEFAULT_MIN_CT_PER_UNIT_HEXAGON})")
    d.add_argument("--always-single-molecule", action="store_true",
                   help="Force single-molecule ON for both pixel FICTURE and cell decode")
    d.add_argument("--never-single-molecule", action="store_true",
                   help="Force single-molecule OFF for both pixel FICTURE and cell decode "
                        "(default: ON for pixel FICTURE, OFF for cell decode)")

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
    # Collection defaults to the run id (the out-dir basename), matching the
    # <multi_id>-<sample_id> per-sample naming; override with --collection.
    if args.s3_upload and not args.collection:
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
