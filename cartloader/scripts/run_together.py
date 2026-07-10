import sys, os, argparse, inspect, json, copy, csv

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import execute_makefile

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROFILE_DIR = os.path.join(repo_dir, "assets", "run_together_profiles")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def deep_merge(base, override):
    """Recursively merge ``override`` into ``base`` (returns a new dict).

    Dict values are merged key-by-key; every other type (including lists) is
    replaced wholesale by the override value.
    """
    if not isinstance(base, dict) or not isinstance(override, dict):
        return copy.deepcopy(override)
    out = copy.deepcopy(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_json(path):
    with open(path) as f:
        return json.load(f)


def load_builtin_profile(platform):
    path = os.path.join(PROFILE_DIR, f"{platform}.json")
    if not os.path.exists(path):
        avail = sorted(
            os.path.splitext(f)[0] for f in os.listdir(PROFILE_DIR) if f.endswith(".json")
        )
        sys.exit(
            f"ERROR: no built-in profile for platform '{platform}'. "
            f"Available: {', '.join(avail)}. Provide one via --profile."
        )
    return load_json(path)


def first_existing(in_dir, candidates):
    """Return the first path under ``in_dir`` that exists, else None."""
    for rel in candidates:
        p = os.path.join(in_dir, rel)
        if os.path.exists(p):
            return p
    return None


def infer_meta_from_outdir(out_dir):
    """Infer (id, collection, batch) from a `.../batch/collection/id` layout."""
    out_dir = os.path.normpath(out_dir)
    parts = out_dir.split(os.sep)
    sid = parts[-1] if parts else out_dir
    collection = parts[-2] if len(parts) >= 2 else None
    batch = parts[-3] if len(parts) >= 3 else None
    return sid, collection, batch


def read_sheet(path):
    """Read a TSV sample-sheet into a list of dicts (one per row).

    A value of '-' (or empty) means "unset / use default".
    """
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            clean = {}
            for k, v in row.items():
                if k is None:
                    continue
                v = (v or "").strip()
                if v and v != "-":
                    clean[k.strip()] = v
            if clean:
                rows.append(clean)
    return rows


# ---------------------------------------------------------------------------
# Configuration resolution
# ---------------------------------------------------------------------------

def build_config(args):
    """Assemble a normalized run configuration from CLI / --config / --sheet."""
    cfg = {}

    if args.config:
        cfg = load_json(args.config)

    # Overlay TSV sheet (tier 3): each row becomes/updates a sample entry.
    if args.sheet:
        cfg.setdefault("samples", [])
        for row in read_sheet(args.sheet):
            if "id" not in row:
                sys.exit("ERROR: --sheet must include an 'id' column")
            cfg["samples"].append(row)

    # Overlay CLI-only single-sample mode (tier 1).
    if args.in_dir:
        cfg.setdefault("samples", [])
        sample = {"in_dir": args.in_dir}
        if args.id:
            sample["id"] = args.id
        cfg["samples"].append(sample)

    # Top-level CLI overrides
    if args.platform:
        cfg["platform"] = args.platform
    if args.out_dir:
        cfg["out_dir"] = args.out_dir
    if args.out_root:
        cfg["out_root"] = args.out_root
    if args.n_factor:
        cfg.setdefault("defaults", {}).setdefault("ficture", {})["n_factor"] = args.n_factor
    if args.width:
        cfg.setdefault("defaults", {}).setdefault("ficture", {})["width"] = args.width
    if args.pretrained_model:
        cfg.setdefault("defaults", {}).setdefault("ficture", {})["pretrained_model"] = args.pretrained_model

    if not cfg.get("samples"):
        sys.exit("ERROR: no samples defined. Use --in-dir, --config, or --sheet.")

    # Resolve out_dir (single) vs out_root (batch of many samples).
    if not cfg.get("out_dir") and not cfg.get("out_root"):
        sys.exit("ERROR: --out-dir (or --out-root for --sheet batches) is required.")

    cfg.setdefault("resources", {})
    cfg["resources"].setdefault("n_jobs", args.n_jobs)
    cfg["resources"].setdefault("threads", args.threads)

    return cfg


def resolve_sample(raw, cfg, external_profile):
    """Produce the fully-resolved settings for a single sample."""
    platform = raw.get("platform") or cfg.get("platform")
    if not platform:
        sys.exit(f"ERROR: sample {raw.get('id', '?')} has no platform (set --platform or a 'platform' field).")

    # Merge chain: builtin -> external --profile -> inline cfg 'profile' -> defaults -> sample
    prof = load_builtin_profile(platform)
    if external_profile:
        prof = deep_merge(prof, external_profile)
    if cfg.get("profile"):
        prof = deep_merge(prof, cfg["profile"])
    if cfg.get("defaults"):
        prof = deep_merge(prof, cfg["defaults"])
    # sample-level overrides limited to known override subsections + scalars
    sample_overrides = {k: v for k, v in raw.items()
                        if k in ("ficture", "cartload", "annotate", "exclude_feature_regex")}
    # Convenience: flat sample-level scalars (handy as TSV sheet columns) map
    # into the FICTURE subsection.
    flat_ficture = {k: raw[k] for k in ("pretrained_model", "model_id", "n_factor", "width")
                    if k in raw}
    if flat_ficture:
        sample_overrides.setdefault("ficture", {})
        sample_overrides["ficture"] = deep_merge(sample_overrides["ficture"], flat_ficture)
    prof = deep_merge(prof, sample_overrides)

    prof["platform"] = platform

    # Determine per-sample id and out_dir.
    if cfg.get("out_dir") and len(cfg["samples"]) == 1 and "id" not in raw:
        sid, _, _ = infer_meta_from_outdir(cfg["out_dir"])
    else:
        sid = raw.get("id")
    if not sid:
        sid = os.path.basename(os.path.normpath(raw["in_dir"]))

    if cfg.get("out_dir"):
        out_dir = cfg["out_dir"]
    else:
        out_dir = os.path.join(cfg["out_root"], sid)

    prof["id"] = sid
    prof["in_dir"] = raw.get("in_dir")
    prof["out_dir"] = out_dir
    # Separate H&E: the profile 'hne' holds settings (a dict); the sample 'hne'
    # holds the actual image path. Keep them under distinct keys.
    if isinstance(prof.get("hne"), dict):
        prof["_hne_settings"] = prof.pop("hne")
    prof["hne"] = raw.get("hne")
    return prof


# ---------------------------------------------------------------------------
# Command builders (one per stage) -> return list[str] of shell commands
# ---------------------------------------------------------------------------

def cmd_sge_convert(s, sge_dir, res, exclude_regex):
    ing = s.get("ingest", {})
    parts = [
        "cartloader", "sge_convert",
        f"--platform {ing.get('sge_platform', s['platform'])}",
        f"--out-dir {sge_dir}",
        f"--n-jobs {res['n_jobs']}",
        f"--pigz-threads {res['threads']}",
        "--gzip pigz",
    ]
    if exclude_regex:
        parts.append(f"--exclude-feature-regex \"{exclude_regex}\"")

    # Input resolution: either autodetect (single file -> flag) or fixed inputs.
    if "autodetect" in ing:
        chosen = None
        for cand in ing["autodetect"]:
            p = os.path.join(s["in_dir"], cand["file"])
            if os.path.exists(p):
                chosen = (cand["flag"], p)
                break
        if chosen is None:
            # Fall back to first candidate path so the (missing) prereq is explicit.
            cand = ing["autodetect"][0]
            chosen = (cand["flag"], os.path.join(s["in_dir"], cand["file"]))
        parts.append(f"{chosen[0]} {chosen[1]}")
    if "inputs" in ing:
        for flag, rel in ing["inputs"].items():
            parts.append(f"{flag} {os.path.join(s['in_dir'], rel)}")
    if ing.get("csv_colnames_others"):
        parts.append("--csv-colnames-others " + " ".join(ing["csv_colnames_others"]))
    parts.extend(ing.get("extra_flags", []))
    return " ".join(p for p in parts if p)


def cmd_ficture(in_list, fic_dir, s, res, exclude_regex):
    fic = s.get("ficture", {})
    parts = [
        "cartloader", "run_ficture2_multi",
        f"--in-list {in_list}",
        f"--out-dir {fic_dir}",
        f"--width {fic['width']}",
        f"--threads {res['threads']}",
        f"--n-jobs {res['n_jobs']}",
        "--gzip pigz",
    ]
    if fic.get("pretrained_model"):
        parts.append(f"--pretrained-model {fic['pretrained_model']}")
        if fic.get("model_id"):
            parts.append(f"--model-id {fic['model_id']}")
    else:
        parts.append(f"--n-factor {fic['n_factor']}")
    if fic.get("min_ct_per_unit_hexagon") is not None:
        parts.append(f"--min-ct-per-unit-hexagon {fic['min_ct_per_unit_hexagon']}")
    if fic.get("decode_scale"):
        parts.append(f"--decode-scale {fic['decode_scale']}")
    if fic.get("single_molecule"):
        parts.append("--single-molecule")
    if exclude_regex:
        parts.append(f"--exclude-feature-regex \"{exclude_regex}\"")
    return " ".join(parts)


def model_path(fic_dir, s):
    fic = s.get("ficture", {})
    width = fic["width"].split(",")[0]
    nfs = [int(x) for x in str(fic.get("n_factor", "0")).split(",") if x]
    largest = max(nfs) if nfs else fic.get("n_factor")
    return os.path.join(fic_dir, f"t{width}_f{largest}.model.tsv")


def cmd_ficture_cells(prefix, list_flags, fic_dir, model, s, res, exclude_regex):
    parts = [
        "cartloader", "run_ficture2_multi_cells", "--all",
        f"--out-prefix {prefix}",
        f"--out-dir {fic_dir}",
        f"--threads {res['threads']}",
        f"--n-jobs {res['n_jobs']}",
        f"--pretrained-model {model}",
        "--gzip pigz",
    ]
    parts.extend(list_flags)
    if s.get("ficture", {}).get("single_molecule"):
        parts.append("--single-molecule")
    if exclude_regex:
        parts.append(f"--exclude-feature-regex \"{exclude_regex}\"")
    return " ".join(parts)


def cmd_cartload(fic_sample_dir, cart_dir, s, res, cell_params):
    cart = s.get("cartload", {})
    parts = [
        "cartloader", "run_cartload2",
        f"--fic-dir {fic_sample_dir}",
        f"--out-dir {cart_dir}",
        f"--id {s['id']}",
        f"--n-jobs {res['n_jobs']}",
        f"--threads {res['threads']}",
        "--colname-count count",
        "--gzip pigz",
    ]
    if cell_params:
        parts.append("--in-cell-params " + " ".join(cell_params))
    if cart.get("use_pmpoint"):
        parts.append("--use-pmpoint")
    if cart.get("sge_scale"):
        parts.append(f"--sge-scale {cart['sge_scale']}")
    if cart.get("bin_count") is not None:
        parts.append(f"--bin-count {cart['bin_count']}")
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Makefile assembly
# ---------------------------------------------------------------------------

def add_sample_targets(mm, samples, cfg, args):
    """Add all per-stage targets for every sample. Joint FICTURE fans in/out."""
    res = cfg["resources"]
    stages = set(args.only.split(",")) if args.only else None
    skip = set(args.skip.split(",")) if args.skip else set()

    def enabled(stage):
        if stage in skip:
            return False
        if stages is not None and stage not in stages:
            return False
        return True

    # All samples in one run share out_dir (single) or are independent (batch).
    # Group by out_dir so each out_dir gets one joint FICTURE model.
    groups = {}
    for s in samples:
        groups.setdefault(s["out_dir"], []).append(s)

    for out_dir, grp in groups.items():
        sge_root = os.path.join(out_dir, "tsv")
        fic_dir = os.path.join(out_dir, "fic")
        cart_root = os.path.join(out_dir, "cartl", "samples")
        mkdir = os.path.join(out_dir, "mk")
        os.makedirs(mkdir, exist_ok=True)
        os.makedirs(sge_root, exist_ok=True)

        # --- ingest (per sample) ---
        sge_flags = []
        sge_tsvs = {}
        for s in grp:
            sge_dir = os.path.join(sge_root, s["id"])
            tsv = os.path.join(sge_dir, "transcripts.unsorted.tsv.gz")
            sge_tsvs[s["id"]] = tsv
            flag = os.path.join(mkdir, f"sge.{s['id']}.done")
            sge_flags.append(flag)
            if enabled("ingest"):
                regex = s.get("exclude_feature_regex")
                cmds = [
                    f"mkdir -p {sge_dir}",
                    cmd_sge_convert(s, sge_dir, res, regex),
                    f"touch {flag}",
                ]
                mm.add_target(flag, [], cmds)

        # Write the joint in_list for FICTURE.
        in_list = os.path.join(sge_root, "in_list.tsv")
        with open(in_list, "w") as f:
            for s in grp:
                f.write(f"{s['id']}\t{sge_tsvs[s['id']]}\n")

        s0 = grp[0]  # FICTURE params come from the shared/first sample
        regex0 = s0.get("exclude_feature_regex")

        # --- ficture (joint model; fan-in on all ingests) ---
        fic_flag = os.path.join(mkdir, "ficture.done")
        if enabled("ficture"):
            mm.add_target(fic_flag, list(sge_flags), [
                cmd_ficture(in_list, fic_dir, s0, res, regex0),
                f"touch {fic_flag}",
            ])

        # --- ficture cells (segmentation-based decode; joint over samples) ---
        cells_flag = os.path.join(mkdir, "cells.done")
        cell_prefixes = plan_cells(grp, sge_root, s0)
        if cell_prefixes and enabled("cells"):
            model = model_path(fic_dir, s0)
            cmds = []
            for prefix, list_flags in cell_prefixes:
                cmds.append(cmd_ficture_cells(prefix, list_flags, fic_dir, model, s0, res, regex0))
            cmds.append(f"touch {cells_flag}")
            mm.add_target(cells_flag, [fic_flag], cmds)

        # --- cartload + imports (per sample; fan-out) ---
        cart_prereq = cells_flag if (cell_prefixes and enabled("cells")) else fic_flag
        for s in grp:
            fic_sample_dir = os.path.join(fic_dir, "samples", s["id"])
            cart_dir = os.path.join(cart_root, s["id"])
            cart_flag = os.path.join(mkdir, f"cartload.{s['id']}.done")
            cell_params = cell_param_paths(fic_sample_dir, cell_prefixes)
            if enabled("cartload"):
                mm.add_target(cart_flag, [cart_prereq], [
                    f"mkdir -p {cart_dir}",
                    cmd_cartload(fic_sample_dir, cart_dir, s, res, cell_params),
                    f"touch {cart_flag}",
                ])

            # --- images / squares / cells imports ---
            img_flag = os.path.join(mkdir, f"images.{s['id']}.done")
            img_cmds = plan_images(s, cart_dir, res)
            if img_cmds and enabled("images"):
                mm.add_target(img_flag, [cart_flag], img_cmds + [f"touch {img_flag}"])

            # --- publish (opt-in) ---
            if args.publish and cfg.get("publish") and enabled("publish"):
                pub_flag = os.path.join(mkdir, f"publish.{s['id']}.done")
                prereq = img_flag if (img_cmds and enabled("images")) else cart_flag
                mm.add_target(pub_flag, [prereq],
                              plan_publish(s, cart_dir, cfg["publish"]) + [f"touch {pub_flag}"])


def plan_cells(grp, sge_root, s0):
    """Write per-list TSVs and return [(prefix, [--list-* flags]), ...]."""
    cells = s0.get("cells")
    if not cells:
        return []
    prefixes = []

    def write_list(name, key, colkey=None):
        rows = []
        for s in grp:
            spec = s.get("cells", {}).get(key)
            if not spec:
                continue
            p = os.path.join(s["in_dir"], spec["file"])
            rows.append((s["id"], p))
        if not rows:
            return None
        path = os.path.join(sge_root, name)
        with open(path, "w") as f:
            for sid, p in rows:
                f.write(f"{sid}\t{p}\n")
        return path

    boundaries = write_list("in_boundaries.tsv", "boundaries")
    xy = write_list("in_xy.tsv", "xy")
    clusters = write_list("in_clust.tsv", "clusters")

    xy_spec = cells.get("xy", {})
    xy_flags = []
    if xy:
        xy_flags = [f"--list-xy {xy}",
                    f"--xy-colname-x {xy_spec.get('colname_x', 'X')}",
                    f"--xy-colname-y {xy_spec.get('colname_y', 'Y')}"]

    # Default "cartloader" prefix: boundaries (+ xy).
    if boundaries or xy:
        flags = []
        if boundaries:
            flags.append(f"--list-boundaries {boundaries}")
        flags.extend(xy_flags)
        prefixes.append(("cartloader", flags))

    # Named cluster-based prefix (e.g. xeniumranger).
    if clusters:
        prefix = cells["clusters"].get("prefix", "clusters")
        flags = [f"--list-cluster {clusters}"]
        if boundaries:
            flags.append(f"--list-boundaries {boundaries}")
        flags.extend(xy_flags)
        prefixes.append((prefix, flags))

    return prefixes


def cell_param_paths(fic_sample_dir, cell_prefixes):
    return [os.path.join(fic_sample_dir, f"ficture.{prefix}.params.json")
            for prefix, _ in cell_prefixes]


def plan_images(s, cart_dir, res):
    """Return shell commands importing morphology images / H&E / squares / cells."""
    cmds = []
    catalog = os.path.join(cart_dir, "catalog.yaml")

    # Morphology images (Xenium): first matching spec per id.
    seen = set()
    for spec in s.get("images", []):
        if spec["id"] in seen:
            continue
        img = first_existing(s["in_dir"], [spec["match"]])
        if not img:
            continue
        seen.add(spec["id"])
        extra = " ".join(spec.get("extra_flags", []))
        method = spec.get("method", "ome2png")
        conv = "--ome2png " if method == "ome2png" else ""
        cmds.append(
            f"cartloader import_image {conv}--png2pmtiles --georeference "
            f"--in-img {img} --out-dir {cart_dir} --img-id {spec['id']} "
            f"--upper-thres-quantile 0.95 --level 0 --colorize {spec['color']} "
            f"--transparent-below 5 {extra}".strip()
        )
        cmds.append(f"echo '    {spec['id']}: {spec['id']}.pmtiles' >> {catalog}")

    # Visium HD square-bin imports.
    for sq in s.get("squares", []):
        sub = os.path.join(s["in_dir"], sq["subdir"])
        if not os.path.exists(sub):
            continue
        cmds.append(
            f"cartloader import_visiumhd_square --in-dir {sub} "
            f"--outprefix {os.path.join(cart_dir, sq['suffix'])} "
            f"--bin-size {sq['bin_size']} --update-catalog"
        )

    # Visium HD segmented-cell import.
    ci = s.get("cell_import")
    if ci and os.path.exists(os.path.join(s["in_dir"], ci["detect_dir"])):
        cmds.append(
            f"cartloader import_visiumhd_cell --in-dir {s['in_dir']} "
            f"--outprefix {os.path.join(cart_dir, ci['suffix'])} --all --update-catalog"
        )

    # H&E image (Visium HD): profile settings live in '_hne_settings', the image
    # path in 'hne'.
    hne_path = s.get("hne")
    settings = s.get("_hne_settings", {})
    if hne_path and settings:
        out_id = settings.get("out_id", "hne")
        prefix = os.path.join(cart_dir, out_id)
        upp = ""
        jrel = settings.get("um_per_pixel_json")
        if jrel:
            jpath = os.path.join(s["in_dir"], jrel)
            key = settings.get("um_per_pixel_key", "microns_per_pixel")
            cmds.append(f'UPP=$(jq -r ".{key}" {jpath})')
            upp = '--um-per-pixel "$UPP"'
        plain = "--georef-plain" if settings.get("georef_plain") else ""
        cmds.append(
            f"cartloader image_png2pmtiles --in-img {hne_path} --out-prefix {prefix} "
            f"--geotif2mbtiles --mbtiles2pmtiles --georeference {plain} {upp}".strip()
        )
        cmds.append(f"echo '    {out_id}: {out_id}.pmtiles' >> {catalog}")

    return cmds


def plan_publish(s, cart_dir, publish):
    cmds = []
    anno = publish.get("annotate")
    if anno:
        parts = [f"cartloader anno_cartostore_folder --cartl-dir {cart_dir}"]
        if anno.get("tissue"):
            parts.append(f"--tissue \"{anno['tissue']}\"")
        if anno.get("organism"):
            parts.append(f"--organism {anno['organism']}")
        for k in ("api_type", "model", "threads"):
            if anno.get(k):
                parts.append(f"--{k.replace('_', '-')} {anno[k]}")
        cmds.append(" ".join(parts))
    up = publish.get("upload")
    if up:
        s3 = up["s3_prefix"].rstrip("/")
        col = publish.get("collection", "")
        batch = publish.get("batch", "")
        dest = f"{s3}/batch={batch}/{col}/{s['id']}"
        catalog = os.path.join(cart_dir, "catalog.yaml")
        profile = f"--profile {up['profile']}" if up.get("profile") else ""
        aws = up.get("aws", "aws")
        cmds.append(
            "(grep -E \"\\.\" " + catalog + " | perl -lane 'print $F[$#F]' | sort | uniq; "
            "echo catalog.yaml;) | xargs -I {} " + aws + " s3 cp " + cart_dir + "/{} "
            + dest + "/{} " + profile
        )
    return cmds


# ---------------------------------------------------------------------------
# Argument parsing / entry point
# ---------------------------------------------------------------------------

def parse_arguments(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="End-to-end multi-platform pipeline (ingest -> FICTURE -> cartload -> images -> publish) "
                    "orchestrated as a single Makefile.",
    )
    run = parser.add_argument_group("Run Options")
    run.add_argument("--dry-run", action="store_true", help="Write the Makefile and print commands without executing")
    run.add_argument("--restart", action="store_true", help="Ignore existing outputs and re-run all steps (make -B)")
    run.add_argument("-j", "--n-jobs", type=int, default=1, help="Parallel jobs for make and per-subcommand (default: 1)")
    run.add_argument("--threads", type=int, default=4, help="Threads per job (default: 4)")
    run.add_argument("--makefn", type=str, default="run_together.mk", help="Master Makefile name (default: run_together.mk)")
    run.add_argument("--only", type=str, help="Comma-separated stages to run exclusively (ingest,ficture,cells,cartload,images,publish)")
    run.add_argument("--skip", type=str, help="Comma-separated stages to skip")
    run.add_argument("--publish", action="store_true", help="Enable the opt-in publish stage (requires a 'publish' config block)")

    io = parser.add_argument_group("Input/Output")
    io.add_argument("--config", type=str, help="Path to JSON run configuration (tier 2)")
    io.add_argument("--sheet", type=str, help="Path to TSV sample-sheet (tier 3)")
    io.add_argument("--profile", type=str, help="Path to an external JSON profile that overrides the built-in one")
    io.add_argument("--platform", type=str, help="Platform preset (tier 1), e.g. 10x_xenium, 10x_visium_hd")
    io.add_argument("--in-dir", type=str, help="Input directory for a single-sample tier-1 run")
    io.add_argument("--out-dir", type=str, help="Output directory (single run / shared joint model)")
    io.add_argument("--out-root", type=str, help="Output root under which each --sheet sample gets its own dir")
    io.add_argument("--id", type=str, help="Sample id for a single-sample tier-1 run")

    key = parser.add_argument_group("Common FICTURE overrides")
    key.add_argument("--n-factor", type=str, help="Comma-separated factor counts (overrides profile)")
    key.add_argument("--width", type=str, help="Hexagon width in um (overrides profile)")
    key.add_argument("--pretrained-model", type=str, help="Pretrained model TSV -> projection instead of de-novo LDA")

    return parser.parse_args(_args)


def run_together(_args):
    args = parse_arguments(_args)

    external_profile = load_json(args.profile) if args.profile else None
    cfg = build_config(args)

    samples = []
    for raw in cfg["samples"]:
        if not raw.get("in_dir"):
            sys.exit(f"ERROR: sample {raw.get('id', '?')} has no 'in_dir'.")
        samples.append(resolve_sample(raw, cfg, external_profile))

    # Determine where to write the master Makefile: the out-root for batches,
    # otherwise the (single) shared out_dir.
    out_dirs = sorted({s["out_dir"] for s in samples})
    anchor = cfg.get("out_root") or out_dirs[0]
    os.makedirs(anchor, exist_ok=True)

    mm = minimake()
    add_sample_targets(mm, samples, cfg, args)

    if not mm.targets:
        sys.exit("ERROR: no targets generated. Check --only/--skip and your configuration.")

    # Prune prerequisites to only flags that are actually emitted. When a stage
    # is excluded via --only/--skip, its flag is dropped from downstream prereqs
    # (the user asserts that stage's outputs already exist), so `make` won't fail
    # with "No rule to make target".
    valid = set(mm.targets.keys())
    for tgt, (srcs, cmds) in mm.targets.items():
        mm.targets[tgt] = ([s for s in srcs if s in valid], cmds)

    # Provenance: write the resolved configuration.
    resolved = {"resources": cfg["resources"], "samples": [
        {k: v for k, v in s.items() if not k.startswith("_")} for s in samples
    ]}
    with open(os.path.join(anchor, "run_together.resolved.json"), "w") as f:
        json.dump(resolved, f, indent=2)

    make_f = os.path.join(anchor, args.makefn)
    mm.write_makefile(make_f)
    print(f"Wrote master Makefile: {make_f}", flush=True)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)


if __name__ == "__main__":
    run_together(sys.argv[1:])
