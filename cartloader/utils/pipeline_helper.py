import sys, os, argparse, logging, subprocess, inspect, re, shlex
from pathlib import Path
import hashlib
from typing import Iterable, List, Optional

from cartloader.utils.utils import add_param_to_cmd, cmd_separator, run_command_w_preq, load_file_to_dict, assert_unique

#  =============
#   args
#  =============
def validate_general_args(parser, args):
    """Group general cross-flag validations in one place."""
    # Mutual exclusion: --run-ficture2 vs --import-ext-ficture2
    if args.run_ficture2 and args.import_ext_ficture2:
        parser.error("--run-ficture2 and --import-ext-ficture2 are mutually exclusive")

    # run_ficture2 requires FICTURE params
    if args.run_ficture2:
        if not args.width:
            parser.error("--width is required when --run-ficture2 is set")
        if not args.n_factor:
            parser.error("--n-factor is required when --run-ficture2 is set")

    # import_ext_ficture2 requires external FICTURE dir
    if args.import_ext_ficture2 and not args.ext_fic_dir:
        parser.error("--ext-fic-dir is required when --import-ext-ficture2 is set")

    # cartload2 requires dataset id
    if args.run_cartload2 and not args.id:
        parser.error("--id is required when --run-cartload2 is set")

    # import cells requires cell_id
    if args.import_cells and not args.cell_id:
        parser.error("--cell-id is required when --import-cells is set")

    # Upload to AWS requires bucket and id
    if args.upload_aws:
        if not args.s3_dir:
            parser.error("--s3-dir is required when --upload-aws is set")
        if not args.id:
            parser.error("--id is required when --upload-aws is set")

    # Upload to Zenodo requires a valid token file
    if args.upload_zenodo:
        if not (args.zenodo_token and os.path.exists(args.zenodo_token)):
            parser.error("--zenodo-token must point to an existing file when --upload-zenodo is set")

    # Validate ID format when provided
    if args.id is not None:
        if re.search(r"\s", args.id):
            parser.error("--id must not contain whitespace; use hyphens if needed")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.id):
            parser.error("--id must match [A-Za-z0-9][A-Za-z0-9_-]* (alphanumeric, hyphens/underscores)")

#  =============
#   images
#  =============

def validate_imageid_args(image_ids, ranger_assets, all_images=False, dry_run=True):
    # Determine image IDs when use --use-json and --all-images
    if all_images:
        print(" * --all-images are enabled", flush=True)
        if dry_run:
            print(" * Dry run: skip image detection; using default values --image-ids ", flush=True)
            assert image_ids, ("In dry-run mode with --all-images, please pass --image-ids explicitly.")
        else:
            print(f"* Detecting all existing images from {ranger_assets}")
            data = load_file_to_dict(ranger_assets)
            images_data = data.get("IMAGES", None)
            if images_data is None:
                raise ValueError("The IMAGES section is missing in the input assets JSON; rerun detect step or provide --image-ids explicitly")
            images_data = {k: v for k, v in images_data.items() if v is not None}
            image_ids = list(images_data.keys())
    else:
        assert image_ids, "--image-ids is required when --import-images"

    # image_ids
    assert_unique(image_ids, "--image-ids")
    image_ids_out = list(image_ids)
    return image_ids_out

def validate_imagecol_args(image_ids, image_colors):
    # Colors
    if image_colors:
        assert_unique(image_colors, "--image-colors", normalize=lambda c: str(c).lstrip('#').lower())
    else:
        default_colors = [        ## A list of colors that work both for dark and light backgrounds
            "#1f77b4",  # Blue
            "#ff7f0e",  # Orange
            "#2ca02c",  # Green
            "#d62728",  # Red
            "#9467bd",  # Purple
            "#8c564b",  # Brown
            "#e377c2",  # Pink
            "#7f7f7f",  # Gray
            "#bcbd22",  # Olive
            "#17becf",  # Cyan
        ]
        image_colors = default_colors[:len(image_ids)]
    assert len(image_ids) <= len(image_colors), "Please specify image colors more or equal to the number of color image IDs"
    image_colors_out = list(image_colors)
    return image_colors_out

def validate_imageloc_args(ranger_dir, tifs_loc, loc_label, ranger_dir_label):
    tif_paths = []
    assert tifs_loc, f"{loc_label} is required when --import-images is set with manual input mode. Otherwise, enable auto-detection mode."
    assert_unique(tifs_loc, loc_label)
    # all items, should exist:
    for tif_loc in tifs_loc:
        tif_path = os.path.join(ranger_dir, tif_loc)
        assert os.path.exists(tif_path), f"File not found: {tif_path} (provided by {ranger_dir_label} and {loc_label})"
        tif_paths.append(tif_path)
    return tif_paths


def resolve_image_plan(image_ids, image_colors, use_json, ranger_assets, img_paths):
    """
    Build an image import plan with IDs, colors, inputs, and per-image prerequisites.
    Returns (plan)
    plan: list of dicts with keys {id, color, img_path, prereq}
    """
    # Build plan
    plan = []
    for i, iid in enumerate(image_ids):
        if not use_json:
            img_path = img_paths[i]
        plan.append({
                'id': iid,
                'color': image_colors[i].replace('#','') if len(image_colors) > 0 else None,
                'use_json': use_json,
                'img_path': None if use_json else img_path,
                'ranger_assets': ranger_assets if use_json else None,
                'prereq': [ranger_assets] if use_json else [img_path]
            })
    return plan

def stage_import_images(cart_dir, args, image_plans, ome2png=False, transparent_below=None, georef_detect=None, update_catalog=False):

    for spec in image_plans:
        image_id = spec['id']
        image_color = spec['color']

        use_json = spec['use_json']

        print(f" * Image ID: {image_id}", flush=True)
        print(f" * Image color: {image_color}", flush=True)

        import_image_cmd = " ".join([
            "cartloader", "import_image",
            f"--in-json {spec.get('ranger_assets', None)}" if use_json else "",
            f"--in-img {spec.get('img_path', None)}" if not use_json else "",
            f"--out-dir {cart_dir}",
            f"--img-id {image_id}",
            # ome2png
            "--ome2png" if ome2png else "",
            f"--colorize \"{image_color}\"" if ome2png and image_color else "",
            f"--upper-thres-quantile 0.95" if ome2png else "",
            f"--lower-thres-quantile 0.5" if ome2png else "",
            f"--transparent-below {transparent_below}" if ome2png and transparent_below is not None else "",
            # png2pmtiles
            "--png2pmtiles",
            f"--georeference --georef-detect {georef_detect}" if georef_detect else "",
            # update catalog
            f"--update-catalog" if update_catalog else "",
            f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
        ])
        import_image_cmd = add_param_to_cmd(import_image_cmd, args, ["pmtiles", "gdaladdo"])
        import_image_cmd = add_param_to_cmd(import_image_cmd, args, ["restart","n_jobs"])

        run_command_w_preq(import_image_cmd, prerequisites=spec.get('prereq', []), dry_run=args.dry_run, flush=True)

#  =============
#   squares
#  =============

def resolve_square_plan(ranger_dir, use_json, ranger_assets, square_input):
    square_plans = []

    if use_json:
        if not ranger_assets or not os.path.exists(ranger_assets):
            raise ValueError("No space_ranger assets JSON found; run with --load-space-ranger or provide --square-input entries")
        data = load_file_to_dict(ranger_assets)
        
        for key, spec in data.items():
            if not isinstance(key, str) and (key.startswith("GRID_") or key.startswith("SQUARE_")) and key.endswith("um"):
                continue
            try:
                bin_token = key.split("_", 1)[1]
                bin_size = int(bin_token.rstrip("um"))
            except (IndexError, ValueError):
                continue
            required = ["POSITION", "CLUSTER", "DE", "UMAP_PROJ"]
            missing = [field for field in required if not spec.get(field)]
            if missing:
                print(f" * Skipping {key} due to missing fields: {', '.join(missing)}", flush=True)
                continue
            square_plans.append({
                "bin_size": bin_size,
                "use_json": True,
                "ranger_assets": ranger_assets,
                "prereq": [ranger_assets],
            })
            print(f" * Detected square bin size {bin_size}um from {key}", flush=True)
    elif square_input:
        if not ranger_dir:
            raise ValueError("--space-ranger-dir is required when --square-input is provided")
        for square_input in square_input:
            inputs = [tok.strip() for tok in square_input.split(",")]
            if len(inputs) != 6:
                raise ValueError(f"Invalid --square-input entry: {square_input}. Must provide 6 comma-separated values: bin_size,position_parquet,csv_cluster,csv_diffexp,scale_json,csv_umap")
            bin_size_str, bin_pos_parquet, bin_csv_cluster, bin_csv_diffexp, bin_scale_json, bin_csv_umap = inputs
            try:
                bin_size = int(bin_size_str)
            except ValueError as exc:
                raise ValueError(f"Invalid bin size in --square-input entry: {square_input}") from exc
            prereq_paths = [bin_pos_parquet, bin_csv_cluster, bin_csv_diffexp, bin_scale_json, bin_csv_umap]
            prereq = [os.path.join(ranger_dir, f) for f in prereq_paths if f]
            square_plans.append({
                "bin_size": bin_size,
                "use_json": False,
                "ranger_assets": None,
                "prereq": prereq,
                "ranger_dir": ranger_dir,
                "POSITION": bin_pos_parquet,
                "CLUSTER": bin_csv_cluster,
                "DE": bin_csv_diffexp,
                "SCALE": bin_scale_json,
                "UMAP_PROJ": bin_csv_umap,
            })
            # print(f" * Added manual square bin size {bin_size}um from --square-input entry", flush=True)
    else:
        raise ValueError("No square binning input provided; use --load-space-ranger or provide --square-input entries")
    
    if not square_plans:
        raise ValueError("No square bin inputs detected; provide --square-input entries or regenerate --space-ranger-assets.")

    return square_plans


def stage_import_squares(cart_dir, args, square_plans, update_catalog=False):
    for spec in square_plans:   
        ranger_assets = spec.get("ranger_assets", None)         
        ranger_dir = spec.get("ranger_dir", None)
        
        bin_size = spec["bin_size"]
        outprefix = os.path.join(cart_dir, f"{args.square_id}-sq{bin_size:03d}")
        square_assets = f"{outprefix}_assets.json" # output 
        
        import_square_cmd = " ".join(filter(None, [
            "cartloader", "import_visiumhd_square",
            f"--outprefix {outprefix}",
            f"--bin-size {bin_size}",
            f"--in-json {ranger_assets}" if spec["use_json"] else "",
            f"--in-dir {ranger_dir}" if not spec["use_json"] else "",
            f"--parquet {spec.get('POSITION')}" if not spec["use_json"] and spec.get('POSITION') else "",
            f"--csv-clust {spec.get('CLUSTER')}" if not spec["use_json"] and spec.get('CLUSTER') else "",
            f"--csv-diffexp {spec.get('DE')}" if not spec["use_json"] and spec.get('DE') else "",
            f"--scale-json {spec.get('SCALE')}" if not spec["use_json"] and spec.get('SCALE') else "",
            f"--csv-umap {spec.get('UMAP_PROJ')}" if not spec["use_json"] and spec.get('UMAP_PROJ') else "",
            f"--tippecanoe {args.tippecanoe}" if args.tippecanoe else "",
            f"--tsv-cmap {args.tsv_cmap}" if args.tsv_cmap else "",
            f"--update-catalog" if update_catalog else "",
            f"--threads {args.threads}" if args.threads else "",
            f"--use-parquet-tools" if args.use_parquet_tools else "",
            f"--parquet-tools" if args.use_parquet_tools and args.parquet_tools else "",
            f"--pigz {args.pigz}" if args.pigz else "",
            f"--pigz-threads {args.pigz_threads}" if args.pigz and args.pigz_threads else "",
        ]))

        import_square_cmd = add_param_to_cmd(import_square_cmd, args, ["tippecanoe", "python", "spatula", "ficture2", "gzip", "sort", "R"])

        prereq = spec.get("prereq", [])
        if os.path.exists(square_assets) and not args.restart:
            print(f" * Skip --import-squares for bin size {bin_size} since the assets file ({square_assets}) already exists. Use --restart to force execution.", flush=True)
            # print("\n", flush=True)
            print(import_square_cmd, flush=True)
            continue

        run_command_w_preq(import_square_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)

#  =============
#   ficture
#  =============
def define_ficture_flags(fic_dir, width, n_factor, anchor_res: int = 6):
    flags = []
    for tw in width.split(","):
        for nf in n_factor.split(","):
            fw = tw
            dc_prefix = os.path.join(fic_dir, f"t{tw}_f{nf}_p{fw}_a{anchor_res}")
            flags.append(f"{dc_prefix}.done")  # decoding step
            flags.append(f"{dc_prefix}_summary.done")  # decoding info and de
            flags.append(f"{dc_prefix}.png")  # decoding visual
    return flags

def stage_run_ficture2(fic_dir, sge_assets, args, prereq, aux_args=None):
    print("="*10, flush=True)
    print("Executing --run-ficture2 (execute via make)", flush=True)
    print("="*10, flush=True)

    os.makedirs(fic_dir, exist_ok=True)

    print(f" * train width {args.width}", flush=True)
    print(f" * n factors {args.n_factor}", flush=True)
    
    # Build ficture2 command
    ficture2_cmd = " ".join([
        "cartloader", "run_ficture2",
        # actions
        "--main",
        f"--out-dir {fic_dir}",
        f"--in-json {sge_assets}",
        f"--width {args.width}",
        f"--n-factor {args.n_factor}",
        f"--colname-feature {args.colname_feature}" if args.colname_feature else "",
        f"--colname-count {args.colname_count}" if args.colname_count else "",
        f"--include-feature-regex \"{args.fic_include_feature_regex}\"" if args.fic_include_feature_regex else "",
        f"--exclude-feature-regex \"{args.fic_exclude_feature_regex}\"" if args.fic_exclude_feature_regex else "",
    ])
    
    # Add optional params
    ficture2_cmd = add_param_to_cmd(ficture2_cmd, args, ["spatula", "ficture2"])
    ficture2_cmd = add_param_to_cmd(ficture2_cmd, args, ["restart","n_jobs","threads"])
    if aux_args:
        ficture2_cmd += " " + aux_args

    # Execute
    run_command_w_preq(
        ficture2_cmd, 
        prerequisites=prereq,
        dry_run=getattr(args, "dry_run", False),
        flush=True,
    )

#  =============
#   cartload2
#  =============
def stage_run_cartload2(cart_dir, fic_dir, sge_dir, 
                        cell_assets, square_assets, background_assets, 
                        args, prereq):

    # Banner
    print("=" * 10, flush=True)
    print("Executing --run-cartload2 (execute via make; function: run_cartload2)", flush=True)
    print("=" * 10, flush=True)

    # Handle prerequisites based on args
    if args.run_ficture2:
        print(f" * Importing FICTURE results: {fic_dir}", flush=True)
        # Note: new runs may be added to existing runs
        # In make, a rule runs only if the target is missing or older than its prerequisites — changing the commands alone doesn’t trigger re-execution.
        # So, should be a flag file relevant to all requested runs
        ficture_flags = define_ficture_flags(fic_dir, args.width, args.n_factor, anchor_res=6)
        prereq.extend(ficture_flags)
        fic_params =os.path.join(fic_dir, "ficture.params.json")
        prereq.append(fic_params)
    elif args.import_ext_ficture2:
        print(f" * Importing external FICTURE results from {args.ext_fic_dir}", flush=True)
        fic_params = os.path.join(args.ext_fic_dir, "ficture.params.json")
        prereq.append(fic_params)

    if args.import_cells:
        prereq.append(cell_assets)
    if args.import_images:
        prereq.extend(background_assets)

    if getattr(args, "import_squares", False):
        prereq.extend(square_assets)

    # Build cartload2 command
    cartload_cmd = " ".join([
        "cartloader", "run_cartload2",
        f"--out-dir {cart_dir}",
        f"--fic-dir {fic_dir}" if args.run_ficture2 else "",
        f"--ext-fic-dir {args.ext_fic_dir}" if args.import_ext_ficture2 else "",
        f"--sge-dir {sge_dir}" if not (args.run_ficture2 or args.import_ext_ficture2) else "",
        f"--cell-assets {cell_assets}" if args.import_cells else "",
        f"--background-assets {' '.join(background_assets)}" if background_assets else "",
        f"--square-assets  {' '.join(square_assets)}" if square_assets else "",
        f"--id {args.id}",
        f"--title {shlex.quote(args.title)}" if args.title else "",
        f"--desc {shlex.quote(args.desc)}" if args.desc else "",
        f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
    ])

    # Add optional params
    cartload_cmd = add_param_to_cmd(cartload_cmd, args, ["pmtiles", "gdaladdo", "spatula", "tippecanoe"])
    cartload_cmd = add_param_to_cmd(cartload_cmd, args, ["restart", "n_jobs", "threads"])

    # Deduplicate prerequisites to avoid noisy repeats
    prereq = list(dict.fromkeys(prereq))
    # Execute
    run_command_w_preq(
        cartload_cmd,
        prerequisites=prereq,
        dry_run=getattr(args, "dry_run", False),
        flush=True,
    )

#  =============
#   aws
#  =============
def stage_upload_aws(cart_dir, args, prereq):
    # Banner
    print("=" * 10, flush=True)
    print("Executing --upload-aws (execute via make)", flush=True)
    print("=" * 10, flush=True)

    # Build command
    aws_cmd = " ".join([
        "cartloader", "upload_aws",
        f"--in-dir {cart_dir}",
        f"--s3-dir {args.s3_dir}",
    ])

    # Inject optional params (preserves your original calls)
    aws_cmd = add_param_to_cmd(aws_cmd, args, ["aws", "profile"])
    aws_cmd = add_param_to_cmd(aws_cmd, args, ["restart", "n_jobs"])
    
    # Execute
    run_command_w_preq(
        aws_cmd,
        prerequisites=prereq,
        dry_run=getattr(args, "dry_run", False),
        flush=True,
    )

#  =============
#   Zenodo
#  =============
def stage_upload_zenodo(cart_dir, args, prereq, flag):
    # Banner
    print("=" * 10, flush=True)
    print("Executing --upload-zenodo (execute directly)", flush=True)
    print("=" * 10, flush=True)

    # Set default title if needed
    if args.zenodo_title is None:
        args.zenodo_title = args.title

    # Build command
    zenodo_cmd = " ".join([
        "cartloader", "upload_zenodo",
        f"--in-dir {cart_dir}",
        f"--upload-method catalog",
        f"--zenodo-token {args.zenodo_token}",
        f"--zenodo-deposition-id {args.zenodo_deposition_id}" if args.zenodo_deposition_id else "",
        f"--title {shlex.quote(args.zenodo_title)}" if args.zenodo_title else "",
        ("--creators " + " ".join(shlex.quote(c) for c in args.creators)) if args.creators else "",
        f"--upload-type dataset" if not args.zenodo_deposition_id else "", # add type for new deposition
        f"--overwrite" if args.restart else ""
    ])

    # Skip if all flags exist and not restarting
    if (all(os.path.exists(f) for f in flag) and not args.restart):
        print(" * Skip --upload-zenodo since all flags exist. Use --restart to force execution.\n", flush=True)
        print(zenodo_cmd, flush=True)
    else:
        run_command_w_preq(zenodo_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)
