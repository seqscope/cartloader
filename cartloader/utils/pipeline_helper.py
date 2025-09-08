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
        if not args.s3_bucket:
            parser.error("--s3-bucket is required when --upload-aws is set")
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

def validate_imageid_args(args, ranger_assets):
    # Determine image IDs when use --use-json and --all-images
    if args.all_images:
        print(" * --all-images are enabled", flush=True)
        if args.dry_run:
            print(" * Dry run: skip image detection; using default values --image-ids ", flush=True)
            assert args.image_ids, ("In dry-run mode with --all-images, please pass --image-ids explicitly.")
        else:
            print(f"* Detecting all existing images from {ranger_assets}")
            data = load_file_to_dict(ranger_assets)
            images_data = data.get("IMAGES", None)
            if images_data is None:
                raise ValueError("The IMAGES section is missing in the input assets JSON; rerun detect step or provide --image-ids explicitly")
            images_data = {k: v for k, v in images_data.items() if v is not None}
            args.image_ids = list(images_data.keys())
    else:
        assert args.image_ids, "--image-ids is required when --import-images"

    # image_ids
    assert_unique(args.image_ids, "--image-ids")
    image_ids = list(args.image_ids)
    print(f"    - image IDs (N={len(image_ids)}): {image_ids}")
    return image_ids

def validate_imagecol_args(args):
    # Colors
    if args.image_colors:
        assert_unique(args.image_colors, "--image-colors", normalize=lambda c: str(c).lstrip('#').lower())
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
        args.image_colors = default_colors[:len(args.image_ids)]
    assert len(args.image_ids) <= len(args.image_colors), "Please specify image colors more or equal to the number of color image IDs"
    print(f"    - image colors (N={len(args.image_colors)}): {args.image_colors}")
    image_colors = list(args.image_colors)
    return image_colors

def resolve_image_plan(args, ranger_assets, use_json, img_paths):
    """
    Build an image import plan with IDs, colors, inputs, and per-image prerequisites.
    Returns (plan)
    plan: list of dicts with keys {id, color, ome_path, prereq}
    """
    # Build plan
    plan = []
    for i, iid in enumerate(args.image_ids):
        if not use_json:
            img_path = img_paths[i]
        plan.append({
                'id': iid,
                'color': args.image_colors[i].replace('#',''),
                'ome_path': None if use_json else img_path,
                'prereq': [ranger_assets] if use_json else [img_path],
            })
    return plan

def stage_import_images(cart_dir, args, ranger_assets, use_json, tif_paths, catalog_yaml):
    plan = resolve_image_plan(args, ranger_assets, use_json, tif_paths)

    for spec in plan:
        image_id = spec['id']
        image_color = spec['color']
        prereq = spec['prereq']
        ome_path = spec['ome_path']

        print(f" * Image ID: {image_id}", flush=True)
        print(f" * Image color: {image_color}", flush=True)

        import_image_cmd = " ".join([
            "cartloader", "import_image",
            "--ome2png",
            "--png2pmtiles",
            f"--update-catalog" if not args.run_cartload2 and os.path.exists(catalog_yaml) else "",
            f"--in-json {ranger_assets}" if use_json else "",
            f"--in-img {ome_path}" if not use_json else "",
            f"--out-dir {cart_dir}",
            f"--img-id {image_id}",
            f"--colorize \"{image_color}\"",
            f"--upper-thres-quantile 0.95",
            f"--lower-thres-quantile 0.5",
            f"--transparent-below {args.transparent_below}" if args.transparent_below is not None else "",
            f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
        ])
        import_image_cmd = add_param_to_cmd(import_image_cmd, args, ["pmtiles", "gdaladdo"])
        import_image_cmd = add_param_to_cmd(import_image_cmd, args, ["restart","n_jobs"])

        run_command_w_preq(import_image_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)


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

def stage_run_ficture2(fic_dir, sge_assets, args, prereq):
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
    ])
    
    # Add optional params
    ficture2_cmd = add_param_to_cmd(ficture2_cmd, args, ["spatula", "ficture2"])
    ficture2_cmd = add_param_to_cmd(ficture2_cmd, args, ["restart","n_jobs","threads"])

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
def stage_run_cartload2(cart_dir, fic_dir, sge_dir, cell_assets, background_assets, args, prereq):

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

    # Build cartload2 command
    cartload_cmd = " ".join([
        "cartloader", "run_cartload2",
        f"--out-dir {cart_dir}",
        f"--fic-dir {fic_dir}" if args.run_ficture2 else "",
        f"--ext-fic-dir {args.ext_fic_dir}" if args.import_ext_ficture2 else "",
        f"--sge-dir {sge_dir}" if not (args.run_ficture2 or args.import_ext_ficture2) else "",
        f"--cell-assets {cell_assets}" if args.import_cells else "",
        f"--background-assets {' '.join(background_assets)}" if background_assets else "",
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
        f"--s3-dir s3://{args.s3_bucket}/{args.id}",
    ])

    # Inject optional params (preserves your original calls)
    aws_cmd = add_param_to_cmd(aws_cmd, args, ["aws"])
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
        print(" * Skip --upload-zenodo since all flags exist. Use --restart to force execution.", flush=True)
        print(zenodo_cmd, flush=True)
    else:
        run_command_w_preq(zenodo_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)
