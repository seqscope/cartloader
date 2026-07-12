import os, sys, argparse, subprocess, inspect
import yaml
from cartloader.utils.utils import create_custom_logger

logger = create_custom_logger(__name__)


def parse_arguments(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="AI-annotate the factors of a CartoStore output folder (or an uploaded S3 dataset), "
                    "writing per-factor alias TSVs and recording them in the catalog.")
    parser.add_argument('--s3-dir', type=str, help='S3 directory to annotate in place (mutually exclusive with --cartl-dir)')
    parser.add_argument('--cartl-dir', type=str, help='Local cartloader output directory to annotate')
    parser.add_argument('--tmp-dir', type=str, help='Temporary directory (S3 mode only; defaults under the current dir)')
    parser.add_argument('--tissue', type=str, required=True, help='Tissue name (passed to the annotation prompt)')
    parser.add_argument('--organism', type=str, required=True, help='Organism/species name (passed to the annotation prompt)')
    parser.add_argument('--api-type', type=str, default='claude', help='Generative-AI API type (default: claude)')
    parser.add_argument('--model', type=str, default='claude-opus-4-8', help='Model name (default: claude-opus-4-8)')
    parser.add_argument('--threads', type=int, default=1, help='Threads for annotation (default: 1)')
    # multi-sample
    parser.add_argument('--multi-sample', action='store_true', help='Annotate the shared factors in the multi-catalog, then propagate (reuse) the annotations into every per-sample sub-folder')
    parser.add_argument('--multi-catalog', type=str, default='multi-catalog.yaml', help='Multi-sample catalog file name under --cartl-dir (default: multi-catalog.yaml)')
    parser.add_argument('--reuse-results-from', type=str, help='Reuse existing alias TSVs from this directory instead of re-annotating')
    # naming / catalog keys
    parser.add_argument('--catalog', type=str, default='catalog.yaml', help='Per-sample catalog file name (default: catalog.yaml)')
    parser.add_argument('--alias-suffix', type=str, default='-alias-ai.tsv', help='Suffix for alias files (default: -alias-ai.tsv)')
    parser.add_argument('--backup-suffix', type=str, default='.bak', help='Suffix for the catalog backup (default: .bak)')
    parser.add_argument('--yaml-key-skip', type=str, nargs='+', default=['alias', 'alias_ai'], help='Factor keys that mark an existing alias (skip if present)')
    parser.add_argument('--yaml-key-store', type=str, default='alias_ai', help='Factor key under which the alias file is recorded (default: alias_ai)')
    # S3
    parser.add_argument('--profile', type=str, default='default', help='AWS profile for S3 upload (default: cartostore)')
    parser.add_argument('--skip-upload', action='store_true', help='S3 mode: do not upload alias files / updated catalog back to S3')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)


def _run(cmd):
    logger.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def iter_factors(catalog):
    """Yield (factor_id, factor_dict) pairs for either a multi-catalog (top-level
    `factors` map) or a per-sample catalog (`assets.factors` list)."""
    if isinstance(catalog.get("factors"), dict):
        return list(catalog["factors"].items())
    return [(f["id"], f) for f in catalog.get("assets", {}).get("factors", [])]


def annotate_catalog(cat_dir, cat_name, args, reuse_from=None):
    """Annotate every factor of one catalog in place. When ``reuse_from`` is given,
    alias files are copied from there instead of being (re)generated. Returns the
    loaded catalog dict (so the caller can read e.g. its `samples`)."""
    cat_path = os.path.join(cat_dir, cat_name)
    with open(cat_path) as f:
        catalog = yaml.safe_load(f)

    for factor_id, factor in iter_factors(catalog):
        if any(k in factor for k in args.yaml_key_skip):
            logger.info(f"Factor {factor_id} already annotated; skipping")
            continue
        de = factor.get("de")
        if not de:
            logger.info(f"Factor {factor_id} has no 'de' entry; skipping")
            continue

        # Normalize to hyphens so a shared (multi) alias matches the per-sample
        # factor id (t12_f48 -> t12-f48) and can be reused across sub-folders.
        alias_tsv = f"{factor_id.replace('_', '-')}{args.alias_suffix}"
        alias_path = os.path.join(cat_dir, alias_tsv)

        if reuse_from is not None:
            src = os.path.join(reuse_from, alias_tsv)
            if not os.path.exists(src):
                raise FileNotFoundError(f"Cannot reuse annotation for factor '{factor_id}': {src} not found")
            logger.info(f"Reusing annotation for {factor_id} from {src}")
            _run(f"cp -f {src} {alias_path}")
        else:
            logger.info(f"Annotating factor {factor_id} from {de}")
            _run(f"cartloader annotate_bulk_de_with_ai --de {os.path.join(cat_dir, de)} "
                 f"--out {alias_path} --api-type {args.api_type} --model-name {args.model} "
                 f"--template-args tissue=\"{args.tissue}\" organism=\"{args.organism}\" --threads {args.threads}")

        factor[args.yaml_key_store] = alias_tsv

    # Back up and rewrite the catalog with the recorded aliases.
    _run(f"cp -f {cat_path} {cat_path}{args.backup_suffix}")
    with open(cat_path, "w") as f:
        yaml.dump(catalog, f, Dumper=yaml.SafeDumper, default_flow_style=False, sort_keys=False)
    return catalog


def anno_cartload_folder(_args):
    args = parse_arguments(_args)

    if args.cartl_dir and args.s3_dir:
        raise ValueError("Specify either --cartl-dir or --s3-dir, not both.")

    if args.cartl_dir:
        cartl = args.cartl_dir.rstrip('/')
        if args.multi_sample:
            # 1) annotate the shared factors recorded in the multi-catalog
            logger.info(f"Annotating shared factors in {os.path.join(cartl, args.multi_catalog)}")
            mc = annotate_catalog(cartl, args.multi_catalog, args, reuse_from=args.reuse_results_from)
            # 2) propagate (reuse) the shared annotations into every sample sub-folder
            for sid, rel in (mc.get("samples") or {}).items():
                sub_dir = os.path.join(cartl, os.path.dirname(rel))
                logger.info(f"Propagating shared annotations to sample '{sid}' ({sub_dir})")
                annotate_catalog(sub_dir, os.path.basename(rel), args, reuse_from=cartl)
        else:
            annotate_catalog(cartl, args.catalog, args, reuse_from=args.reuse_results_from)
        logger.info("Annotation complete (local; no upload).")
        return

    if not args.s3_dir:
        raise ValueError("Specify either --cartl-dir or --s3-dir.")

    # --- S3 mode: download the dataset, annotate locally, upload results back ---
    import boto3
    s3_dir = args.s3_dir.rstrip('/')
    if s3_dir.startswith('s3://'):
        s3_dir = s3_dir[5:]
    bucket = s3_dir.split('/')[0]
    prefix = '/'.join(s3_dir.split('/')[1:])
    tmp = args.tmp_dir or s3_dir.split('/')[-1]
    os.makedirs(tmp, exist_ok=True)

    s3 = boto3.client('s3')
    logger.info(f"Downloading catalog s3://{bucket}/{prefix}/{args.catalog} -> {tmp}")
    s3.download_file(Bucket=bucket, Key=f"{prefix}/{args.catalog}", Filename=os.path.join(tmp, args.catalog))
    with open(os.path.join(tmp, args.catalog)) as f:
        catalog = yaml.safe_load(f)
    for factor_id, factor in iter_factors(catalog):
        de = factor.get("de")
        if de and not any(k in factor for k in args.yaml_key_skip):
            logger.info(f"Downloading {de}")
            s3.download_file(Bucket=bucket, Key=f"{prefix}/{de}", Filename=os.path.join(tmp, de))

    annotate_catalog(tmp, args.catalog, args, reuse_from=args.reuse_results_from)

    if args.skip_upload:
        logger.info("Skipping upload (--skip-upload).")
        return
    up = boto3.Session(profile_name=args.profile).client('s3')
    with open(os.path.join(tmp, args.catalog)) as f:
        updated = yaml.safe_load(f)
    for _, factor in iter_factors(updated):
        alias = factor.get(args.yaml_key_store)
        if alias:
            logger.info(f"Uploading {alias}")
            up.upload_file(Bucket=bucket, Key=f"{prefix}/{alias}", Filename=os.path.join(tmp, alias))
    logger.info(f"Uploading updated {args.catalog}")
    up.upload_file(Bucket=bucket, Key=f"{prefix}/{args.catalog}", Filename=os.path.join(tmp, args.catalog))


if __name__ == "__main__":
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], script_name)
    func(sys.argv[1:])
