import sys, os, argparse, logging, inspect, json, glob
import yaml

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, add_param_to_cmd, execute_makefile
from cartloader.utils.cartload_helper import render_umap_cmd, umap_tippecanoe_cmd

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# run_cartload2 options forwarded to each per-sample invocation
aux_args = {
    "params": [
        "in_fic_params", "out_fic_assets", "out_catalog",
        "rename_x", "rename_y", "colname_feature", "colname_count",
        "out_molecules_id", "max_join_dist_um", "join_tile_size", "bin_count",
        "preserve_point_density_thres",
        "umap_colname_factor", "umap_colname_x", "umap_colname_y", "umap_min_zoom", "umap_max_zoom",
        "sge_scale", "use_pmpoint", "tile_format_pmpoint",
        "skip_umap", "skip_raster",
        "tmp_dir", "keep_intermediate_files",
        "transparent_below", "transparent_above",
    ],
    "env": ["gzip", "pmtiles", "gdal_translate", "gdaladdo", "tippecanoe", "spatula", "pmpoint", "ficture2"],
    "run": ["restart", "threads", "log", "log_suffix"],
}


def parse_arguments(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Package a multi-sample FICTURE run into per-sample PMTiles (one run_cartload2 per sample, "
                    "in parallel) plus a shared multi-catalog YAML.")

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate the Makefile but do not execute it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and re-run all steps')
    run_params.add_argument('--makefn', type=str, default="run_cartload2_multi.mk", help='Name of the generated Makefile')
    run_params.add_argument('-j', '--n-jobs', type=int, default=1, help='Number of samples to package in parallel (default: 1)')
    run_params.add_argument('--threads', type=int, default=None, help='Threads per job (forwarded to run_cartload2)')
    run_params.add_argument('--log', action='store_true', default=False, help='Write logs to a file under the output directory')
    run_params.add_argument('--log-suffix', type=str, default=None, help='Suffix for the log filename')

    inout_params = parser.add_argument_group("Input/Output Parameters")
    inout_params.add_argument('--out-dir', type=str, required=True, help='Output directory. Holds one self-contained sub-directory per sample plus the multi-catalog YAML')
    inout_params.add_argument('--fic-dir', type=str, required=True, help='Multi-sample FICTURE directory produced by run_ficture2_multi; must contain the multi manifest (see --in-multi-params)')
    inout_params.add_argument('--id', type=str, default=None, help='Identifier for the multi-sample run; per-sample directories are named <id>-<sample_id> (default: basename of --out-dir)')
    inout_params.add_argument('--in-multi-params', type=str, default="ficture.multi.params.json", help='File name of the shared multi manifest under --fic-dir (default: ficture.multi.params.json)')
    inout_params.add_argument('--multi-catalog', type=str, default="multi-catalog.yaml", help='File name of the output multi-sample catalog under --out-dir (default: multi-catalog.yaml)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Forwarded to each per-sample run_cartload2; defaults work for most cases")
    aux_params.add_argument('--in-fic-params', type=str, help='FICTURE params JSON under each sample dir (default: ficture.params.json)')
    aux_params.add_argument('--out-fic-assets', type=str, help='Output FICTURE assets JSON name')
    aux_params.add_argument('--out-catalog', type=str, default="catalog.yaml", help='Per-sample catalog YAML name (default: catalog.yaml)')
    aux_params.add_argument('--colname-feature', type=str, help='Column name for feature/gene')
    aux_params.add_argument('--colname-count', type=str, help='Column name for molecule counts')
    aux_params.add_argument('--rename-x', type=str, help='Column rename mapping for X (old:new)')
    aux_params.add_argument('--rename-y', type=str, help='Column rename mapping for Y (old:new)')
    aux_params.add_argument('--out-molecules-id', type=str, help='Base name for output molecules PMTiles')
    aux_params.add_argument('--max-join-dist-um', type=float, help='Max distance (µm) to associate molecules with pixels')
    aux_params.add_argument('--join-tile-size', type=float, help='Tile size (µm) when joining molecules with pixels')
    aux_params.add_argument('--bin-count', type=int, help='Number of bins when splitting input molecules')
    aux_params.add_argument('--preserve-point-density-thres', type=int, help='Tippecanoe point-density preservation threshold')
    aux_params.add_argument('--sge-scale', type=int, help='Scale factor from input coordinates to output pixels')
    aux_params.add_argument('--use-pmpoint', action='store_true', default=False, help='Use pmpoint/MLT instead of tippecanoe for point PMTiles')
    aux_params.add_argument('--tile-format-pmpoint', choices=['MLT', 'MVT'], help='Tile format when --use-pmpoint is enabled')
    aux_params.add_argument('--umap-colname-factor', type=str, help='UMAP dominant-factor column name')
    aux_params.add_argument('--umap-colname-x', type=str, help='UMAP X column name')
    aux_params.add_argument('--umap-colname-y', type=str, help='UMAP Y column name')
    aux_params.add_argument('--umap-min-zoom', type=int, help='Minimum zoom for UMAP PMTiles')
    aux_params.add_argument('--umap-max-zoom', type=int, help='Maximum zoom for UMAP PMTiles')
    aux_params.add_argument('--skip-umap', action='store_true', default=False, help='Skip UMAP PMTiles')
    aux_params.add_argument('--skip-raster', action='store_true', default=False, help='Skip raster image generation')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate files')
    aux_params.add_argument('--transparent-below', type=int, help='Transparent below this pixel value')
    aux_params.add_argument('--transparent-above', type=int, help='Transparent above this pixel value')

    env_params = parser.add_argument_group("Env Parameters")
    env_params.add_argument('--gzip', type=str, help='Path to gzip-compatible binary')
    env_params.add_argument('--pmtiles', type=str, help='Path to pmtiles binary')
    env_params.add_argument('--gdal_translate', type=str, help='Path to gdal_translate binary')
    env_params.add_argument('--gdaladdo', type=str, help='Path to gdaladdo binary')
    env_params.add_argument('--tippecanoe', type=str, help='Path to tippecanoe binary')
    env_params.add_argument('--spatula', type=str, help='Path to spatula binary')
    env_params.add_argument('--pmpoint', type=str, help='Path to pmpoint binary')
    env_params.add_argument('--ficture2', type=str, help='Path to punkst (ficture2) repository')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)


def build_multi_catalog(mm, args, manifest, cells_manifests, samples, multi_id):
    """Assemble the multi-sample catalog and add targets that materialize the shared
    factor files at the out_dir root, sourced from the FICTURE output (via the
    manifests) with run_cartload2-consistent naming. Returns (catalog, shared_flags)."""
    out_catalog = args.out_catalog
    fic, root = args.fic_dir, args.out_dir
    gzip = args.gzip or "gzip"
    factors = {}
    shared_flags = []

    def cp(src_rel, dst):
        return f'cp -f "{os.path.join(fic, src_rel)}" "{os.path.join(root, dst)}"'

    def rgb(src_rel, dst):   # normalize the FICTURE cmap into an rgb.tsv (as run_cartload2 does)
        s = os.path.join(fic, src_rel); d = os.path.join(root, dst)
        return (f'python3 -c "from cartloader.utils.cartload_helper import copy_rgb_tsv; '
                f"copy_rgb_tsv('{s}', '{d}', restart=True)\"")

    def umap_cmds(src_tsv_rel, src_png_rel, oid):
        tip = args.tippecanoe or f"{repo_dir}/submodules/tippecanoe/tippecanoe"
        ndjson = os.path.join(args.tmp_dir, f"{oid}-umap.ndjson")
        out_tsv = os.path.join(root, f"{oid}-umap.tsv.gz")
        out_pmt = os.path.join(root, f"{oid}-umap.pmtiles")
        return [
            cp(src_tsv_rel, f"{oid}-umap.tsv.gz"),
            cp(src_png_rel, f"{oid}-umap.png"),
            render_umap_cmd(out_tsv, ndjson,
                            args.umap_colname_factor or "topK",
                            args.umap_colname_x or "UMAP1",
                            args.umap_colname_y or "UMAP2"),
            umap_tippecanoe_cmd(out_pmt, ndjson, tip, args.tmp_dir,
                                threads=args.threads or 4,
                                min_zoom=args.umap_min_zoom if args.umap_min_zoom is not None else 0,
                                max_zoom=args.umap_max_zoom if args.umap_max_zoom is not None else 18,
                                preserve_thres=args.preserve_point_density_thres or 1024),
            f"rm -f {ndjson}",
        ]

    # --- model-derived factors (from the joint LDA models) ---
    for tp in manifest.get("shared", {}).get("train_params", []):
        oid = tp["model_id"].replace("_", "-")
        cmds = cmd_separator([], f"Materializing shared factor {tp['model_id']}")
        # `post` keeps its provenance suffix: -model.tsv for an LDA model.
        cmds += [cp(tp["model_path"], f"{oid}-model.tsv"),
                 rgb(tp["cmap"], f"{oid}-rgb.tsv"),
                 cp(tp["de_path"], f"{oid}-de.tsv"),
                 cp(tp["info_path"], f"{oid}-info.tsv")]
        entry = {"post": f"{oid}-model.tsv", "rgb": f"{oid}-rgb.tsv",
                 "de": f"{oid}-de.tsv", "info": f"{oid}-info.tsv"}
        u = tp.get("umap")
        if u:
            cmds += umap_cmds(u["tsv"], u["png"], oid)
            entry["umap"] = {"pmtiles": f"{oid}-umap.pmtiles", "png": f"{oid}-umap.png",
                             "tsv": f"{oid}-umap.tsv.gz"}
        flag = os.path.join(root, f"{oid}.shared.done")
        cmds.append(f"touch {flag}")
        mm.add_target(flag, [os.path.join(fic, tp["model_path"])], cmds)
        shared_flags.append(flag)
        factors[tp["model_id"]] = entry

    # --- cell-derived factors (same keys) ---
    for cm in cells_manifests:
        prefix = cm.get("out_prefix"); oid = prefix.replace("_", "-")
        csh = cm.get("shared", {})
        cmds = cmd_separator([], f"Materializing shared cell factor {prefix}")
        entry, prereqs = {}, []
        if "shared_cluster_pseudobulk" in csh:
            src = csh["shared_cluster_pseudobulk"]
            # `post` keeps its provenance suffix: -pseudobulk.tsv.gz for cell clusters.
            cmds.append(f'{gzip} -c "{os.path.join(fic, src)}" > "{os.path.join(root, oid + "-pseudobulk.tsv.gz")}"')
            entry["post"] = f"{oid}-pseudobulk.tsv.gz"; prereqs.append(os.path.join(fic, src))
        if "shared_cmap" in csh:
            cmds.append(rgb(csh["shared_cmap"], f"{oid}-rgb.tsv")); entry["rgb"] = f"{oid}-rgb.tsv"
        if "shared_cluster_de" in csh:
            cmds.append(cp(csh["shared_cluster_de"], f"{oid}-de.tsv")); entry["de"] = f"{oid}-de.tsv"
        if "shared_cluster_info" in csh:
            cmds.append(cp(csh["shared_cluster_info"], f"{oid}-info.tsv")); entry["info"] = f"{oid}-info.tsv"
        um = csh.get("manifolds", {}).get("umap")
        if um:
            cmds += umap_cmds(um["tsv"], um["png"], oid)
            entry["umap"] = {"pmtiles": f"{oid}-umap.pmtiles", "png": f"{oid}-umap.png",
                             "tsv": f"{oid}-umap.tsv.gz"}
        if "shared_cluster_model_heatmap_pdf" in csh:
            cmds.append(cp(csh["shared_cluster_model_heatmap_pdf"], f"{oid}-heatmap.pdf"))
            cmds.append(cp(csh["shared_cluster_model_heatmap_tsv"], f"{oid}-heatmap.tsv"))
            entry["heatmap"] = {"pdf": f"{oid}-heatmap.pdf", "tsv": f"{oid}-heatmap.tsv"}
        flag = os.path.join(root, f"{oid}.shared.done")
        cmds.append(f"touch {flag}")
        mm.add_target(flag, prereqs, cmds)
        shared_flags.append(flag)
        factors[prefix] = entry

    catalog = {
        "id": multi_id,
        "analysis_type": "multi-sample",
        "n_samples": len(samples),
        "samples": {sid: os.path.join(f"{multi_id}-{sid}", out_catalog) for sid in samples},
        "factors": factors,
    }
    return catalog, shared_flags


def run_cartload2_multi(_args):
    args = parse_arguments(_args)
    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(args.out_dir, "tmp")

    os.makedirs(args.out_dir, exist_ok=True)

    # Read the shared multi manifest written by run_ficture2_multi
    manifest_path = os.path.join(args.fic_dir, args.in_multi_params)
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(
            f"Multi manifest not found: {manifest_path}. Run run_ficture2_multi first (it writes {args.in_multi_params}).")
    with open(manifest_path) as f:
        manifest = json.load(f)
    samples = list(manifest.get("samples", {}).keys())
    if not samples:
        raise ValueError(f"No samples found in {manifest_path}")

    # Any cell-analysis multi manifests (ficture.multi.<prefix>.params.json)
    cells_manifests = []
    for p in sorted(glob.glob(os.path.join(args.fic_dir, "ficture.multi.*.params.json"))):
        if os.path.basename(p) == args.in_multi_params:
            continue
        with open(p) as f:
            cells_manifests.append(json.load(f))

    multi_id = args.id if args.id is not None else os.path.basename(os.path.normpath(args.out_dir))

    mm = minimake()
    sample_catalogs = {}
    for sid in samples:
        sample_rel = manifest["samples"][sid]                     # e.g. samples/<sid>/ficture.params.json
        sample_fic_dir = os.path.join(args.fic_dir, os.path.dirname(sample_rel))
        out_id = f"{multi_id}-{sid}"
        cart_dir = os.path.join(args.out_dir, out_id)
        catalog_yaml = os.path.join(cart_dir, args.out_catalog)
        sample_catalogs[sid] = catalog_yaml

        # Auto-detect per-sample cell analyses (ficture.<prefix>.params.json, excluding the base)
        cell_params = [p for p in sorted(glob.glob(os.path.join(sample_fic_dir, "ficture.*.params.json")))
                       if os.path.basename(p) != "ficture.params.json"]

        prereqs = [manifest_path, os.path.join(args.fic_dir, sample_rel)]

        cmds = cmd_separator([], f"Packaging sample {sid} (run_cartload2)")
        cmd = " ".join([
            "cartloader", "run_cartload2",
            f"--out-dir {cart_dir}",
            f"--fic-dir {sample_fic_dir}",
            f"--id {out_id}",
            ("--in-cell-params " + " ".join(cell_params)) if cell_params else "",
            "--makefn run_cartload2.mk",
        ])
        cmd = add_param_to_cmd(cmd, args, aux_args["params"])
        cmd = add_param_to_cmd(cmd, args, aux_args["env"])
        cmd = add_param_to_cmd(cmd, args, aux_args["run"])
        cmds.append(cmd)
        mm.add_target(catalog_yaml, prereqs, cmds)

    if len(mm.targets) == 0:
        raise ValueError("No tasks were generated. Check inputs and parameters.")

    os.makedirs(args.tmp_dir, exist_ok=True)

    # Build the multi-catalog and add targets that materialize the shared factor
    # files at the out_dir root — sourced from the FICTURE output (via the
    # manifests) and processed with run_cartload2-consistent naming.
    multi_catalog, _shared_flags = build_multi_catalog(mm, args, manifest, cells_manifests, samples, multi_id)
    with open(os.path.join(args.out_dir, args.multi_catalog), "w") as f:
        yaml.safe_dump(multi_catalog, f, sort_keys=False)

    make_f = os.path.join(args.out_dir, args.makefn)
    mm.write_makefile(make_f)
    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)


if __name__ == "__main__":
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], script_name)
    func(sys.argv[1:])
