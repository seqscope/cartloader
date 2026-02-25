import sys, os, argparse, logging, subprocess, inspect
import pandas as pd
from pathlib import Path

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, add_param_to_cmd, execute_makefile, valid_and_touch_cmd


aux_args = {
    "params": [
        "in_fic_params", "out_fic_assets", "out_catalog",
        "rename_x","rename_y",
        "colname_feature", "colname_count",
        "out_molecules_id", "max_join_dist_um", "join_tile_size", "bin_count", "max_tile_bytes", "max_feature_counts", "preserve_point_density_thres",
        "umap_colname_factor", "umap_colname_x", "umap_colname_y", "umap_min_zoom", "umap_max_zoom",
        "skip_raster",
        "tmp_dir", "keep_intermediate_files",
        "transparent_below", "transparent_above"
    ],
    "env":[
        "gzip", "pmtiles", "gdaladdo", "tippecanoe", "spatula"
    ],
    "run":[
        "restart", "n_jobs", "threads", "log", "log_suffix"
    ]
}

def parse_arguments(_args):
    """
    Package SGE and optional FICTURE outputs into PMTiles and a catalog for CartoScope (sge → optional ficture2/joins → pmtiles → catalog)
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Package SGE and optional FICTURE outputs into PMTiles and a catalog for CartoScope (sge → optional ficture2/joins → pmtiles → catalog)"
    )

    run_params = parser.add_argument_group("Run Options", "Execution controls for generating and running the Makefile")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate the Makefile but do not execute it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and re-run all steps')
    run_params.add_argument('--makefn', type=str, default="run_cartload2_multi.mk", help='Name of the generated Makefile (default: run_cartload2.mk)')
    run_params.add_argument('--n-jobs', type=int, default=None, help='Number of parallel jobs to run (default: 1)') # pass down to sub-makefiles
    run_params.add_argument('--threads', type=int, default=None, help='Maximum number of threads per job for tippecanoe') # pass down to sub-makefiles
    run_params.add_argument('--log', action='store_true', default=False, help='Write logs to a file under the output directory') # pass down to sub-makefiles
    run_params.add_argument('--log-suffix', type=str, default=None, help='Suffix for the log filename; final path is <out_dir>_cartload<suffix> (default: .log)') # pass down to sub-makefiles

    inout_params = parser.add_argument_group("Input/Output Parameters")
    inout_params.add_argument('--out-dir', type=str, required=True, help='Output directory (PMTiles, assets JSON, and catalog YAML)')
    inout_params.add_argument('--fic-dir', type=str, required=True, help='Path tp FICTURE results directory produced by "cartloader run_ficture"; must include FICTURE results and a parameter JSON (see --in-fic-params)')
    inout_params.add_argument('--in-list', type=str, required=True, help='Path to input TSV without header. One sample per row with 2 required columns [sample id, path to transcript] and 2 optional columns [title, description]. If id is omitted it defaults to sample (use lowercase and replace _ by -); quote title/description if they contain spaces.')

    ## Show the default value in help text but not set it here. All defaults will be applied in run_cartload2()
    aux_params = parser.add_argument_group("Auxiliary Parameters", "Advanced settings; defaults work for most cases") 
    aux_params.add_argument('--in-fic-params', type=str, help='File name of FICTURE params JSON/YAML under --fic-dir, providing FICTURE paramaters (default: ficture.params.json)')
    aux_params.add_argument('--out-fic-assets', type=str, help='File name of output JSON/YAML for FICTURE asset metadata under --out-dir (default: ficture_assets.json)')
    aux_params.add_argument('--out-catalog', type=str,  help='File name of output catalog YAML under --out-dir (default: catalog.yaml)')
    aux_params.add_argument('--rename-x', type=str, help='Column rename mapping for X axis in tippecanoe, format old:new (default: x:lon)')  
    aux_params.add_argument('--rename-y', type=str, help='Column rename mapping for Y axis in tippecanoe, format old:new (default: y:lat)')  
    aux_params.add_argument('--colname-feature', type=str, help='Column name for feature/gene (default: gene)')
    aux_params.add_argument('--colname-count', type=str, help='Column name for molecule counts (default: count)')
    aux_params.add_argument('--out-molecules-id', type=str, help='Base name for output molecules PMTiles files (no directory)')
    aux_params.add_argument('--max-join-dist-um', type=float, help='Max distance (in µm) to associate molecules with decoded pixels (default: 0.1)')
    aux_params.add_argument('--join-tile-size', type=float, help='Tile size (in µm) when joining molecules with decoded pixels (default: 500)')
    aux_params.add_argument('--bin-count', type=int, help='Number of bins when splitting input molecules (default: 50)')
    aux_params.add_argument('--max-tile-bytes', type=int, help='Maximum tile size in bytes for tippecanoe/PMTiles (default: 5000000)')
    aux_params.add_argument('--max-feature-counts', type=int, help='Maximum features per tile for tippecanoe/PMTiles (default: 500000)')
    aux_params.add_argument('--preserve-point-density-thres', type=int, help='Tippecanoe point-density preservation threshold (default: 1024)')
    aux_params.add_argument('--umap-colname-factor', type=str, help='Column name encoding the dominant factor assignment in a UMAP TSV (default: topK)')
    aux_params.add_argument('--umap-colname-x', type=str, help='Column name for the UMAP X coordinate (default: UMAP1)')
    aux_params.add_argument('--umap-colname-y', type=str, help='Column name for the UMAP Y coordinate (default: UMAP2)')
    aux_params.add_argument('--umap-min-zoom', type=int, help='Minimum zoom for generated UMAP PMTiles (default: 0)')
    aux_params.add_argument('--umap-max-zoom', type=int, help='Maximum zoom for generated UMAP PMTiles (default: 18)')
    aux_params.add_argument('--skip-raster', action='store_true', help='Skip raster image generation (no GDAL/go-pmtiles required)')
    aux_params.add_argument('--tmp-dir', type=str, help='Temporary directory (default: <out_dir>/tmp)')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', help='Keep intermediate files instead of cleaning up')
    aux_params.add_argument('--transparent-below', type=int, help='Set pixels below this value to transparent for dark background (range: 0~255)')
    aux_params.add_argument('--transparent-above', type=int, help='Set pixels above this value to transparent for light background (range: 0~255)')

    ## Show the default value in help text but not set it here. All defaults will be applied in run_cartload2()
    env_params = parser.add_argument_group("Env Parameters", "Tool paths (override defaults if needed)")
    env_params.add_argument('--gzip', type=str, help='Path to gzip-compatible binary (tip: use "pigz -p4" for speed)')
    env_params.add_argument('--pmtiles', type=str, help='Path to pmtiles binary from go-pmtiles')
    env_params.add_argument('--gdal_translate', type=str, help='Path to gdal_translate binary')
    env_params.add_argument('--gdaladdo', type=str, help='Path to gdaladdo binary')
    env_params.add_argument('--tippecanoe', type=str, help='Path to tippecanoe binary') # default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", 
    env_params.add_argument('--spatula', type=str, help='Path to spatula binary') # default=f"{repo_dir}/submodules/spatula/bin/spatula",

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    
    args=parser.parse_args(_args)

    # dir
    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(args.out_dir, "tmp")
    
    return args

def run_cartload2_multi(_args):
    """
    Build resources for CartoScope
    """

    # parse argument
    args=parse_arguments(_args)

    # start mm
    mm = minimake()

    # files/dirs
    os.makedirs(args.out_dir, exist_ok=True)
    ficture2_flag = os.path.join(args.fic_dir, "multi_json_each.done")

    # read the in_list
    df=pd.read_csv(args.in_list, sep="\t", dtype=str, header=None)
    
    if df.shape[1] < 2:
        raise ValueError(f"Input list {args.in_list} must have at least 2 columns: sample id and path to transcript")
    
    df.columns = ["id", "transcript_path"] + (["title"] if df.shape[1] >= 3 else []) + (["desc"] if df.shape[1] >= 4 else [])
    
    # process by row
    for idx, row in df.iterrows():
        id = row["id"]
        transcript_path = row["transcript_path"]
        title = row["title"] if "title" in df.columns else None
        desc = row["desc"] if "desc" in df.columns else None
        out_id = id.replace("_", "-").lower()   
        
        fic_dir = os.path.join(args.fic_dir, "samples", id)
        fic_param_assets = os.path.join(fic_dir, "ficture.params.json")

        cart_dir = os.path.join(args.out_dir, id)
        os.makedirs(cart_dir, exist_ok=True)
        catalog_yaml = os.path.join(cart_dir, args.out_catalog if args.out_catalog is not None else "catalog.yaml")

        prerequisites = [ficture2_flag, fic_param_assets]

        cmds = cmd_separator([], f"Running run_cartload2 for sample {id}")
        cmd = " ".join([
            "cartloader", "run_cartload2",
            f"--out-dir {cart_dir}",
            f"--fic-dir {fic_dir}",
            f"--id {out_id}",
            f"--title '{title}'" if title is not None else "",
            f"--desc '{desc}'" if desc is not None else "",
            f"--gdal_translate '{args.gdal_translate}'" if args.gdal_translate is not None else "",
            f"--makefn run_cartload2.mk",
        ])
        cmd = add_param_to_cmd(cmd, args, aux_args["params"])
        cmd = add_param_to_cmd(cmd, args, aux_args["env"])
        cmd = add_param_to_cmd(cmd, args, aux_args["run"])
        cmds.append(cmd)
        mm.add_target(catalog_yaml, prerequisites, cmds)

    if len(mm.targets) == 0:
        raise ValueError("No tasks were generated. Check inputs and parameters.")

    ## write makefile
    make_f = os.path.join(args.out_dir, args.makefn)
    mm.write_makefile(make_f)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=1)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
