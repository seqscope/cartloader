import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, re
import pandas as pd
import numpy as np

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader run_fig2pmtiles", description="Convert a figure to pmtiles")

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--main', action='store_true', default=False, help='Run main commands (geotif2mbtiles, mbtiles2pmtiles, upload-aws, update-yaml)')
    cmd_params.add_argument('--geotif2mbtiles', action='store_true', default=False, help='Convert a geotiff file to mbtiles')
    cmd_params.add_argument('--mbtiles2pmtiles', action='store_true', default=False, help='Convert mbtiles to pmtiles')
    cmd_params.add_argument('--upload-aws', action='store_true', default=False, help='Upload the new pmtiles to the AWS S3 bucket')
    cmd_params.add_argument('--update-yaml', action='store_true', default=False, help='Update the catalog.yaml file with the new pmtiles and upload to AWS')
    cmd_params.add_argument('--georeference', action='store_true', default=False, help='Plus function. Create a geotiff file from PNG or TIF file. If enabled, the user must provide georeferenced bounds using --in-tsv or --in-bounds.')
    cmd_params.add_argument('--flip-vertical', action='store_true', default=False, help='Plus function. Flip the image vertically')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Define the input file according to the user's needs.")
    inout_params.add_argument('--in-fig', type=str, help='The input figure file (PNG or TIF) to be converted to pmTiles')
    inout_params.add_argument('--in-tsv', type=str, default=None, help='If --georeference is required, use the *.pixel.sorted.tsv.gz from run_ficture to provide georeferenced bounds.')
    inout_params.add_argument('--in-bounds', type=str, default=None, help='If --georeference is required, provide the bounds in the format of "<ulx>,<uly>,<lrx>,<lry>", which represents upper-left X, upper-left Y, lower-right X, lower-right Y.')
    inout_params.add_argument('--out-prefix', required=True, type=str, help='The output prefix. TNew directory will be created if needed')
    inout_params.add_argument('--basemap-key', type=str, default=None, help='The key to use for updating the basemap in the catalog.yaml file')

    key_params = parser.add_argument_group("Key parameters")
    key_params.add_argument('--srs', type=str, default='EPSG:3857', help='For the georeference and geo2tiff steps, define the spatial reference system (default: EPSG:3857)')
    key_params.add_argument('--resample', type=str, default='cubic', help='For the geo2tiff step, define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.')
    key_params.add_argument('--blocksize', type=int, default='512', help='For the geo2tiff step, define the blocksize (default: 512)')
    key_params.add_argument('--aws-bucket', type=str, default=None, help='For the update-aws step, define the path to the AWS S3 bucket')
    key_params.add_argument('--yaml', type=str, default=None, help='For the update-aws step, define the yaml file to update (default: catalog.yaml in the output directory specified by --out-prefix)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles')
    aux_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    aux_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binary')
    aux_params.add_argument('--keep-files', action='store_true', default=False, help='Keep intermediate files')

    
    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default=None, help='The file name of the Makefile to generate (default: {out-prefix}.mk)')
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it (default: False)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def run_fig2pmtiles(_args):

    args=parse_arguments(_args)

    if args.main:
        args.geotif2mbtiles = True
        args.mbtiles2pmtiles = True
        args.upload_aws = True
        args.update_yaml = True

    if args.update_yaml and args.basemap_key is None:
        raise ValueError("Please provide the basemap key using --basemap-key for updating the catalog.yaml file")

    scheck_app(args.pmtiles)
    scheck_app(args.gdal_translate)
    scheck_app(args.gdaladdo)

    # files
    out_dir = os.path.dirname(args.out_prefix)
    os.makedirs(out_dir, exist_ok=True)
    # - geotiff for 
    geotif_f = args.in_fig 
    # - mbtiles
    mbtile_f=f"{args.out_prefix}.pmtiles.mbtiles"
    mbtile_f_ann=f"{args.out_prefix}.pmtiles.{args.resample}.mbtiles"
    # - pmtiles
    pmtiles_f=f"{args.out_prefix}.pmtiles"
    catalog_f = f"{out_dir}/catalog.yaml" if args.yaml is None else args.yaml

    # start mm
    mm = minimake()

    cmds_rm = cmd_separator([], f"Removing intermediate files")

    # 0. Create a geotiff file from PNG or TIF file
    if args.georeference:
        georef_f =f"{args.out_prefix}.pmtiles.tif"
        assert args.in_fig is not None and os.path.exists(args.in_fig), "Please provide a valid input figure file using --in-figure"
        if args.in_bounds is not None:
            ulx,uly,lrx,lry = args.in_bounds.split(",")
        elif args.in_tsv is not None:
            with gzip.open(args.in_tsv, 'rt') as f:
                for i in range(3):
                    line = f.readline()
            ann2val = {x.split("=")[0]:x.split("=")[1] for x in line.strip().replace("##", "").split(";")}
            ulx = float(ann2val["OFFSET_X"]) # Upper left X-coordinate
            uly = float(ann2val["OFFSET_Y"]) # Upper left Y-coordinate
            lrx = float(ann2val["SIZE_X"])+1+ulx # Lower Right X-coordinate
            lry = float(ann2val["SIZE_Y"])+1+uly # Lower Right Y-coordinate 
        else:
            raise ValueError("Please provide either --in-bounds or --in-tsv to georeference the figure")
        
        cmds = cmd_separator([], f"Geo-referencing {args.in_fig} to {georef_f}")
        cmds.append(f"{args.gdal_translate} -of GTiff -a_srs {args.srs} -a_ullr {ulx} {uly} {lrx} {lry} {args.in_fig} {georef_f}")
        mm.add_target(georef_f, [args.in_fig], cmds)
        # update geotif_f
        geotif_f = georef_f

        if args.mbtiles2pmtiles:
            cmds_rm.append(f"rm -f {georef_f}")

    # 1. Flip the image when required 
    if args.flip_vertical:
        vflip_f = f"{args.out_prefix}.pmtiles.vflip.tif"
        cmds = cmd_separator([], f"Flipping the geotif vertically: {geotif_f}")        
        cmd = f'''
        INFO=$(gdalinfo "{geotif_f}" 2>&1)
        if [[ $INFO =~ Size\\ is\\ ([0-9]+),\\ ([0-9]+) ]]; then
            WIDTH=${{BASH_REMATCH[1]}}
            HEIGHT=${{BASH_REMATCH[2]}}
            echo "width: ${{WIDTH}}, height: ${{HEIGHT}}"
        else
            echo "Failed to extract image dimensions."
        fi
        '''
        cmds.append(cmd)
        cmd = " ".join([
                "gdalwarp",
                f'"{geotif_f}"',  # Add quotes around file names to handle spaces
                f'"{vflip_f}"',
                "-b", "1",
                "-b", "2",
                "-b", "3",
                "-ct", "\"+proj=pipeline +step +proj=axisswap +order=1,-2\"",
                "-overwrite",
                "-ts", "$WIDTH", "$HEIGHT"  # Use the WIDTH and HEIGHT from the previous command
            ])
        cmds.append(cmd)
        # update geotif_f
        geotif_f = vflip_f

        if args.mbtiles2pmtiles:
            cmds_rm.append(f"rm -f {geotif_f}")

    # 2. Convert a geotiff file to mbtiles
    if args.geotif2mbtiles:
        cmds = cmd_separator([], f"Converting from geotif to mbtiles: {geotif_f}")
        cmd = " ".join([
            args.gdal_translate, 
            "-b", "1", 
            "-b", "2",
            "-b", "3",
            "-strict",
            "-co", "\"ZOOM_LEVEL_STRATEGY=UPPER\"",
            "-co", f"\"RESAMPLING={args.resample}\"",
            "-co", f"\"BLOCKSIZE={args.blocksize}\"",
            "-ot", "Byte",
            "-scale", "0", "255",
            "-of", "mbtiles",
            "-a_srs", args.srs,
            geotif_f, 
            mbtile_f
        ])
        cmds.append(cmd)
        #cmds.append(f"gdal_translate -b 1 -b 2 -b 3 -strict -co \"ZOOM_LEVEL_STRATEGY=UPPER\" -co \"RESAMPLING={args.resample}\" -co \"BLOCKSIZE={args.blocksize}\" -ot Byte -scale 0 255 -of mbtiles -a_srs {args.srs} {geotif_f} {mbtile_f}")
        mm.add_target(mbtile_f, [geotif_f], cmds)
        if args.mbtiles2pmtiles:
            cmds_rm.append(f"rm -f {mbtile_f}")
    
    # 3. Convert mbtiles to pmtiles
    if args.mbtiles2pmtiles:
        cmds = cmd_separator([], f"Resampling mbtiles and converting to pmtiles: {geotif_f}")
        cmds.append(f"cp {mbtile_f} {mbtile_f_ann}")
        cmds.append(f"{args.gdaladdo} {mbtile_f_ann} -r {args.resample} 2 4 8 16 32 64 128 256")
        cmds.append(f"{args.pmtiles} convert --force {mbtile_f_ann} {pmtiles_f}")
        cmds_rm.append(f"rm -f {mbtile_f_ann}")

        if not args.keep_files:
            cmds += cmds_rm ## remove intermediate files

        mm.add_target(pmtiles_f, [mbtile_f], cmds)

    # 4. Upload new PMtiles to AWS
    if args.upload_aws:
        assert args.aws_bucket is not None, "Please provide the AWS S3 bucket path using --aws-bucket"
        cmds = cmd_separator([], f"Uploading pmtiles to AWS: {geotif_f}")
        pmtiles_fn=os.path.basename(pmtiles_f)
        cmds.append(f"aws s3 cp {pmtiles_f} {args.aws_bucket}/{pmtiles_fn}")
        cmds.append(f"touch {pmtiles_f}.done")
        mm.add_target(f"{pmtiles_f}.done", [pmtiles_f], cmds)

    # 5. Update the catalog.yaml file with the new pmtiles and upload to AWS
    if args.update_yaml:
        cmds = cmd_separator([], f"Updating yaml and uploading to AWS: {geotif_f}")
        cmds.append(f"cartloader update_yaml_for_basemap --yaml {catalog_f} --pmtiles {args.basemap_key}:{pmtiles_f}")
        cmds.append(f"aws s3 cp {catalog_f} {args.aws_bucket}/catalog.yaml")
        cmds.append(f"touch {pmtiles_f}.yaml.done")
        mm.add_target(f"{pmtiles_f}.yaml.done", [pmtiles_f, catalog_f], cmds)
    
    ## write makefile
    make_f=os.path.join(out_dir, args.makefn) if args.makefn is not None else f"{args.out_prefix}.mk"
    mm.write_makefile(make_f)

    # run makefile
    # if args.dry_run:
    #     os.system(f"make -f {args.out_dir}/{args.makefn} -n")
    #     print(f"To execute the pipeline, run the following command:\nmake -f {args.out_dir}/{args.makefn} -j {args.n_jobs}")
    # else:
    #     os.system(f"make -f {args.out_dir}/{args.makefn} -j {args.n_jobs}")

    if args.dry_run:
        os.system(f"make -f {make_f} -n")
        print(f"To execute the pipeline, run the following command:\nmake -f {make_f} -j {args.n_jobs}")
    else:
        result = subprocess.run(f"make -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}", shell=True)
        if result.returncode != 0:
            print(f"Error in converting the figure ({args.in_fig}) to pmtiles ({pmtiles_f})")
            sys.exit(1)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])