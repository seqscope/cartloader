import os
import sys
import subprocess
import argparse
import yaml
import datetime
from types import SimpleNamespace

from cartloader.utils.utils import add_param_to_cmd
from cartloader.utils.image_helper import update_orient

histology_suffixes=[".tif", ".tiff", ".png"]
histology_suffixes.extend([suffix.upper() for suffix in histology_suffixes])

aux_env_args = {
        "run_ficture": ['spatula', 'bgzip', "tabix", "gzip", "sort", "sort_mem"],
        "run_cartload_join": ['magick', 'pmtiles', 'gdal_translate', 'gdaladdo', 'tippecanoe', 'spatula'],
        "run_fig2pmtiles": ['pmtiles', 'gdal_translate', 'gdaladdo']
}

aux_params_args = {
     "run_ficture": ['anchor_res', 'radius_buffer', 
                     'csv_colidx_x', 'csv_colidx_y',
                     'hexagon_n_move', 'hexagon_precision', 'min_ct_per_unit_hexagon',
                     'minibatch_size', 'minibatch_buffer',
                     'train_epoch', 'train_epoch_id_len', 'lda_rand_init', 'lda_plot_um_per_pixel',
                     'fit_width', 'fit_precision', 'min_ct_per_unit_fit', 'fit_plot_um_per_pixel',
                     'decode_top_k', 'decode_block_size', 'decode_scale', 'decode_precision', 'decode_plot_um_per_pixel',
                     'merge_max_dist_um', 'merge_max_k', 'merge_max_p',
                     'min_ct_per_feature', 'de_max_pval', 'de_min_fold']
    }

def merge_config(base_config, args, keys, prefix=None):
    """
    Merges parameters from a base configuration dictionary and command-line arguments.
    Args:
        base_config (dict): Dictionary containing default values.
        args (argparse.Namespace): Parsed command-line arguments.
        keys (list): Keys to be merged from args and base_config.
        prefix (str, optional): Prefix to be added to keys in args.
    Returns:
        SimpleNamespace: Merged configuration.
    """
    config = base_config.get(prefix, {}).copy() if prefix else base_config.copy()
    for key in keys:
        val = getattr(args, f"{prefix}_{key}" if prefix else key, None)
        if val is not None:
            config[key] = val
    return SimpleNamespace(**config)

def check_file(file_path):
    if not (os.path.isfile(file_path) or os.path.islink(file_path)):
        print(f"Input file not found: {file_path}")
        sys.exit(1)
    else:
        print(f"Checked: {file_path}")

def define_sge(sgelab, out_dir, sge_dir):
    # file names by sge lab
    tsvfn    = "transcripts.unsorted.tsv.gz" if sgelab == "raw" else "filtered.transcripts.unsorted.tsv.gz"
    cstsvfn  = "transcripts.sorted.tsv.gz" if sgelab == "raw" else "filtered.transcripts.sorted.tsv.gz"
    ftrfn    = "feature.clean.tsv.gz" if sgelab == "raw" else "filtered.feature.lenient.tsv.gz"
    minmaxfn = "coordinate_minmax.tsv" if sgelab == "raw" else "filtered.coordinate_minmax.tsv"

    if sge_dir is None:
        sge_dir = out_dir

    # link files 
    fic_dir=os.path.join(out_dir, "ficture")
    os.makedirs(fic_dir, exist_ok=True)
    if sge_dir != fic_dir:
        for infn in [cstsvfn, tsvfn, ftrfn, minmaxfn]:
            src = os.path.join(sge_dir, infn) # source
            dst = os.path.join(fic_dir, infn) # destination
            if infn == cstsvfn and not os.path.exists(os.path.abspath(src)):
                continue
            if not os.path.isfile(dst) and not os.path.exists(dst):
                os.symlink(src, dst)
                print(f"Creating symlink for {infn}")

    tsv = os.path.join(fic_dir, tsvfn)
    cstsv = os.path.join(fic_dir, cstsvfn)
    ftr = os.path.join(fic_dir, ftrfn)
    minmax = os.path.join(fic_dir, minmaxfn)

    # check if the files exist
    for infile in [ftr, minmax]:
        check_file(infile)
    # tsv or cstsv
    if os.path.isfile(cstsv) or os.path.islink(cstsv):
        check_file(cstsv)
        sge_arg=f"--in-cstranscript {cstsv} --in-feature {ftr} --in-minmax {minmax}"
    elif os.path.isfile(tsv) or os.path.islink(tsv):
        check_file(tsv)
        sge_arg=f"--in-transcript {tsv} --in-feature {ftr} --in-minmax {minmax}" 
    else:
        print(f"Input file not found: {tsv} or {cstsv}")
        sys.exit(1)
    return sge_arg

def define_args_rgb(rgbpath):
    if rgbpath is None:
        rgb_arg = ""
    else:
        check_file(rgbpath)
        rgb_arg = f"--cmap-static --static-cmap-file {rgbpath}"
    return rgb_arg

# def submit_cmd(mkpref, jobname, slurm, submit):
#     slurm_account=slurm.get("account", None)
#     slurm_partition=slurm.get("partition", None)
#     slurm_mail_user=slurm.get("mail_user", None)
#     slurm_cpus_per_task=slurm.get("cpus_per_task", None)
#     slurm_mem=slurm.get("mem", None)
#     slurm_hours=slurm.get("hours", None)
#     if submit:
#         submit_cmd = " ".join([
#             "sbatch",
#             f"--account={slurm_account}",
#             f"--partition={slurm_partition}",
#             f"--mail-user={slurm_mail_user}",
#             "--mail-type=END,FAIL,REQUEUE",
#             f"--job-name={jobname}",
#             f"--output={mkpref}_%j.out",
#             f"--error={mkpref}_%j.err",
#             f"--time={slurm_hours}:00:00",
#             "--nodes=1",
#             "--ntasks=1",
#             f"--cpus-per-task={slurm_cpus_per_task}",
#             f"--mem={slurm_mem}",
#             f"--wrap=\"make -f {mkpref}.mk -j 10\""
#         ])
#         print("Submitting job with command:", submit_cmd)
#         subprocess.run(submit_cmd, shell=True)
#     # else:
#     #     print("debug:", submit_cmd)


def submit_job(jobname, jobpref, cmds, slurm, args):
    """
    Create and optionally submit a SLURM job.
    """
    # Generate SLURM script content
    if args.submit_mode == "local":
        job_content = """#!/bin/bash\n""" + "\n".join(cmds)
    else:
        # slurm_options = {
        #     "account": f"#SBATCH --account={slurm.get('account')}" if slurm.get("account") else "",
        #     "partition": f"#SBATCH --partition={slurm.get('partition')}" if slurm.get("partition") else "",
        #     "mail_user": f"#SBATCH --mail-user={slurm.get('mail_user')}" if slurm.get("mail_user") else "",
        #     "mail_type": "#SBATCH --mail-type=END,FAIL,REQUEUE" if slurm.get("mail_user") else "",
        #     "job_name": f"#SBATCH --job-name={jobname}",
        #     "output": f"#SBATCH --output={jobpref}_%j.out",
        #     "error": f"#SBATCH --error={jobpref}_%j.err",
        #     "time": f"#SBATCH --time={slurm.get('hours', "5")}:00:00",
        #     "nodes": "#SBATCH --nodes=1",
        #     "ntasks": "#SBATCH --ntasks=1",
        #     "cpus_per_task": f"#SBATCH --cpus-per-task={slurm.get('cpus_per_task', 1)}",
        #     "mem": f"#SBATCH --mem={slurm.get('mem', '6500mb')}"
        # }
        # slurm_header = ["#!/bin/bash"] + [value for value in slurm_options.values() if value] # Filter out empty entries and construct the header
        # job_content = "\n".join(slurm_header + cmds)
        # Generate SLURM options
        slurm_options = {
            "account": f"#SBATCH --account={slurm.account}" if getattr(slurm, "account", None) else "",
            "partition": f"#SBATCH --partition={slurm.partition}" if getattr(slurm, "partition", None) else "",
            "mail_user": f"#SBATCH --mail-user={slurm.mail_user}" if getattr(slurm, "mail_user", None) else "",
            "mail_type": "#SBATCH --mail-type=END,FAIL,REQUEUE" if getattr(slurm, "mail_user", None) else "",
            "job_name": f"#SBATCH --job-name={jobname}",
            "output": f"#SBATCH --output={jobpref}_%j.out",
            "error": f"#SBATCH --error={jobpref}_%j.err",
            "time": f"#SBATCH --time={getattr(slurm, 'hours', '5')}:00:00",
            "nodes": "#SBATCH --nodes=1",
            "ntasks": "#SBATCH --ntasks=1",
            "cpus_per_task": f"#SBATCH --cpus-per-task={getattr(slurm, 'cpus_per_task', 1)}",
            "mem": f"#SBATCH --mem={getattr(slurm, 'mem', '6500mb')}"
        }

        slurm_header = ["#!/bin/bash"] + [value for value in slurm_options.values() if value]  # Filter out empty entries and construct the header
        job_content = "\n".join(slurm_header + cmds)

    # Write SLURM script
    job_file = f"{jobpref}.job"
    try:
        with open(job_file, "w") as f:
            f.write(job_content)
        print(f"A job file created: {job_file}")
    except Exception as e:
        print(f"Error writing a job file: {e}")
        return

    # Submit job if requested
    if args.submit:
        try:
            if args.submit_mode == "local":
                result = subprocess.run(f"bash {job_file}", shell=True, check=True, capture_output=True, text=True)
                print("Job executed successfully:\n", result.stdout)
            elif args.submit_mode == "slurm":
                # ?? The minimum required SLURM parameters
                missing_params = [k for k in [ "cpus_per_task", "mem"] if not getattr(slurm, k, None)]
                if missing_params:
                    print(f"Error: Missing required SLURM parameters: {', '.join(missing_params)}")
                    return
                result = subprocess.run(f"sbatch {job_file}", shell=True, check=True, capture_output=True, text=True)
                print("Job submitted successfully:\n", result.stdout)
        except subprocess.CalledProcessError as e:
            print("Error submitting job:\n", e.stderr)

def cmd_run_ficture(run_i, args, env):
    # makefile 
    mkbn = "run_ficture" if run_i["ext_path"] is None else f"run_ficture_{run_i['ext_id']}"

    # sge and rgb
    fic_dir=os.path.join(run_i["out_dir"], "ficture")
    sge_arg = define_sge(run_i["sge_label"], run_i["out_dir"], run_i.get("sge_dir", None))
    cmap_arg = define_args_rgb(run_i["cmap"])
    
    # model & parameters
    assert run_i["train_width"] is not None, "Error: --train-width is required"
    if run_i["ext_path"] is not None:
        assert run_i["ext_id"] is not None, "Error: --ext-id is required when running ficture with an external model"
        assert os.path.exists(run_i["ext_path"]), f"Error: --ext-path is defined with a missing file: {run_i['ext_path']}"
        mod_args = [
            "--main-ext" if not args.init_ext else "--init-ext",
            f"--ext-path {run_i['ext_path']}",
            f"--ext-id {run_i['ext_id']}",
            "--copy-ext-model" if run_i.get("copy_ext_model") else "",
        ]
    else:
        assert run_i['n_factor'] is not None, "Error: --n-factor is required when running standard ficture"
        mod_args = [
            "--main",
            f"--n-factor {run_i['n_factor']}",
        ]

    mod_args.append(f"--train-width {run_i['train_width']}")
    if run_i.get("anchor_res") is not None:
        mod_args.append(f"--anchor-res {run_i['anchor_res']}")
    if run_i.get("radius_buffer") is not None:
        mod_args.append(f"--radius-buffer {run_i['radius_buffer']}")

    # cmd
    ficture_cmd = " ".join([
        "cartloader", "run_ficture",
        f"--makefn {mkbn}.mk",
        f"--out-dir {fic_dir}",
        sge_arg,
        f"--major-axis {run_i.get('major_axis', 'X')}",
        f"--colname-count {run_i['colname_count']}" if run_i['colname_count'] else "",
        " ".join(mod_args),
        "--skip-coarse-report" if args.skip_coarse_report else "",
        cmap_arg,
        "--threads 10",
        "--n-jobs 10",
    ])
    # add aux tools
    ficture_cmd = add_param_to_cmd(ficture_cmd, env, aux_env_args["run_ficture"])
    # add aux parameters
    ficture_aug = merge_config(run_i, args, aux_params_args["run_ficture"], prefix=None)  # merge auxiliary parameters
    ficture_cmd = add_param_to_cmd(ficture_cmd, ficture_aug, aux_params_args["run_ficture"])

    # dry-run
    if args.dry_run:
        ficture_cmd_dry = f"{ficture_cmd} --dry-run"
        print(f"Dry-run: {ficture_cmd_dry}")
        subprocess.run(ficture_cmd_dry, shell=True)
    return ficture_cmd

def cmd_run_cartload_join(run_i, args, env):
    fic_dir = os.path.join(run_i["out_dir"], "ficture")
    cartload_dir = os.path.join(run_i["out_dir"], "cartload")
    os.makedirs(cartload_dir, exist_ok=True)
    cartload_cmd=" ".join([
        "cartloader", "run_cartload_join", 
        "--makefn run_cartload_join.mk", 
        f"--fic-dir {fic_dir}", 
        f"--out-dir {cartload_dir}",
        f"--id {run_i['id']}",
        f"--major-axis {run_i.get('major_axis', 'X')}",
        f"--colname-count {run_i['colname_count']}" if run_i['colname_count'] else "",
        "--n-jobs 10"])
    # add aux tools
    # remove "gdal_translate" from the list of aux_env_args["run_cartload_join"]
    cartload_cmd = add_param_to_cmd(cartload_cmd, env, aux_env_args["run_cartload_join"], underscore2dash=False)
    # dry-run
    if args.dry_run:
        cartload_cmd_dry = f"{cartload_cmd} --dry-run"
        print(f"Dry-run: {cartload_cmd_dry}")
        subprocess.run(cartload_cmd_dry, shell=True)

    return cartload_cmd

def parse_histology_args(histology_args):
    hist_input = []
    for histology in histology_args:
        hist_info = histology.split(";")
        if len(hist_info) < 2:
            raise ValueError("Error: --histology should at least have <histology_type> and <histology_path>.")
        
        keys = [
            "hist_id", "path",
            "transform", 
            "georeference", "georef_tsv", "georef_bounds", 
            "rotate", "flip"
        ]
        
        hist_dict = {key: hist_info[i] if i < len(hist_info) else None for i, key in enumerate(keys)}
        hist_input.append(hist_dict)

    return hist_input

def update_orient_in_histology(histology):
    hist_path = histology["path"]

    rotation = histology.get("rotate", None)
    
    flip = histology.get("flip", None)
    flip_vertical = flip in [ "vertical", "both"]
    flip_horizontal = flip in ["horizontal", "both"]

    new_rot, new_vflip, new_hflip = update_orient(rotation, flip_vertical, flip_horizontal, hist_path)

    # Update the histology dictionary with the best solution
    histology["rotate"] = new_rot
    if new_vflip and new_hflip:
        histology["flip"] = "both"
    elif new_vflip:
        histology["flip"] = "vertical"
    elif new_hflip:
        histology["flip"] = "horizontal"
    else:
        histology["flip"] = None

    return histology

# def fig2pmtiles_transform(histology, hist_prefix):
#     # inpath
#     hist_path = histology["path"]
#     # rotation
#     if histology.get("rotate", None) is not None:
#         if histology.get("rotate") == "90" or histology.get("rotate") == 90:
#             rot_args =" --rotate-clockwise"
#         elif histology.get("rotate") == "270" or histology.get("rotate") == 270:
#             rot_args =" --rotate-counter"        
#     else:
#         rot_args=""
#     # cmd
#     hist_cmd= " ".join([
#             "cartloader", "transform_aligned_histology", 
#             f"--tif {hist_path}",
#             f"--out-prefix {hist_prefix}",
#             rot_args
#         ])
#     return hist_cmd

def cmd_run_fig2pmtiles(run_i, args, env):
    assert len(run_i.get("histology", [])) > 0, "Error: --histology is required when running fig2pmtiles"
    hist_cmds=[]

    cartload_dir=os.path.join(run_i["out_dir"], "cartload")
    catalog_yaml=os.path.join(cartload_dir, "catalog.yaml")
    for histology in run_i.get("histology", []):
        # 1. histology tif to pmtiles
        hist_path = histology["path"]
        # prefix
        hist_inname = os.path.basename(hist_path).replace("_", "-")
        for suffix in histology_suffixes:
            if hist_inname.endswith(suffix):
                hist_inname = hist_inname[:-len(suffix)]
                break
        hist_prefix = os.path.join(cartload_dir, hist_inname)
        # basemap key
        # update the orientation
        histology = update_orient_in_histology(histology)
        # cmds
        hist_cmd= " ".join([
            "cartloader", "run_fig2pmtiles", 
            "--transform" if histology.get("transform", False) else "",
            "--georeference" if histology.get("georeference", False) else "",
            "--geotif2mbtiles", 
            "--mbtiles2pmtiles", 
            f"--update-catalog --basemap-key {histology["hist_id"]}" if os.path.exists(catalog_yaml) else "",
            f"--in-fig {hist_path}",
            f"--out-prefix {hist_prefix}",
            "--flip-vertical" if histology["flip"] in ["vertical", "both"] else "",
            "--flip-horizontal" if histology["flip"] in ["horizontal", "both"] else "",
            f"--rotate {histology['rotate']}" if histology.get("rotate", None) is not None else "",
            f"--in-tsv {histology.get('georef_tsv', None)}" if histology.get("georef_tsv", None) is not None else "",
            f"--in-bounds {histology.get('georef_bounds', None)}" if histology.get("georef_bounds", None) is not None else "",
        ])
        hist_cmd = add_param_to_cmd(hist_cmd, env, aux_env_args["run_fig2pmtiles"], underscore2dash=False)
        hist_cmds.append(hist_cmd)

    if args.dry_run:
        for hist_cmd in hist_cmds:
            hist_cmd_dry = f"{hist_cmd} --dry-run"
            print(f"Dry-run: {hist_cmd_dry}")
            subprocess.run(hist_cmd_dry, shell=True)
    return hist_cmds
    #return hist_cmds + update_catalog_cmds

def cmd_upload_aws(run_i, args, env):
    cartload_dir=os.path.join(run_i["out_dir"], "cartload")
    # Option 1: use "." to locate the file and handle by perl -lane
    # aws_cmd="\n".join([
    #     f"out={cartload_dir}",
    #     f"id={run_i['id']}",
    #     f"aws_bucket={args.aws_bucket}",
    #     "aws s3 cp ${out}/catalog.yaml s3://${aws_bucket}/${id}/catalog.yaml",
    #     "grep -E '\.' ${out}/catalog.yaml | perl -lane 'print $F[$#F]' | xargs -I {} aws s3 cp ${out}/{} s3://${aws_bucket}/${id}/{}"
    # ])
    # Option 2: use cartloader upload_aws_by_catalog.py
    aws_cmd=" ".join([
        "cartloader", "upload_aws_by_catalog",
        f"--in-dir {cartload_dir}",
        f"--s3-dir \"s3://{args.aws_bucket}/{run_i['id']}\""
    ])
    # dry-run
    if args.dry_run:
        aws_cmd_dry = f"{aws_cmd} --dry-run"
        subprocess.run(aws_cmd_dry, shell=True)
    return aws_cmd

def stepinator(_args):
    parser = argparse.ArgumentParser(description="""
    Stepinator: A tool to run steps in cartloader.
    The input, actions, parameters, and tools can be provided in two ways: 1) using a YAML file or 2) using individual command-line arguments.
    It also allows the job execution in two modes: local and slurm.
                                     """)
    # * actions
    action_params = parser.add_argument_group("Actions", "Actions to be performed")
    action_params.add_argument("--run-ficture", action="store_true", help="Run run-ficture in cartloader. Only the main function is executed.")
    action_params.add_argument("--run-cartload-join", action="store_true", help="Run run-cartload-join in cartloader")
    action_params.add_argument("--run-fig2pmtiles", action="store_true", help="Run run-fig2pmtiles in cartloader. This requires provide histology using --in-yaml or --histology.")
    action_params.add_argument("--upload-aws", action="store_true", help="Upload files to AWS S3 bucket. This requires provide AWS bucket name using --in-yaml or --aws-bucket.")
    action_params.add_argument("--copy-ext-model", action="store_true", help="Auxiliary action parameters for run-ficture. Copy external model when running FICTURE with an external model")
    action_params.add_argument("--init-ext", action="store_true", help="Auxiliary action parameters for run-ficture. Only Initialize external model without run main-ext")
    action_params.add_argument("--skip-coarse-report", action="store_true", help="Auxiliary action parameters for run-ficture. Skip coarse report")

    # * mode
    run_params = parser.add_argument_group("Run", "Run mode")
    run_params.add_argument("--submit", action="store_true",  help="Submit the job")
    run_params.add_argument('--submit-mode', type=str,  default="local", choices=["slurm", "local"], help='Specify how the job should be executed. Choose "slurm" for SLURM, "local" for local execution')
    run_params.add_argument("--dry-run", action="store_true", help="Perform a dry run")

    input_params = parser.add_argument_group("Input Configuration", 
        "Provide input data and settings. You can choose one of these methods to define inputs:\n"
        "1. Use --in-yaml and --ids-from-yaml to provide inputs in an input YAML file. This allows to handle more than one runs at once.\n"
        "2. Use individual command-line arguments to specify input files and settings."
    )
    input_params.add_argument("--in-yaml", '-i', type=str, default=None, help="Input yaml file")
    input_params.add_argument("--ids-from-yaml", nargs="*", default=[], help="Define one or more run ID.") #  We suggest to name a run ID as <data_id>_<version_id>.
    input_params.add_argument("--id", type=str, default=None, help="Data ID")
    input_params.add_argument("--out-dir", '-o', type=str, help="The directory")
    input_params.add_argument("--sge-dir",  type=str, default=None, help="SGE directory (Default: <out-dir>/sge)")
    input_params.add_argument("--sge-label", '-s', type=str, default="raw", help="SGE label")
    input_params.add_argument("--colname-count", type=str, default="count", help="")
    input_params.add_argument("--major-axis", type=str, default="X", help="Major axis")
    input_params.add_argument("--ext-path", type=str, default=None, help="(Optional) The path for the external model. Required when running FICTURE with an external model")
    input_params.add_argument("--ext-id", type=str, default=None, help="(Optional) The ID for the external model. Required when running FICTURE with an external model")
    input_params.add_argument("--train-width", '-w', type=str, default=None, help="Train width. Required when running FICTURE")
    input_params.add_argument("--n-factor", '-n', type=str, default=None, help="Number of factors. Only required when running FICTURE with LDA")
    input_params.add_argument("--histology", type=str, nargs="?", default=[], help="""
                              (Optional) Provide histology info as <hist_id>;<path>;<transform>;<georeference>;<georef_tsv>;<georef_bounds>;<rotate_degree>;<flip_direction>. 
                              Only <hist_id> and <path> are required. Supports multiple files; use <hist_id> to distinguish them. (Default: [])
                              Define <transform> and <georeference> by True or False. If georeference=True, specify <georef_tsv> or <georef_bounds>. 
                              Orientation allows rotation (90, 180, 270) and flipping (vertical, horizontal, both). 
                            """)
    input_params.add_argument('--aws-bucket', type=str, default=None, help='AWS bucket name')

    # * tools
    env_params = parser.add_argument_group("Environment and Tools", 
    "Specify paths to environment tools and settings using one of the following methods:\n" 
    "1. Provide all environment settings in the input YAML file (--in-yaml).\n" 
    "2. Provide individual tool paths and environment variables directly as command-line arguments.\n"  
    )
    env_params.add_argument('--slurm-account', type=str, default=None, help='If --submit-mode slurm, provide a SLURM account')
    env_params.add_argument('--slurm-partition', type=str, default=None, help='If --submit-mode slurm, provide a SLURM partition')
    env_params.add_argument('--slurm-mail-user', type=str, default=None, help='If --submit-mode slurm, provide a SLURM mail user')
    env_params.add_argument('--slurm-cpus-per-task', type=int, default=10, help='If --submit-mode slurm, provide a SLURM cpus per task')
    env_params.add_argument('--slurm-mem', type=str, default="65000mb", help='If --submit-mode slurm, provide a SLURM memory')
    # modules & env
    env_params.add_argument('--hpc-modules', type=str, default=None, help='(Optional) HPC modules to load, separated by comma. When a version is required, use the format: <module>/<version>')
    env_params.add_argument('--conda-base', type=str, default=None, help='(Optional) Conda base path')
    env_params.add_argument('--conda-env', type=str, default=None, help='(Optional) Conda environment to activate. If it is installed in the default location, i.e., within the base path, provide the environment name. Otherwise, provide the full path to the environment.')
    # tools
    env_params.add_argument('--spatula', type=str,  default=None,  help='Path to spatula binary. When not provided, it will use the spatula from the submodules.')    
    env_params.add_argument('--pmtiles', type=str,  default=None,  help='Path to pmtiles binary. When not provided, it will use the pmtiles from the submodules.')
    env_params.add_argument('--tippecanoe', type=str,  default=None,  help='Path to tippecanoe binary. When not provided, it will use the tippecanoe from the submodules.')
    env_params.add_argument('--bgzip', type=str, default=None, help='Path to bgzip binary. For faster processing, use "bgzip -@ 4')
    env_params.add_argument('--tabix', type=str, default=None, help='Path to tabix binary')
    env_params.add_argument('--gzip', type=str, default=None, help='Path to gzip binary. For faster processing, use "pigz -p 4"')
    env_params.add_argument('--sort', type=str, default=None, help='Path to sort binary. For faster processing, you may add arguments like "sort -T /path/to/new/tmpdir --parallel=20 -S 10G"')
    env_params.add_argument('--sort-mem', type=str, default=None, help='Memory size for each process')
    env_params.add_argument('--magick', type=str, default=None, help='Path to magick binary')
    env_params.add_argument('--gdal_translate', type=str, default=None, help='Path to gdal_translate binary')
    env_params.add_argument('--gdaladdo', type=str, default=None, help='Path to gdaladdo binary')

    ficture_aux_params = parser.add_argument_group(
        "Auxiliary Parameters for run_ficture", 
        "Parameters for run_ficture, required only if --run-ficture is used with non-default values. Default values are recommended."
    )    
    ficture_aux_params.add_argument("--cmap", '-c', type=str, default=None, help="Only required if the user prefers to use a pre-built color map (Default: None)")
    # input column indexes
    ficture_aux_params.add_argument('--csv-colidx-x',  type=int, default=1, help='Column index for X-axis in the --in-transcript (default: 1)')
    ficture_aux_params.add_argument('--csv-colidx-y',  type=int, default=2, help='Column index for Y-axis in the --in-transcript (default: 2)')
    # segmentation - ficture
    ficture_aux_params.add_argument('--hexagon-n-move', type=int, default=None, help='Level of hexagonal sliding when creating hexagon-indexed SGE in FICTURE compatible format')
    ficture_aux_params.add_argument('--hexagon-precision', type=float, default=None, help='Output precision of hexagon coordinates for FICTURE compatible format')
    ficture_aux_params.add_argument('--min-ct-per-unit-hexagon', type=int, default=None, help='Minimum count per hexagon in hexagon segmentation in FICTURE compatible format')
    # minibatch
    ficture_aux_params.add_argument('--minibatch-size', type=int, default=None, help='Batch size used in minibatch processing')
    ficture_aux_params.add_argument('--minibatch-buffer', type=int, default=None, help='Batch buffer used in minibatch processing')
    # train 
    ficture_aux_params.add_argument('--train-epoch', type=int, default=None, help='Training epoch for LDA model')
    ficture_aux_params.add_argument('--train-epoch-id-len', type=int, default=None, help='Training epoch ID length')
    ficture_aux_params.add_argument('--lda-rand-init', type=int, default=None, help='Number of random initialization during model training')
    ficture_aux_params.add_argument('--lda-plot-um-per-pixel', type=float, default=None, help='Image resolution for LDA plot')
    # fit 
    ficture_aux_params.add_argument('--fit-width',  type=str, default=None, help='Hexagon flat-to-flat width (in um) during model fitting (default: same to train-width)')
    ficture_aux_params.add_argument('--fit-precision', type=float, default=None, help='Output precision of model fitting')
    ficture_aux_params.add_argument('--min-ct-per-unit-fit', type=int, default=None, help='Minimum count per hexagon unit during model fitting')
    ficture_aux_params.add_argument('--fit-plot-um-per-pixel', type=float, default=None, help='Image resolution for fit coarse plot')   # in Scopeflow, this is set to 2
    # decode
    ficture_aux_params.add_argument('--anchor-res', type=int, default=None, help='Anchor resolution for decoding. If absent, run_ficture will use the default value: 4.')
    ficture_aux_params.add_argument('--radius-buffer', type=int, default=None, help='Buffer to radius(=anchor_res + radius_buffer) for pixel-level decoding. If absent, run_ficture will use the default value: 1.')
    ficture_aux_params.add_argument('--decode-top-k', type=int, default=None, help='Top K columns to output in pixel-level decoding results')
    ficture_aux_params.add_argument('--decode-block-size', type=int, default=None, help='Block size for pixel decoding output')
    ficture_aux_params.add_argument('--decode-scale', type=int, default=None, help='Scale parameters for pixel decoding output')
    ficture_aux_params.add_argument('--decode-precision', type=float, default=None, help='Precision of pixel level decoding')
    ficture_aux_params.add_argument('--decode-plot-um-per-pixel', type=float, default=None, help='Image resolution for pixel decoding plot')
    # merge_by_pixel
    ficture_aux_params.add_argument('--merge-max-dist-um', type=float, default=None, help='Maximum distance in um for merging pixel-level decoding results') 
    ficture_aux_params.add_argument('--merge-max-k', type=int, default=None, help='Maximum number of K columns to output in merged pixel-level decoding results')
    ficture_aux_params.add_argument('--merge-max-p', type=int, default=None, help='Maximum number of P columns to output in merged pixel-level decoding results')
    # others parameters shared across steps
    ficture_aux_params.add_argument('--min-ct-per-feature', type=int, default=None, help='Minimum count per feature during LDA training, transform and decoding')
    ficture_aux_params.add_argument('--de-max-pval', type=float, default=None, help='p-value cutoff for differential expression')
    ficture_aux_params.add_argument('--de-min-fold', type=float, default=None, help='Fold-change cutoff for differential expression')
    args = parser.parse_args(_args)

    assert args.run_ficture or args.run_cartload_join or args.upload_aws or args.run_fig2pmtiles, "Error: at least one action is required"

    if args.in_yaml is not None:
        with open(args.in_yaml, "r") as f:
             #yml.update(yaml.safe_load(f))
            yml = yaml.safe_load(f)
        # run info
        if args.ids_from_yaml==["all"]:
            in_ids = [item.get("id") for item in yml.get("id2param", []) if "id" in item] 
        else:
            in_ids = args.ids_from_yaml
        assert len(in_ids) > 0, "Error: --id or --all-ids is required"
        runinfo=[]
        for run_i in yml.get("id2param", []):
            if run_i.get("id") in in_ids:
                runinfo.append(run_i)
    else:
        yml={}
        args.out_dir=os.path.abspath(args.out_dir)
        # create a dictionary for the input    
        runinfo=[
            {   
                # run ID
                "id": args.id,
                # in/out directory
                "sge_dir": args.sge_dir,
                "out_dir": args.out_dir,
                # sge info
                "sge_label": args.sge_label,
                "major_axis": args.major_axis,
                "colname_count": args.colname_count,
                # external model when applied
                "ext_id": args.ext_id,
                "ext_path": args.ext_path,
                "copy_ext_model": args.copy_ext_model,
                "cmap": args.cmap,
                # analysis parameters
                "train_width": args.train_width,
                "n_factor": args.n_factor,
                "histology":[],
            }
        ] 
        # histology
        if len(args.histology) > 0:
            runinfo[0]["histology"].extend(parse_histology_args(args.histology))

    # env
    # * tools (collect from yaml and args)
    aux_env_args = []
    for action, aux_env_args_i in aux_env_args:
        aux_env_args.extend(aux_env_args_i)
    env   = merge_config(yml, args, aux_env_args,  prefix="env")  
    # * tools with submodules
    if env.spatula is None:
        env.spatula  = os.path.join(cartloader_repo, "submodules", "spatula", "bin", "spatula")
    if env.tippecanoe is None:
        env.tippecanoe = os.path.join(cartloader_repo, "submodules", "tippecanoe", "bin", "tippecanoe")

    # * slurm (collect from yaml and args)
    slurm = merge_config(yml, args, ["account", "partition", "mail_user", "cpus_per_task", "mem", "hours"], prefix="slurm")
    
    # * hpc modules
    hpc_modules =  yml.get("env", {}).get("hpc_modules", [])
    if args.hpc_modules:
        hpc_modules = args.hpc_modules.split(",")

    if len(hpc_modules) > 0:
        module_load_cmds="module load "+ " ".join(hpc_modules)
    else:
        module_load_cmds=""

    # * conda env
    conda_base = args.conda_base if args.conda_base else yml.get("env", {}).get("conda", {}).get("base_path", None)
    conda_env  = args.conda_env  if args.conda_env  else yml.get("env", {}).get("conda", {}).get("env_prefix", None)

    conda_cmd=[]
    if conda_base:
        conda_cmd.append(f"source {conda_base}/etc/profile.d/conda.sh")
    if conda_env:
        conda_cmd.append(f"conda activate {conda_env}")

    # * aws 
    if args.aws_bucket is None:
        args.aws_bucket = yml.get("aws_bucket", None)
    
    for run_i in runinfo:
        print("====================================")
        print(f"  ID: {run_i['id']} ")
        print("====================================")
        
        # create out_dir if not exist
        os.makedirs(run_i["out_dir"], exist_ok=True)

        # cmds
        cmds =  []
        cmds.append(module_load_cmds) 
        cmds.extend(conda_cmd)
        if args.run_ficture:
            cmds.append(cmd_run_ficture(run_i, args, env))
        if args.run_cartload_join:
            cmds.append(cmd_run_cartload_join(run_i, args, env))
        if args.run_fig2pmtiles:
            cmds.extend(cmd_run_fig2pmtiles(run_i, args, env))
        if args.upload_aws:
            cmds.append(cmd_upload_aws(run_i, args, env))
        
        # job/err/out file names
        sjobfns=[]
        if args.run_ficture:
            sjobfns.append("ficture")
        if args.run_cartload_join:
            sjobfns.append("cartload")
        if args.run_fig2pmtiles:
            sjobfns.append("fig2pmtiles")
        if args.upload_aws:
            sjobfns.append("aws")

        sjobfn="_".join(sjobfns)
        
        os.makedirs(os.path.join(run_i["out_dir"], "worklog"), exist_ok=True)
        jobpref=os.path.join(run_i["out_dir"], "worklog", sjobfn) 
        submit_job(run_i["id"], jobpref, cmds, slurm, args)

if __name__ == "__main__":
    # Get the path to the cartloader repository
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    #print(sys.argv)
    func(sys.argv[1:])