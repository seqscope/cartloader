import sys, os, gzip, argparse, logging, warnings, shutil, subprocess
import pandas as pd

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, find_major_axis, add_param_to_cmd, create_symlink

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader run_ficture", 
                                     description="""
                                     Helper script to prepare a directory for running FICTURE with external model and target data.
                                     """)
    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--model', type=str, default=None, help='(Optional) Path to the external model file')
    key_params.add_argument('--train-width', type=str, default="0", help='(Optional) If the external model has a training width, specify it here. Default: 0')
    key_params.add_argument('--target-dir', type=str, default=None, help='Path to the directory containing the projection files. We assume the cartloader run_ficture has been applied to the target data, including sorttsv, minibatch')
    key_params.add_argument('--target-cstranscript', type=str, default=None, help='Path to the coordinate-sorted transcript-indexed SGE file in TSV format')
    key_params.add_argument('--target-minmax', type=str, default=None, help='Path to the coordinate minmax TSV file')
    key_params.add_argument('--out-dir', required= True, type=str, help='Output directory')
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)


def prepare_external_modelmat(_args):
    # args
    args=parse_arguments(_args)

    # input/output
    # 1. input external model
    if args.model is not None:
        # * start with read the number of factors from the model file
        print(f"Reading the number of factors from the model file {args.model}...")
        model_header = pd.read_csv(args.model, sep="\t", nrows=1)
        n_factor = len(model_header.columns) - 1
        print(f"Number of factors in the model: {n_factor}")
        # * create a soft link from the model file to the output directory
        model_prefix=f"t{args.train_width}_f{n_factor}"
        model_stdfn = f"{model_prefix}.model_matrix.tsv.gz"
        model_ln = os.path.join(args.out_dir, model_stdfn)
        create_symlink(args.model, model_ln)

    # 2. output dirs and files
    os.makedirs(args.out_dir, exist_ok=True)

    # 3. link the target files for projection from the target directory to the output directory
    # projection files: assume it has performed sorttsv, minibatch, 
    if args.target_cstranscript is None:
        args.target_cstranscript = os.path.join(args.target_dir, "transcripts.sorted.tsv.gz")

    if args.target_minmax is None:
        args.target_minmax = os.path.join(args.target_dir, "coordinate_minmax.tsv")

    in_cstranscript=f"{args.out_dir}/transcripts.sorted.tsv.gz"
    in_minmax=f"{args.out_dir}/coordinate_minmax.tsv"

    # create symlinks from the target directory to the output directory
    target2input ={
        args.target_cstranscript:in_cstranscript,
        args.target_minmax:in_minmax,
        f"{args.target_dir}/batched.matrix.tsv.gz":f"{args.out_dir}/batched.matrix.tsv.gz"          
    }

    for source, dest in target2input.items():
        create_symlink(source, dest)
    
    # create the done flag file to indicate the completion of the model matrix preparation
    done_flag = os.path.join(args.out_dir, f"{model_prefix}.done")
    with open(done_flag, "w") as f:
        f.write("")

if __name__ == "__main__":
    # Get the path to the cartloader repository
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])