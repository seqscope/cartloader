import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, csv, yaml, inspect
import pandas as pd
import numpy as np
from openTSNE import TSNE

from cartloader.utils.utils import scheck_app, create_custom_logger, flexopen, unquote_str, smartsort, write_dict_to_file, load_file_to_dict, run_command

from cartloader.utils.color_helper import normalize_rgb, rgb_to_hex

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Perform Leiden Clustering on LDA fit results")

    run_params = parser.add_argument_group("Run Options", "Run options for GNU Make")
    run_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    run_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    inout_params = parser.add_argument_group("Input output parameters", "Input and output parameters")
    inout_params.add_argument('--tsv', type=str, required=True, help='Input TSV file in the LDA output')
    inout_params.add_argument('--out', type=str, required=True, help='Output CSV file containing clustering results')
    inout_params.add_argument('--offset-data', type=int, default=3, help='Offset for the results file (default: 2)')
    inout_params.add_argument('--random-seed', type=int, help='Random seed (default: 42)')
    inout_params.add_argument('--key-ids', type=str, nargs='+',default=['cell_id'], help='Key ID for the results file (default: cell_id)')
    inout_params.add_argument('--drop-columns', type=str, nargs='*', default=['topK', 'topP'], help='Columns to drop from the results file (default: ["topK", "topP"])')
    inout_params.add_argument('--perplexity', type=int, default=50, help='Perplexity for fit-SNE (default: 50)')
    inout_params.add_argument('--exaggeration', type=float, help='Exaggeration for fit-SNE (default: None)')
    inout_params.add_argument('--initialization', type=str, default='spectral', choices=['spectral', 'pca', 'random'], help='Initialization for fit-SNE (default: spectral)')
    inout_params.add_argument('--sep', type=str, default="\t", help='Separator for the results file (default: ,)')
    inout_params.add_argument('--n-jobs', type=int, default=8, help='Number of jobs for fit-SNE (default: 8)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(_args)

    return args

def lda_tsne(_args):
    """
    Perform Leiden clustering on LDA fit results
    """

    # parse argument
    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out + "_lda_tsne" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    ## Read the results file
    logger.info(f"Reading results file {args.tsv}...")
    df_tsv = pd.read_csv(args.tsv, sep="\t", compression="gzip" if args.tsv.endswith('.gz') else None)
    shape_original = df_tsv.shape 

    logger.info(f"Dropping columns {args.drop_columns}...")
    df_out = df_tsv[args.key_ids].copy()
    mat_tsv = np.sqrt(df_tsv.iloc[:, (args.offset_data-1):].drop(columns=args.drop_columns, errors='ignore').to_numpy())
    logger.info(f"Results file has shape {shape_original}, The extracted matrix has shape {mat_tsv.shape}")

    logger.info(f"Perform fit-SNE dimensionality reduction...")
    tsne = TSNE(
        n_components=2,
        perplexity=args.perplexity,
        metric="cosine",
        n_jobs=args.n_jobs,
        random_state=args.random_seed,
        negative_gradient_method="fft",  # FIt-SNE interpolation
        initialization=args.initialization,
        learning_rate="auto",
        exaggeration=args.exaggeration,
    )
    embedding = tsne.fit(mat_tsv)

    df_out['TSNE1'] = embedding[:, 0]
    df_out['TSNE2'] = embedding[:, 1]

    ## Save the results, if the file name ends with .gz, compress with gzip
    if args.out.endswith(".gz"):
        df_out.to_csv(args.out, sep=args.sep, index=False, compression="gzip")
    else:
        df_out.to_csv(args.out, sep=args.sep, index=False)
    
    logger.info(f"Analysis completed. Output saved to {args.out}")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
