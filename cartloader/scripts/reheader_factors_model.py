import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json, subprocess
from collections import defaultdict, OrderedDict
from cartloader.utils.utils import create_custom_logger

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="Reheader the TSV file output generated with external models")
    parser.add_argument('--input-model', required=True, type=str, help='TSV file for fit_result.tsv.gz')
    parser.add_argument('--out', required=True, type=str, help='Prefix for the output TSV files')
    parser.add_argument('--model-suffix', type=str, default='.model_matrix.tsv.gz', help='Suffix for the output model_matrix TSV file')
    parser.add_argument('--factormap-suffix', type=str, default='.factormap.tsv', help='Suffix for the output factormap TSV file')
    parser.add_argument('--gzip', type=str, default='gzip', help='Path to gzip binary. May be replaced, e.g. to "pigz -p 10", for faster compression')
    parser.add_argument('--log', action='store_true', default=False, help='Write log to file')
    parser.add_argument('--log-suffix', type=str, default=".reheader.log", help='The suffix for the log file (appended to the output directory). Default: .log')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args) 


def reheader_factors_model(_args):
    args = parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    # Read the header of the input model
    logger.info(f"Reading factors from the input model: {args.input_model}")
    with gzip.open(args.input_model, "rt") as f:
        lines = f.readlines()

    header = lines[0].strip()
    factor_list = header.split("\t")[1:]  # Extract factor list (excluding "gene")
    logger.info(f"Factors: {factor_list}")

    if all(x.isdigit() for x in factor_list):
        logger.info(f"Copying the input model to {args.out}{args.model_suffix}...")
        shutil.copyfile(args.input_model, args.out + args.model_suffix)
    else:
        logger.info("Converting to factor indices.")
        # Create a mapping from factor names to indices
        factor_dict = {x: str(i) for i, x in enumerate(factor_list)}

        # Rewrite the file with updated factors in a single pass
        with gzip.open(args.out + args.model_suffix, "wt") as f_out:
            updated_header = "gene\t" + "\t".join(factor_dict.values()) + "\n"
            f_out.write(updated_header)
            f_out.writelines(lines[1:])

        ## store the factormap file
        logger.info(f"Writing the factormap file {args.out}{args.factormap_suffix}...")
        with open(f"{args.out}{args.factormap_suffix}", 'w') as f:
            f.write("Factor\tFactorName\n")
            for factor, index in factor_dict.items():
                f.write(f"{index}\t{factor}\n")

    logger.info(f"Analysis Finished")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
