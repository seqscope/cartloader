import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json, subprocess
from collections import defaultdict, OrderedDict
# from cartloader.utils.minimake import minimake
# from cartloader.utils.utils import cmd_separator, scheck_app
from cartloader.utils.utils import create_custom_logger

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="Reheader the TSV file output generated with external models")
    parser.add_argument('--fit-tsv', required=True, type=str, help='TSV file for fit_result.tsv.gz')
    parser.add_argument('--postcount-tsv', required=True, type=str, help='TSV file for posterior_count.tsv.gz')
    parser.add_argument('--out', required=True, type=str, help='Prefix for the output TSV files')
    parser.add_argument('--fit-suffix', type=str, default='.fit_result.tsv.gz', help='Suffix for the output fit_results TSV file')
    parser.add_argument('--model-suffix', type=str, default='.model_matrix.tsv.gz', help='Suffix for the output model_matrix TSV file')
    parser.add_argument('--postcount-suffix', type=str, default='.posterior.count.tsv.gz', help='Suffix for the output posterior_count TSV file')
    parser.add_argument('--factormap-suffix', type=str, default='.factormap.tsv', help='Overwrite the existing TSV files')
    parser.add_argument('--gzip', type=str, default='gzip', help='Path to gzip binary. May be replaced, e.g. to "pigz -p 10", for faster compression')
    parser.add_argument('--log', action='store_true', default=False, help='Write log to file')
    parser.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args) 

def reheader_factors_tsv(_args):
    args = parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    ## read the header of the two files and check the consistency 
    ## postcount tsv contains "gene" + factors
    tmp_hdr_fit = f"{args.out}{args.fit_suffix}.tmp.hdr.tsv"

    ## The fit_result.tsv.gz preserves the original factor names while posteriot_count.tsv.gz does not
    ## We will reorder the posterior counts to match the order of the factors in the fit_results
    
    ## Writing the header files for fit_tsv files
    factor_names = []
    n_xtra_col = 0
    n_factor = 0
    logger.info("Writing the header of the new fit_results file")
    with gzip.open(args.fit_tsv, 'rt') as f:
        with open(tmp_hdr_fit, 'w') as f_out:
            fit_header = f.readline().strip().split('\t')
            n_fit = len(fit_header)
            found_topP = False      ## topP is the last extra column
            for i in range(n_fit):
                if found_topP:
                    f_out.write(f"\t{len(factor_names)}")
                    factor_names.append(fit_header[i])
                    n_factor += 1
                else:
                    if i > 0:
                        f_out.write("\t")
                    f_out.write(fit_header[i])
                    n_xtra_col += 1
                    if fit_header[i] == "topP":
                        found_topP = True
            f_out.write("\n")

    name2idx = {factor_names[i]:i for i in range(n_factor)}
    assert len(name2idx) == n_factor, "Duplicate factor names found"

    ## Reorder and rewrite the posterior count file    
    logger.info("Writing new posterior count and model matrix files with new header")
    with gzip.open(args.postcount_tsv, 'rt') as f:
        with gzip.open(f"{args.out}{args.postcount_suffix}", 'wt') as out1:
            with gzip.open(f"{args.out}{args.model_suffix}", 'wt') as out2:
                idx_orders = []
                for line in f:
                    fields = line.strip().split('\t')
                    if len(idx_orders) == 0:
                        field2idx = {fields[i]:i for i in range(1,len(fields))}
                        for i, name in enumerate(factor_names):
                            if name in field2idx:
                                idx_orders.append(field2idx[name])
                                fields[field2idx[name]] = i
                            else:
                                raise ValueError(f"Factor '{name}' not found in the posterior_count file")
                        #print(idx_orders)
                    out1.write(fields[0])
                    out2.write(fields[0])
                    for i in idx_orders:
                        out1.write(f"\t{fields[i]}")
                        out2.write(f"\t{fields[i]}")
                    out1.write("\n")
                    out2.write("\n")

    ## store the factormap file
    logger.info(f"Writing the factormap file {args.out}{args.factormap_suffix}...")
    with open(f"{args.out}{args.factormap_suffix}", 'w') as f:
        f.write("Factor\tFactorName\n")
        for i in range(n_factor):
            f.write(f"{i}\t{factor_names[i]}\n")

    ## Create the new fit file
    logger.info(f"Writing the new fit_results file {args.out}{args.fit_suffix}...")
    cmd = f"(cat {tmp_hdr_fit}; gzip -cd {args.fit_tsv} | tail -n +2;)| {args.gzip} > {args.out}{args.fit_suffix}"
    res_fit = subprocess.run(cmd, shell=True, check=True)

    ## remove the temporary files
    logger.info(f"Removing temporary files...")
    os.remove(tmp_hdr_fit)

    logger.info(f"Analysis Finished")

if __name__ == "__main__":
    # get the cartloader path
    global cartloader_repo
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
