import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json
from collections import defaultdict, OrderedDict
# from cartloader.utils.minimake import minimake
# from cartloader.utils.utils import cmd_separator, scheck_app

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="""
                                     Write a JSON file to summarize the parameters.""")
    parser.add_argument('--out-dir', required= True, type=str, help='Output directory')
    parser.add_argument('--out-json', type=str, default=None, help='Path to the output JSON file. Default: <out-dir>/ficture.params.json')
    parser.add_argument('--use-external', action='store_true', help='Use external model.')
    parser.add_argument('--external-model', type=str, default=None, help='If --use-external, provide an external matrix model to use for projection.')
    parser.add_argument('--external-model-flag', type=str, default=None, help='(Optional) Specify a flag to mark the output file for projection using an external model. If not provided, the flag will be automatically generated based on the first part of the model file name split by "."')
    parser.add_argument('--sge-type', type=str, default=None, help='(Optional) Type of SGE (options: "raw" or "filtered"). If specified, the script will automatically define the file name for transcript, feature and tsv using the default file name in sge_convert in cartloader. The path will be defined in the out-dir.')
    parser.add_argument('--cstranscript', type=str, default=None, help='If --sge-type is not used, provide path to the transcript file used for run_ficture.')
    parser.add_argument('--feature', type=str, default=None, help='If --sge-type is not used, provide path to the input feature file used for run_ficture.')
    parser.add_argument('--minmax', type=str, default=None, help='If --sge-type is not used, provide path to the input minmax file used for run_ficture.')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args) 

def extract_parameters(filename):
    """Extract parameters from the filename using regular expressions."""
    #pattern = re.compile(r'nF(\d+)\.d_(\d+)(?:\.decode)?(?:\.prj_(\d+)\.r_(\d+)(?:_(\d+))?)?\.bulk_chisq\.tsv')
    pattern = re.compile(r't(\d+)_f(\d+)(?:_p(\d+)_a(\d+))?(?:_r(\d+))?\.done')
    match = pattern.search(filename)
    if match:
        train_width, n_factor, fit_width, anchor_res, radius = match.groups()
        return {
            "n_factor": int(n_factor),
            "train_width": int(train_width),
            "fit_width": int(fit_width) if fit_width else None,
            "anchor_res": int(anchor_res) if anchor_res else None,
            "radius": int(radius) if radius else None,
        }
    return None

def categorize_files(files):
    """Categorize files into train, proj, and decode based on parameters."""
    train_dict = defaultdict(lambda: {"proj_params": defaultdict(lambda: {"decode_params": []})})
    
    for file in files:
        params = extract_parameters(file)
        if params:
            n_factor = params["n_factor"]
            train_width = params["train_width"]
            fit_width = params["fit_width"]
            anchor_res = params["anchor_res"]
            radius = params["radius"]
            
            if fit_width is None and anchor_res is None and radius is None:
                train_dict[(n_factor, train_width)]["n_factor"] = n_factor
                train_dict[(n_factor, train_width)]["train_width"] = train_width
            elif radius is None:
                train_dict[(n_factor, train_width)]["proj_params"][(fit_width, anchor_res)].update({
                    "fit_width": fit_width,
                    "anchor_res": anchor_res,
                    "decode_params": []
                })
            else:
                train_dict[(n_factor, train_width)]["proj_params"][(fit_width, anchor_res)]["decode_params"].append({"radius": radius})
    
    return train_dict

def convert_to_ordered_list_structure(train_dict):
    """Convert the nested dictionary structure to the required ordered list structure."""
    train_params = []
    for (n_factor, train_width), train_entry in train_dict.items():
        proj_params = []
        for (fit_width, anchor_res), proj_entry in train_entry["proj_params"].items():
            proj_params.append(OrderedDict([
                ("fit_width", proj_entry["fit_width"]),
                ("anchor_res", proj_entry["anchor_res"]),
                ("decode_params", proj_entry["decode_params"])
            ]))
        ordered_train_entry = OrderedDict([
            ("train_width", train_entry["train_width"]),
            ("n_factor", train_entry["n_factor"]),
            ("proj_params", proj_params)
        ])
        train_params.append(ordered_train_entry)
    return train_params

# def convert_to_list_structure(train_dict):
#     """Convert the nested dictionary structure to the required list structure."""
#     train_params = []
#     for (n_factor, train_width), train_entry in train_dict.items():
#         proj_params = []
#         for (fit_width, anchor_res), proj_entry in train_entry["proj_params"].items():
#             proj_params.append(proj_entry)
#         train_entry["proj_params"] = proj_params
#         train_params.append(train_entry)
#     return train_params

# def order_files(files):
#     train_files=[x for x in files if "prj" not in x]
#     proj_files=[x for x in files if "prj" in x and "decode" not in x] 
#     decode_files=[x for x in files if "decode" in x]
#     ordered_files=train_files+proj_files+decode_files
#     return ordered_files   

def define_input_sge(args):
    if args.sge_type is not None:
        if args.sge_type == "raw":
            args.cstranscript = os.path.join(args.out_dir, "transcripts.sorted.tsv.gz")
            args.feature = os.path.join(args.out_dir, "feature.clean.tsv.gz")
            args.minmax = os.path.join(args.out_dir, "coordinate_minmax.tsv")
        elif args.sge_type == "filtered":
            args.cstranscript = os.path.join(args.out_dir, "filtered.transcripts.sorted.tsv.gz")
            args.feature = os.path.join(args.out_dir, "filtered.feature.lenient.tsv.gz")
            args.minmax = os.path.join(args.out_dir, "filtered.coordinate_minmax.tsv")
        else:
            raise ValueError(f"Invalid SGE type: {args.sge_type}")
    assert os.path.exists(args.cstranscript), f"Transcript file not found: {args.cstranscript}"
    assert os.path.exists(args.feature), f"Feature file not found: {args.feature}"
    assert os.path.exists(args.minmax), f"Minmax file not found: {args.minmax}"
    sge_dict = {
        "tsv": args.cstranscript,
        "feature": args.feature,
        "minmax": args.minmax
    }
    return sge_dict
        

def rank_file(file_name):
    parts = file_name.split('_')
    if len(parts) == 2:  # Only t*_f*.done
        return 1
    elif len(parts) == 4:  # t*_f*_p*_a*.done
        return 2
    elif len(parts) == 5:  # t*_f*_p*_a*_r*.done
        return 3
    else:
        return 4  # If it doesn't fit the above patterns

def write_json(data, filename):
    """Write the given data to a JSON file."""
    with open(filename, 'w') as file:
        json.dump(data, file, indent=4, sort_keys=False)

def write_json_for_ficture(_args):
    args = parse_arguments(_args)
    if args.out_json is None:
      args.out_json = os.path.join(args.out_dir, "ficture.params.json")
    # Note we used done files here
    if args.use_external:
        flag_paths = os.path.join(args.out_dir, "*.done")
        # append the 
    else:
        flag_paths = os.path.join(args.out_dir, "*.done")
        flag_fn = [os.path.basename(x) for x in glob.glob(flag_paths)]
        flag_fn = [x.replace(".done", "") for x in flag_fn]
    # order the files by train, proj, and decode parameters
    flag_fn_ordered = sorted(flag_fn, key=rank_file) 
    train_dict = categorize_files(flag_fn_ordered)
    train_params = convert_to_ordered_list_structure(train_dict)
    # define the input files
    in_sge=define_input_sge(args)
    json_data = {
        "input_sge": in_sge,
        "train_params": train_params
    }
    write_json(json_data, args.out_json)
    print(f'Data has been written to {args.out_json}')

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
