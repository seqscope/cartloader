# Purpose: This will generate the new json file assuming no external model was used.
# Usage example: 
#   cartloader write_json_for_lda_ficture --sge-type filtered --out-dir /nfs/turbo/sph-hmkang/index/data/weiqiuc/testruns/testcases/umst_mouse_brain/nova_st/out4_filtered

import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json
from collections import defaultdict, OrderedDict

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="""
                                     Write a JSON file to summarize the parameters.""")
    parser.add_argument('--fic-dir', required= True, type=str, help='run_ficture output directory')
    parser.add_argument('--out-json', type=str, default=None, help='Path to the output JSON file. Default: <out-dir>/ficture.params.json')
    parser.add_argument('--sge-type', type=str, default=None, help='(Optional) Type of SGE (options: "raw" or "filtered"). If specified, the script will automatically define the file name for transcript, feature and tsv using the default file name in sge_convert in cartloader. The path will be defined in the out-dir.')
    parser.add_argument('--in-cstranscript', type=str, default=None, help='If --sge-type is not used, provide path to the transcript file used for run_ficture.')
    parser.add_argument('--in-feature', type=str, default=None, help='If --sge-type is not used, provide path to the input feature file used for run_ficture.')
    parser.add_argument('--in-minmax', type=str, default=None, help='If --sge-type is not used, provide path to the input minmax file used for run_ficture.')
    parser.add_argument('--static-cmap-file', type=str, default=None, help='(Optional) If a static cmap was used here, provide the path to the static cmap file.')
    parser.add_argument('--overwrite', action='store_true', default=False, help='Overwrite the existing JSON file.')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args) 

def extract_parameters(filename):
    """Extract parameters from the filename using regular expressions."""
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

def categorize_files(files, out_dir, static_cmap_file):
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
                model_path = os.path.join(out_dir, file.replace(".done", ".model_matrix.tsv.gz"))
                cmap_path = os.path.join(out_dir, file.replace(".done", ".rgb.tsv")) if static_cmap_file is None else static_cmap_file
                train_dict[(n_factor, train_width)].update({
                    "model_type": "lda",
                    "model_id": f"t{train_width}_f{n_factor}",
                    "model_path": model_path,
                    "n_factor": n_factor,
                    "train_width": train_width,
                    "cmap": cmap_path
                })
            elif radius is None:
                train_dict[(n_factor, train_width)]["proj_params"][(fit_width, anchor_res)].update({
                    "proj_id": f"t{train_width}_f{n_factor}_p{fit_width}_a{anchor_res}",
                    "fit_width": fit_width,
                    "anchor_res": anchor_res,
                    "decode_params": []
                })
            else:
                train_dict[(n_factor, train_width)]["proj_params"][(fit_width, anchor_res)]["decode_params"].append({
                    "decode_id": f"t{train_width}_f{n_factor}_p{fit_width}_a{anchor_res}_r{radius}",
                    "radius": radius
                })
    return train_dict

def convert_to_ordered_list_structure(train_dict):
    """Convert the nested dictionary structure to the required ordered list structure."""
    train_params = []
    for (n_factor, train_width), train_entry in train_dict.items():
        proj_params = []
        for (fit_width, anchor_res), proj_entry in train_entry["proj_params"].items():
            proj_params.append(OrderedDict([
                ("proj_id", proj_entry["proj_id"]),
                ("fit_width", proj_entry["fit_width"]),
                ("anchor_res", proj_entry["anchor_res"]),
                ("decode_params", proj_entry["decode_params"])
            ]))
        ordered_train_entry = OrderedDict([
            ("model_type", train_entry["model_type"]),
            ("model_id", train_entry["model_id"]),
            ("model_path", train_entry["model_path"]),
            ("train_width", train_entry["train_width"]),
            ("n_factor", train_entry["n_factor"]),
            ("cmap", train_entry["cmap"]),
            ("proj_params", proj_params)
        ])
        train_params.append(ordered_train_entry)
    return train_params

def define_input_sge(args):
    if args.sge_type is not None:
        if args.sge_type == "raw":
            args.in_cstranscript = os.path.join(args.out_dir, "transcripts.sorted.tsv.gz")
            args.in_feature = os.path.join(args.out_dir, "feature.clean.tsv.gz")
            args.in_minmax = os.path.join(args.out_dir, "coordinate_minmax.tsv")
        elif args.sge_type == "filtered":
            args.in_cstranscript = os.path.join(args.out_dir, "filtered.transcripts.sorted.tsv.gz")
            args.in_feature = os.path.join(args.out_dir, "filtered.feature.lenient.tsv.gz")
            args.in_minmax = os.path.join(args.out_dir, "filtered.coordinate_minmax.tsv")
        else:
            raise ValueError(f"Invalid SGE type: {args.sge_type}")
    assert os.path.exists(args.in_cstranscript), f"Transcript file not found: {args.in_cstranscript}"
    assert os.path.exists(args.in_feature), f"Feature file not found: {args.in_feature}"
    assert os.path.exists(args.in_minmax), f"Minmax file not found: {args.in_minmax}"
    sge_dict = {
        "in_cstranscript": args.in_cstranscript,
        "in_feature": args.in_feature,
        "in_minmax": args.in_minmax
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

def write_json_for_lda_ficture(_args):
    args = parse_arguments(_args)
    args.out_dir = os.path.realpath(args.out_dir)
    if args.out_json is None:
      args.out_json = os.path.join(args.out_dir, "ficture.params.json")
    # Use done files to extract parameters
    flag_paths = os.path.join(args.out_dir, "*.done")
    flag_fn = [os.path.basename(x) for x in glob.glob(flag_paths)]
    # order the files by train, proj, and decode parameters
    flag_fn_ordered = sorted(flag_fn, key=rank_file) 
    train_dict = categorize_files(flag_fn_ordered, args.out_dir, args.static_cmap_file)
    train_params = convert_to_ordered_list_structure(train_dict)
    # define the input files
    in_sge=define_input_sge(args)
    json_data = {
        "in_sge": in_sge,
        "train_params": train_params
    }
    if os.path.exists(args.out_json) and not args.overwrite:
        raise FileExistsError(f"Output JSON file already exists: {args.out_json}. Please use --overwrite to overwrite the file.")
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
