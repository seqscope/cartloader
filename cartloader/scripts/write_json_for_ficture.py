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
    
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args) 

def extract_parameters(filename):
    """Extract parameters from the filename using regular expressions."""
    pattern = re.compile(r'nF(\d+)\.d_(\d+)(?:\.decode)?(?:\.prj_(\d+)\.r_(\d+)(?:_(\d+))?)?\.bulk_chisq\.tsv')
    match = pattern.search(filename)
    if match:
        n_factor, train_width, fit_width, anchor_res, radius = match.groups()
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

def order_files(files):
    train_files=[x for x in files if "prj" not in x]
    proj_files=[x for x in files if "prj" in x and "decode" not in x] 
    decode_files=[x for x in files if "decode" in x]
    ordered_files=train_files+proj_files+decode_files
    return ordered_files   

def write_json(data, filename):
    """Write the given data to a JSON file."""
    with open(filename, 'w') as file:
        json.dump(data, file, indent=4, sort_keys=False)

def write_json_for_ficture(_args):
    args = parse_arguments(_args)
    if args.out_json is None:
      args.out_json = os.path.join(args.out_dir, "ficture.params.json")
    # Note we used de files here
    all_de_paths = os.path.join(args.out_dir, "*.bulk_chisq.tsv")
    all_de_fn = [os.path.basename(x) for x in glob.glob(all_de_paths)]
    # order the files by train, proj, and decode parameters
    all_de_fn_ordered = order_files(all_de_fn)
    train_dict = categorize_files(all_de_fn_ordered)
    train_params = convert_to_ordered_list_structure(train_dict)
    json_data = {
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
