import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json
from collections import defaultdict, OrderedDict
# from cartloader.utils.minimake import minimake
# from cartloader.utils.utils import cmd_separator, scheck_app

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="Write a JSON file to summarize the parameters.")
    parser.add_argument('--out-dir', required=True, type=str, help='Output directory')
    parser.add_argument('--out-json', type=str, default=None, help='Path to the output JSON file. Default: <out-dir>/ficture.params.json')
    parser.add_argument('--in-cstranscript', type=str, default=None, help='Path to the transcript file.')
    parser.add_argument('--in-feature', type=str, default=None, help='Path to the feature file.')
    parser.add_argument('--in-minmax', type=str, default=None, help='Path to the minmax file.')
    parser.add_argument('--lda-model', nargs='*', type=str, default=None, help='LDA Model information: <model_type>,<model_path>,<model_id>,<train_width>,<n_factor>,<cmap>')
    parser.add_argument('--ext-model', nargs='*', type=str, default=None, help='External Model information: <model_type>,<model_path>,<model_id>,<train_width>,<factor_map>,<cmap>')
    parser.add_argument('--projection', nargs='*', type=str, default=None, help='Projection information: <model_type>,<model_id>,<projection_id>,<fit_width>,<anchor_res>')
    parser.add_argument('--decode', nargs='*', type=str, default=None, help='Decode information: <model_type>,<model_id>,<projection_id>,<decode_id>,<radius>')
    parser.add_argument('--merge', action='store_true', default=False, help='Merge with to the existing JSON file.')
    parser.add_argument('--overwrite', action='store_true', default=False, help='Overwrite the existing JSON file.')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args) 

def write_json(data, filename):
    """Write the given data to a JSON file."""
    with open(filename, 'w') as file:
        json.dump(data, file, indent=4, sort_keys=False)

def merge_params(params1, params2, keynames, nodenames):
    out_params = []
    cur_key_name = keynames[0]
    next_key_names = keynames[1:]
    cur_node_name = nodenames[0]
    next_node_names = nodenames[1:]
    ## sanity check: make sure that cur_key_name exists in each entry
    value2idx = {} ## value -> (idx1, idx2)
    for (i, p) in enumerate(params1):
        if cur_key_name not in p:
            raise ValueError(f"Key '{cur_key_name}' does not exist in the parameter entry: {p}")
        value2idx[p[cur_key_name]] = [i, None]
    for (i, p) in enumerate(params2):
        if cur_key_name not in p:
            raise ValueError(f"Key '{cur_key_name}' does not exist in the parameter entry: {p}")
        if p[cur_key_name] in value2idx:
            value2idx[p[cur_key_name]][1] = i
        else:
            value2idx[p[cur_key_name]] = [None, i]
    ## merge the parameters
    for (key, [idx1, idx2]) in value2idx.items():
        if idx1 is not None and idx2 is not None: ## merge the contents
            if len(next_key_names) > 0:
                mrg_params = merge_params(params1[idx1][cur_node_name], params2[idx2][cur_node_name], next_key_names, next_node_names)
                params1[idx1][cur_node_name] = mrg_params
                params2[idx2][cur_node_name] = mrg_params
                if params1[idx1] != params2[idx2]:
                    raise ValueError(f"Conflicting parameters: {params1[idx1]} vs {params2[idx2]}")
                else:
                    out_params.append(params1[idx1])
            else:
                if params1[idx1] != params2[idx2]:
                    raise ValueError(f"Conflicting parameters: {params1[idx1]} vs {params2[idx2]}")
                else:
                    out_params.append(params1[idx1])
        elif idx1 is not None:
            out_params.append(params1[idx1])
        else:
            out_params.append(params2[idx2])
    return out_params

def write_json_for_ficture(_args):
    args = parse_arguments(_args)
    if args.overwrite and args.append:
        raise ValueError("Cannot use both --overwrite and --append options.")
    if args.out_json is None:
      args.out_json = os.path.join(args.out_dir, "ficture.params.json")
    
    # Input SGE data
    sge_data={
        "in_cstranscript": args.in_cstranscript,
        "in_feature": args.in_feature,
        "in_minmax": args.in_minmax
    }
    # Model data structure
    train_params = []
    model_dict = {}
    # Process model data
    if args.lda_model is not None:
        for model in args.lda_model:
            model_type, model_path, model_id, train_width, n_factor, cmap = model.split(',')
            model_entry = {
                "model_type": model_type,
                "model_id": model_id,
                "model_path": model_path,
                "train_width": int(train_width),
                "n_factor": int(n_factor),
                "cmap": cmap,
                "proj_params": []
            }
            model_dict[(model_type, model_id)] = model_entry
            train_params.append(model_entry)

    if args.ext_model is not None:
        for model in args.ext_model:
            model_type, model_path, model_id, train_width, factor_map, cmap = model.split(',')
            model_entry = {
                "model_type": model_type,
                "model_id": model_id,
                "model_path": model_path,
                "train_width": int(train_width),
                "factor_map": factor_map,
                "cmap": cmap,
                "proj_params": []
            }
            model_dict[(model_type, model_id)] = model_entry
            train_params.append(model_entry)

    # Process projection data
    if args.projection is not None:
        for proj in args.projection:
            model_type, model_id, projection_id, fit_width, anchor_res = proj.split(',')
            projection_entry = {
                "proj_id": projection_id,
                "fit_width": int(fit_width),
                "anchor_res": int(anchor_res),
                "decode_params": []
            }
            if (model_type, model_id) in model_dict:
                model_dict[(model_type, model_id)]["proj_params"].append(projection_entry)
    
    # Process decode data
    if args.decode is not None:
        for dec in args.decode:
            model_type, model_id, projection_id, decode_id, radius = dec.split(',')
            decode_entry = {
                "decode_id": decode_id,
                "radius": int(radius)
            }
            # Find the right model and projection to add decode parameters
            if (model_type, model_id) in model_dict:
                for proj in model_dict[(model_type, model_id)]["proj_params"]:
                    if proj["proj_id"] == projection_id:
                        proj["decode_params"].append(decode_entry)

    # Construct final JSON data
    json_data = {
        "in_sge": sge_data,
        "train_params": train_params
    }

    if os.path.exists(args.out_json):
        if args.overwrite:
            print(f'Overwriting the existing JSON file: {args.out_json}')
        elif args.merge:
            print(f'Merging with the existing JSON file: {args.out_json}')
            with open(args.out_json, 'r') as file:
                existing_data = json.load(file)
            if "in_sge" not in existing_data:
                raise ValueError("Existing JSON file does not have the 'in_sge' key.")
            if "train_params" not in existing_data:
                raise ValueError("Existing JSON file does not have the 'train_params' key.")
            ## The in_sge data must be identical
            if existing_data["in_sge"] != sge_data:
                raise ValueError("The 'in_sge' data in the existing JSON file is different from the new data. NOT compartible and --merge option failed")
            ## Merge the existing data with the new data
            json_data["train_params"] = merge_params(existing_data["train_params"], json_data["train_params"], ["model_id", "proj_id", "decode_id"], ["proj_params", "decode_params", ""])
        else:
            raise FileExistsError(f"Output JSON file already exists: {args.out_json}. Please use --overwrite to overwrite the file, and --merge to merge with existing JSON file")

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
