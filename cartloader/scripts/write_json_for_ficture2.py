import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json
from collections import defaultdict, OrderedDict
# from cartloader.utils.minimake import minimake
# from cartloader.utils.utils import cmd_separator, scheck_app

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="Write a JSON file to summarize the parameters.")
    parser.add_argument('--out-json', type=str, default=None, required=True, help='Path to the output JSON file (recommended naming scheme and directory: "ficture.params.json" in the directory hosting FICTURE results)')
    parser.add_argument('--in-transcript', type=str, default=None, help='Path to the transcript file.')
    parser.add_argument('--in-feature', type=str, default=None, help='Path to the feature file.')
    parser.add_argument('--in-minmax', type=str, default=None, help='Path to the minmax file.')
    parser.add_argument('--in-feature-ficture', type=str, default=None, help='(Optional) If FICTURE used a different feature file than the in-feature file, specify the path to the feature file used for FICTURE analysis.')
    parser.add_argument('--lda-model', nargs='*', type=str, default=None, help='LDA Model information: <model_type>,<model_path>,<model_id>,<train_width>,<n_factor>,<cmap>')
    parser.add_argument('--decode', nargs='*', type=str, default=None, help='Projection information: <model_type>,<model_id>,<projection_id>,<fit_width>,<anchor_res>')
    parser.add_argument('--umap', action='store_true', default=False, help='UMAP information if exists: <model_id>,<umap_tsv>,<umap_png>,<umap_single_factor_png>.')
    parser.add_argument('--merge', action='store_true', default=False, help='If enabled and the output JSON already exists, integrates new input into the existing file. If not set, the script will write a new JSON file or overwrite the existing JSON file.')
    parser.add_argument('--merge-override', action='store_true', default=False, help='When used with --merge, allows new input arguments (--in-transcript, --in-feature, --in-minmax, --in-feature-ficture) to override conflicting values in the existing profile. If not set and a conflict occurs, the script will raise an error.')

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

def write_json_for_ficture2(_args):
    args = parse_arguments(_args)

    # Existing files:
    if args.merge_override and not args.merge:
        raise ValueError("--merge-override requires --merge")

    # Input SGE data
    new_sge={
        "in_transcript": args.in_transcript,
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
                "decode_params": []
            }
            model_dict[(model_type, model_id)] = model_entry
            train_params.append(model_entry)

    # Process decode data
    if args.decode is not None:
        for dec in args.decode:
            model_type, model_id, decode_id, fit_width, anchor_res = dec.split(',')
            decode_entry = {
                "decode_id": decode_id,
                "fit_width": int(fit_width),
                "anchor_res": int(anchor_res),
            }
            if (model_type, model_id) in model_dict:
                model_dict[(model_type, model_id)]["decode_params"].append(decode_entry)
            else:
                warnings.warn(f"No matching model for decode entry: {dec}")  # Warning if model not found

    def update_single_value(old_val, new_val, error_msg, override):
        """Update a single scalar value with merge override logic."""
        if old_val == new_val:
            return old_val  # No change

        if override:
            return new_val

        # Only reach here if values differ and override is False
        if (old_val is None) != (new_val is None):
            # One is None and the other is not
            raise ValueError(error_msg)

        raise ValueError(error_msg)

    def finalize_sge(existing_sge, new_sge, override):
        """Resolve and update the in_sge dictionary."""
        final_sge = {}
        for key in ["in_transcript", "in_feature", "in_minmax"]:
            old_val = existing_sge.get(key)
            new_val = new_sge.get(key)
            error_msg = (
                f"The '{key}' in 'in_sge' data in the existing JSON file is different from your input arguments. "
                f"Disable --merge or enable --merge-override to proceed."
            )
            final_sge[key] = update_single_value(old_val, new_val, error_msg, override)
        return final_sge

    # Loading and merging existing JSON data if exists and --merge
    if os.path.exists(args.out_json) and args.merge:
        print(f"Merging with existing JSON: {args.out_json}")
        with open(args.out_json, 'r') as f:
            existing_data = json.load(f)

        # Merge 'in_sge'
        sge_data = finalize_sge( existing_data.get("in_sge", {}), new_sge, args.merge_override)

        # Merge 'in_feature_ficture'
        existing_ftr_ficture = existing_data.get("in_feature_ficture")
        new_ftr_ficture = args.in_feature_ficture
        error_msg = (
            "The 'in_feature_ficture' in the existing JSON file is different from your input arguments. "
            "Disable --merge or enable --merge-override to proceed."
        )
        ftr_ficture_data = update_single_value(existing_ftr_ficture, new_ftr_ficture, error_msg, args.merge_override)

        # Merge model params
        train_params = merge_params(
            existing_data["train_params"],
            train_params,
            ["model_id", "decode_id"],
            ["decode_params", ""]
        )
    else:
        if os.path.exists(args.out_json):
            print(f"Overwriting the existing {args.out_json}.")
        sge_data = new_sge
        ftr_ficture_data = args.in_feature_ficture


    # Construct final JSON data
    if ftr_ficture_data is None:
        json_data = {
            "in_sge": sge_data,
            "train_params": train_params
        }
    else:
        json_data = {
            "in_sge": sge_data,
            "in_feature_ficture": ftr_ficture_data,
            "train_params": train_params
        }
 
    write_json(json_data, args.out_json)
    print(f'Data has been written to {args.out_json}')

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
