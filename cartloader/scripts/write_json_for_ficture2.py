import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json, copy
from collections import defaultdict, OrderedDict
# from cartloader.utils.minimake import minimake

from cartloader.utils.utils import write_json, reconcile_field

# def update_single_value(old_val, new_val, error_msg, override):
#     """Update a single scalar value with merge override logic."""
#     if old_val == new_val:
#         return old_val  # No change

#     if override:
#         return new_val

#     # Only reach here if values differ and override is False
#     if (old_val is None) != (new_val is None):
#         # One is None and the other is not
#         raise ValueError(error_msg)

#     raise ValueError(error_msg)

# def finalize_sge(existing_sge, new_sge, override):
#     """Resolve and update the in_sge dictionary."""
#     final_sge = {}
#     for key in ["in_transcript", "in_feature", "in_minmax"]:
#         old_val = existing_sge.get(key)
#         new_val = new_sge.get(key)
#         error_msg = (
#             f"The '{key}' in 'in_sge' data in the existing JSON file is different from your input arguments. "
#             f"Disable --merge or enable --merge-override to proceed."
#         )
#         final_sge[key] = update_single_value(old_val, new_val, error_msg, override)
#     return final_sge

# def merge_params(params1, params2, keynames, nodenames):
#     out_params = []
#     cur_key_name = keynames[0]
#     next_key_names = keynames[1:]
#     cur_node_name = nodenames[0]
#     next_node_names = nodenames[1:]
#     ## sanity check: make sure that cur_key_name exists in each entry
#     value2idx = {} ## value -> (idx1, idx2)
#     for (i, p) in enumerate(params1):
#         if cur_key_name not in p:
#             raise ValueError(f"Key '{cur_key_name}' does not exist in the parameter entry: {p}")
#         value2idx[p[cur_key_name]] = [i, None]
#     for (i, p) in enumerate(params2):
#         if cur_key_name not in p:
#             raise ValueError(f"Key '{cur_key_name}' does not exist in the parameter entry: {p}")
#         if p[cur_key_name] in value2idx:
#             value2idx[p[cur_key_name]][1] = i
#         else:
#             value2idx[p[cur_key_name]] = [None, i]
    
#     ## merge the parameters
#     for (key, [idx1, idx2]) in value2idx.items():
#         if idx1 is not None and idx2 is not None: ## merge the contents
#             if len(next_key_names) > 0:
#                 mrg_params = merge_params(params1[idx1][cur_node_name], params2[idx2][cur_node_name], next_key_names, next_node_names)
#                 params1[idx1][cur_node_name] = mrg_params
#                 params2[idx2][cur_node_name] = mrg_params
#                 if params1[idx1] != params2[idx2]:
#                     raise ValueError(f"Conflicting parameters: {params1[idx1]} vs {params2[idx2]}")
#                 else:
#                     out_params.append(params1[idx1])
#             else:
#                 if params1[idx1] != params2[idx2]:
#                     raise ValueError(f"Conflicting parameters: {params1[idx1]} vs {params2[idx2]}")
#                 else:
#                     out_params.append(params1[idx1])
#         elif idx1 is not None:
#             out_params.append(params1[idx1])
#         else:
#             out_params.append(params2[idx2])
#     return out_params

def _normalize_path(path):
    return os.path.abspath(path) if path is not None else None

def _needs_feature_entry(feature_path, sge_feature_path):
    if feature_path is None:
        return False
    if sge_feature_path is None:
        return True
    return _normalize_path(feature_path) != _normalize_path(sge_feature_path)

def build_sge_data(args, old_data):
    sge_keys = ["in_transcript", "in_feature", "in_minmax"]
    new_sge = {key: getattr(args, key) for key in sge_keys}

    # --- write mode or append mode without existing json file ---
    if args.mode == "write" or (args.mode == "append" and not os.path.exists(args.out_json)):
        for key in sge_keys:
            path = new_sge[key]
            assert path is not None, f"Path not provided: --{key.replace('_', '-')}"
            assert os.path.exists(path), f"File not found: {path} ( --{key.replace('_', '-')})"
        return new_sge

    # --- append mode with existing json file ---
    old_sge = old_data.get("in_sge", {})

    final_sge = {}

    for key in sge_keys:
        new_path = new_sge[key]
        old_path = old_sge.get(key)

        if new_path is None:
            # User didn't override -> keep old value (ensure exists)
            assert old_path is not None, f"Path for --{key.replace('_', '-')} not provided and not present in existing JSON."
            assert os.path.exists(old_path), f"File not found: {old_path} (from existing JSON for --{key.replace('_', '-')})"
            final_sge[key] = old_path
        else:
            # User provided a path; check it exists
            assert os.path.exists(new_path), ( f"File not found: {new_path} ( --{key.replace('_', '-')})")
            if old_path is not None:
                # Require consistency with existing path
                assert os.path.abspath(new_path) == os.path.abspath(old_path), f"Found inconsistent absolute path for '{key}' SGE between existing json ({old_path}) and the input arguments ({new_path})."
                final_sge[key] = old_path
            else:
                # No old value; accept the new one
                final_sge[key] = new_path

    return final_sge

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="Write a JSON file to summarize the parameters.")
    parser.add_argument('--out-json', type=str, default=None, required=True, help='Path to the output JSON file (recommended naming scheme and directory: "ficture.params.json" in the directory hosting FICTURE results)')
    parser.add_argument('--mode', type=str, default="write", choices=["write", "append"], help='Write mode for the output JSON. Default: write. If write, a new file will be created based on the arguments provided. If append, the new parameters will be merged into the existing JSON file if it exists.')
    parser.add_argument('--in-transcript', type=str, default=None, help='Path to the transcript file. Required if --mode is write.')
    parser.add_argument('--in-feature', type=str, default=None, help='Path to the feature file. Required if --mode is write.')
    parser.add_argument('--in-minmax', type=str, default=None, help='Path to the minmax file. Required if --mode is write.')
    
    parser.add_argument('--in-feature-ficture', type=str, default=None, help='(Optional) If FICTURE used a different feature file than the in-feature file, specify the path to the feature file used for FICTURE analysis.')
    
    parser.add_argument('--lda-model', nargs='*', type=str, default=None, help='LDA Model information: <model_type>,<model_path>,<model_id>,<train_width>,<n_factor>,<cmap>')
    parser.add_argument('--decode', nargs='*', type=str, default=None, help='Projection information: <model_type>,<model_id>,<projection_id>,<fit_width>,<anchor_res>')
    parser.add_argument('--umap', nargs='*', type=str, default=None, help='UMAP information if exists. Each entry: <model_type>,<model_id>,<umap_tsv>,<umap_png>,<umap_single_factor_png>.')    
    # parser.add_argument('--merge', action='store_true', default=False, help='If enabled and the output JSON already exists, integrates new input into the existing file. If not set, the script will write a new JSON file or overwrite the existing JSON file.')
    # parser.add_argument('--merge-override', action='store_true', default=False, help='When used with --merge, allows new input arguments (--in-transcript, --in-feature, --in-minmax, --in-feature-ficture) to override conflicting values in the existing profile. If not set and a conflict occurs, the script will raise an error.')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args) 

def write_json_for_ficture2(_args):

    args = parse_arguments(_args)

    old_data = {}
    old_train_params = []
    old_feature_ficture = None
    model_dict = {}
    
    # (optional) load data when append mode with existing json file
    if args.mode == "append" and os.path.exists(args.out_json):
        with open(args.out_json, "r") as f:
            old_data = json.load(f)
        old_feature_ficture = old_data.get("in_feature_ficture")
        old_train_params = copy.deepcopy(old_data.get("train_params", []))

    # 1. Construct SGE data
    sge_data = build_sge_data(args, old_data)
    sge_feature_path = sge_data.get("in_feature")

    train_params = []
    legacy_feature = old_feature_ficture if _needs_feature_entry(old_feature_ficture, sge_feature_path) else None

    for entry in old_train_params:
        feature_present = "feature" in entry
        feature_value = entry.get("feature")
        if feature_present and not _needs_feature_entry(feature_value, sge_feature_path):
            entry.pop("feature", None)
            feature_present = False
        if not feature_present and legacy_feature is not None:
            entry["feature"] = legacy_feature
        train_params.append(entry)
        model_dict[(entry["model_type"], entry["model_id"])] = entry

    # 2. Model data structure
    # Process LDA model data from CLI (may extend or validate existing models)
    if args.lda_model is not None:
        for model in args.lda_model:
            model_type, model_path, model_id, train_width, n_factor, cmap = model.split(',')
            key = (model_type, model_id)
            existing = model_dict.get(key)
            if existing is None:
                model_entry = {
                    "model_type": model_type,
                    "model_id": model_id,
                    "model_path": model_path,
                    "train_width": int(train_width),
                    "n_factor": int(n_factor),
                    "cmap": cmap,
                    "decode_params": []
                }
                if _needs_feature_entry(args.in_feature_ficture, sge_feature_path):
                    model_entry["feature"] = args.in_feature_ficture
                model_dict[key] = model_entry
                train_params.append(model_entry)
            else:
                # append value to existing model and override fields if necessary
                msg_id = f"model ({model_type}, {model_id})"
                reconcile_field(
                    existing, "model_path", model_path,
                    type="override",
                    msg_id=msg_id,
                    normalize=lambda p: os.path.abspath(p)
                )
                reconcile_field(existing, "train_width", int(train_width), type="override", msg_id=msg_id)
                reconcile_field(existing, "n_factor", int(n_factor), type="override", msg_id=msg_id)
                reconcile_field(existing, "cmap", cmap, type="override", msg_id=msg_id)
                if _needs_feature_entry(args.in_feature_ficture, sge_feature_path):
                    reconcile_field(
                        existing, "feature", args.in_feature_ficture,
                        type="override",
                        msg_id=msg_id,
                        normalize=lambda p: os.path.abspath(p) if p is not None else None
                    )
                elif "feature" in existing and not _needs_feature_entry(existing.get("feature"), sge_feature_path):
                    existing.pop("feature", None)
    # 3. Process decode data    
    # Ensure decode_params list exists on all models
    for m in train_params:
        m.setdefault("decode_params", [])
    
    if args.decode is not None:
        for dec in args.decode:
            model_type, model_id, decode_id, fit_width, anchor_res = dec.split(',')
            key = (model_type, model_id)
            
            model_entry = model_dict.get(key)
            assert model_entry is not None, f"No matching model for decode entry: {dec}"
            
            decode_entry = {
                "decode_id": decode_id,
                "fit_width": int(fit_width),
                "anchor_res": int(anchor_res),
            }
            existing_dec = next((d for d in model_entry["decode_params"] if d["decode_id"] == decode_id), None)
            if existing_dec is None:
                model_entry["decode_params"].append(decode_entry)
            else:
                msg_id = f"model ({model_type}, {model_id}) - decode_id {decode_id}"
                reconcile_field(existing_dec, "fit_width", int(fit_width), type="override", msg_id=msg_id)
                reconcile_field(existing_dec, "anchor_res", int(anchor_res), type="override", msg_id=msg_id)

    # 4. Process UMAP data
    if args.umap is not None:
        for umap in args.umap:
            model_type, model_id, umap_tsv, umap_png, umap_single_factor_png = umap.split(',')
            key = (model_type, model_id)
            model_entry = model_dict.get(key)
            assert model_entry is not None, f"No matching model for UMAP entry: {umap}"
            
            umap_entry = {
                "tsv": umap_tsv,
                "png": umap_png,
                "ind_png": umap_single_factor_png
            }

            if "umap" not in model_entry:
                model_entry["umap"] = umap_entry
            else:
                existing_umap = model_entry["umap"]
                msg_id = f"model ({model_type}, {model_id}) - UMAP"
                reconcile_field(existing_umap, "tsv", umap_tsv, type="override", msg_id=msg_id)
                reconcile_field(existing_umap, "png", umap_png, type="override", msg_id=msg_id)
                reconcile_field(existing_umap, "ind_png", umap_single_factor_png, type="override", msg_id=msg_id)

    # 5. Construct final JSON data
    json_data = {
        "in_sge": sge_data,
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
