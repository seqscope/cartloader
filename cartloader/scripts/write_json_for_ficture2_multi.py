import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json, copy
from collections import defaultdict, OrderedDict
# from cartloader.utils.minimake import minimake
# from cartloader.utils.utils import cmd_separator, scheck_app

from cartloader.utils.utils import write_json, reconcile_field

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="Write a JSON file to summarize the parameters.")
    parser.add_argument('--out-dir', required=True, type=str, help='Output directory')
    parser.add_argument('--out-json', type=str, default=None, help='Path to the output JSON file. Default: <out-dir>/ficture.params.json')
    parser.add_argument('--mode', type=str, default="write", choices=["write", "append"], help='Write mode for the output JSON. Default: write. If write, a new file will be created based on the arguments provided. If append, the new parameters will be merged into the existing JSON file if it exists.')
    parser.add_argument('--in-transcript', type=str, default=None, help='Path to the transcript file.')
    parser.add_argument('--in-feature', type=str, default=None, help='Path to the feature file.')
    parser.add_argument('--in-minmax', type=str, default=None, help='Path to the minmax file.')
    parser.add_argument('--in-feature-ficture', type=str, default=None, help='(Optional) If FICTURE used a different feature file than the in-feature file, specify the path to the feature file used for FICTURE analysis.')
    parser.add_argument('--lda-model', nargs='*', type=str, default=None, help='LDA Model information: <model_type>,<model_path>,<model_id>,<train_width>,<n_factor>,<cmap>')
    parser.add_argument('--decode', nargs='*', type=str, default=None, help='Projection information: <model_type>,<model_id>,<projection_id>,<fit_width>,<anchor_res>')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args) 

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

    if args.mode == "write" or not old_data:
        for key in sge_keys:
            path = new_sge[key]
            assert path is not None, f"Path not provided: --{key.replace('_', '-')}"
            assert os.path.exists(path), f"File not found: {path} ( --{key.replace('_', '-')})"
        return new_sge

    old_sge = old_data.get("in_sge", {})
    final_sge = {}

    for key in sge_keys:
        new_path = new_sge[key]
        old_path = old_sge.get(key)

        if new_path is None:
            assert old_path is not None, f"Path for --{key.replace('_', '-')} not provided and not present in existing JSON."
            assert os.path.exists(old_path), f"File not found: {old_path} (from existing JSON for --{key.replace('_', '-')})"
            final_sge[key] = old_path
        else:
            assert os.path.exists(new_path), ( f"File not found: {new_path} ( --{key.replace('_', '-')})")
            if old_path is not None:
                assert os.path.abspath(new_path) == os.path.abspath(old_path), f"Found inconsistent absolute path for '{key}' SGE between existing json ({old_path}) and the input arguments ({new_path})."
                final_sge[key] = old_path
            else:
                final_sge[key] = new_path
    return final_sge

def write_json_for_ficture2_multi(_args):
    args = parse_arguments(_args)
    if args.out_json is None:
        args.out_json = os.path.join(args.out_dir, "ficture.params.json")
    
    old_data = {}
    old_train_params = []
    old_feature_ficture = None
    model_dict = {}

    if args.mode == "append" and os.path.exists(args.out_json):
        with open(args.out_json, "r") as f:
            old_data = json.load(f)
        old_feature_ficture = old_data.get("in_feature_ficture")
        old_train_params = copy.deepcopy(old_data.get("train_params", []))

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
        entry.setdefault("decode_params", [])
        train_params.append(entry)
        model_dict[(entry["model_type"], entry["model_id"])] = entry

    if args.lda_model is not None:
        for model in args.lda_model:
            model_type, model_id, train_width, n_factor, cmap, model_path, fit_tsv_path, de_tsv_path, info_tsv_path = model.split(',')
            key = (model_type, model_id)
            existing = model_dict.get(key)
            if existing is None:
                model_entry = {
                    "model_type": model_type,
                    "model_id": model_id,
                    "train_width": int(train_width),
                    "n_factor": int(n_factor),
                    "cmap": cmap,
                    "model_path": model_path,
                    "fit_path": fit_tsv_path,
                    "de_path": de_tsv_path,
                    "info_path": info_tsv_path,
                    "decode_params": []
                }
                if _needs_feature_entry(args.in_feature_ficture, sge_feature_path):
                    model_entry["feature"] = args.in_feature_ficture
                model_dict[key] = model_entry
                train_params.append(model_entry)
            else:
                msg_id = f"model ({model_type}, {model_id})"
                reconcile_field(
                    existing, "train_width", int(train_width),
                    type="override",
                    msg_id=msg_id
                )
                reconcile_field(existing, "n_factor", int(n_factor), type="override", msg_id=msg_id)
                reconcile_field(existing, "cmap", cmap, type="override", msg_id=msg_id)
                reconcile_field(
                    existing, "model_path", model_path,
                    type="override",
                    msg_id=msg_id,
                    normalize=lambda p: os.path.abspath(p)
                )
                reconcile_field(existing, "fit_path", fit_tsv_path, type="override", msg_id=msg_id)
                reconcile_field(existing, "de_path", de_tsv_path, type="override", msg_id=msg_id)
                reconcile_field(existing, "info_path", info_tsv_path, type="override", msg_id=msg_id)
                if _needs_feature_entry(args.in_feature_ficture, sge_feature_path):
                    reconcile_field(
                        existing, "feature", args.in_feature_ficture,
                        type="override",
                        msg_id=msg_id,
                        normalize=lambda p: os.path.abspath(p) if p is not None else None
                    )
                elif "feature" in existing and not _needs_feature_entry(existing.get("feature"), sge_feature_path):
                    existing.pop("feature", None)

    # Process decode data
    if args.decode is not None:
        for dec in args.decode:
            model_type, model_id, decode_id, fit_width, anchor_res, decode_pixel_tsv, decode_pixel_png, decode_pseudobulk_tsv, decode_de_tsv, decode_info_tsv  = dec.split(',')
            decode_entry = {
                "decode_id": decode_id,
                "fit_width": int(fit_width),
                "anchor_res": int(anchor_res),
                "pixel_tsv_path": decode_pixel_tsv,
                "pixel_png_path": decode_pixel_png,
                "pseudobulk_tsv_path": decode_pseudobulk_tsv,
                "de_tsv_path": decode_de_tsv,
                "info_tsv_path": decode_info_tsv
            }
            model_entry = model_dict.get((model_type, model_id))
            assert model_entry is not None, f"No matching model for decode entry: {dec}"
            existing_dec = next((d for d in model_entry["decode_params"] if d["decode_id"] == decode_id), None)
            if existing_dec is None:
                model_entry["decode_params"].append(decode_entry)
            else:
                msg_id = f"model ({model_type}, {model_id}) - decode_id {decode_id}"
                reconcile_field(existing_dec, "fit_width", int(fit_width), type="override", msg_id=msg_id)
                reconcile_field(existing_dec, "anchor_res", int(anchor_res), type="override", msg_id=msg_id)
                reconcile_field(existing_dec, "pixel_tsv_path", decode_pixel_tsv, type="override", msg_id=msg_id)
                reconcile_field(existing_dec, "pixel_png_path", decode_pixel_png, type="override", msg_id=msg_id)
                reconcile_field(existing_dec, "pseudobulk_tsv_path", decode_pseudobulk_tsv, type="override", msg_id=msg_id)
                reconcile_field(existing_dec, "de_tsv_path", decode_de_tsv, type="override", msg_id=msg_id)
                reconcile_field(existing_dec, "info_tsv_path", decode_info_tsv, type="override", msg_id=msg_id)
    
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
