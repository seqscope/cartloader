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
    parser.add_argument('--model', nargs='*', type=str, default=None, help='Model information: <model_type>,<model_path>,<model_id>,<train_width>,<n_factor>')
    parser.add_argument('--projection', nargs='*', type=str, default=None, help='Projection information: <model_type>,<model_id>,<projection_id>,<fit_width>,<anchor_res>')
    parser.add_argument('--decode', nargs='*', type=str, default=None, help='Decode information: <model_type>,<model_id>,<projection_id>,<decode_id>,<radius>')
    
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args) 

def write_json(data, filename):
    """Write the given data to a JSON file."""
    with open(filename, 'w') as file:
        json.dump(data, file, indent=4, sort_keys=False)

def write_json_for_ficture(_args):
    args = parse_arguments(_args)
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
    if args.model is not None:
        for model in args.model:
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
