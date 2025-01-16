import logging, os, shutil, sys, importlib, csv, shlex, subprocess, json, yaml, re

def cmd_separator(cmds, info):
    """
    Append messages separating between commands
    """
    cmds.append(rf"$(info --------------------------------------------------------------)")
    cmds.append(rf"$(info {info})")
    cmds.append(rf"$(info --------------------------------------------------------------)")
    return cmds

def scheck_app(app_cmd):
    """
    Check if the specified application is available
    """
    #if not shutil.which(app_cmd.split(" ")[0]):
    if not shutil.which(app_cmd):
        logging.error(f"Cannot find {app_cmd}. Please make sure that the path to specify {app_cmd} is correct")
        sys.exit(1)

def get_func(name):
    """
    Get the function object among the runnable scripts based on the script name
    """
    #print(f"get_func({name}) was called")
    module = importlib.import_module(f"cartloader.scripts.{name}")
    return getattr(module,name)

def create_custom_logger(name, logfile=None, level=logging.INFO):
    """
    Create a custom logger object
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False  # Prevent log messages from being propagated to parent loggers
    
    if logger.hasHandlers():
        logger.handlers.clear()
    
    log_console_handler = logging.StreamHandler()
    log_console_format = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    log_console_handler.setFormatter(log_console_format)
    logger.addHandler(log_console_handler)
    
    ## create a directory for the logger file if needed
    if logfile is not None:
        output_dir = os.path.dirname(logfile)
        if not os.path.exists(output_dir):
            logger.info(f"Creating directory {output_dir} for storing output")
            os.makedirs(output_dir)
        
        log_file_handler = logging.FileHandler(logfile)
        log_file_handler.setLevel(logging.INFO)
        log_file_format = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        log_file_handler.setFormatter(log_file_format)
        logger.addHandler(log_file_handler)
    
    return logger

# copied from NovaScope
def find_major_axis(filename, format):
    # purpose: find the longer axis of the image
    # (1) detect from the "{uid}.coordinate_minmax.tsv"  
    if format == "row":
        data = {}
        with open(filename, 'r') as file:
            for line in file:
                key, value = line.split()
                data[key] = float(value)
        if not all(k in data for k in ['xmin', 'xmax', 'ymin', 'ymax']):
            raise ValueError("Missing one or more required keys (xmin, xmax, ymin, ymax)")
        xmin = data['xmin']
        xmax = data['xmax']
        ymin = data['ymin']
        ymax = data['ymax']
    # (2) detect from the barcodes.minmax.tsv
    elif format == "col":
        with open(filename, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            try:
                row = next(reader)  # Read the first row
            except StopIteration:
                raise ValueError("File is empty")
            try:
                next(reader)
                raise ValueError("Error: More than one row of data found.")
            except StopIteration:
                pass  
            xmin = int(row['xmin'])
            xmax = int(row['xmax'])
            ymin = int(row['ymin'])
            ymax = int(row['ymax'])
    deltaX = xmax - xmin
    deltaY = ymax - ymin
    if deltaX > deltaY:
        return "X"
    else:
        return "Y"


def log_dataframe(df, log_message="DataFrame Info:", indentation=""): 
    ## purpose:format the log messages so that each column value occupies a fixed width

    # Calculate column widths
    col_widths = {col: max(df[col].astype(str).apply(len).max(), len(col)) for col in df.columns}

    # Prepare the header string with column names aligned
    header = ' | '.join([col.ljust(col_widths[col]) for col in df.columns])

    # Log the header
    logging.info(f"{log_message}")
    logging.info(indentation+header)
    logging.info(indentation+"-" * len(header))  # Divider line

    # Iterate over DataFrame rows and log each, maintaining alignment
    for _, row in df.iterrows():
        row_str = ' | '.join([str(row[col]).ljust(col_widths[col]) for col in df.columns])
        logging.info(indentation+row_str)


def read_minmax(filename, format):
    # purpose: find the longer axis of the image
    # (1) detect from the "{uid}.coordinate_minmax.tsv"  
    if format == "row":
        data = {}
        with open(filename, 'r') as file:
            for line in file:
                key, value = line.split()
                data[key] = float(value)
        if not all(k in data for k in ['xmin', 'xmax', 'ymin', 'ymax']):
            raise ValueError("Missing one or more required keys (xmin, xmax, ymin, ymax)")
        xmin = data['xmin']
        xmax = data['xmax']
        ymin = data['ymin']
        ymax = data['ymax']
    # (2) detect from the barcodes.minmax.tsv
    elif format == "col":
        with open(filename, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            try:
                row = next(reader)  # Read the first row
            except StopIteration:
                raise ValueError("File is empty")
            try:
                next(reader)
                raise ValueError("Error: More than one row of data found.")
            except StopIteration:
                pass  
            xmin = int(row['xmin'])
            xmax = int(row['xmax'])
            ymin = int(row['ymin'])
            ymax = int(row['ymax'])
    return {
        "xmin": xmin,
        "xmax": xmax,
        "ymin": ymin,
        "ymax": ymax
    }
    
# def add_param_to_cmd(cmd, args, aux_argset):
#     aux_args = {k: v for k, v in vars(args).items() if k in aux_argset}
#     for arg, value in aux_args.items():
#         if value or isinstance(value, bool):
#             arg_name = arg.replace('_', '-')
#             if isinstance(value, bool) and value:
#                 cmd += f" --{arg_name}"
#             elif isinstance(value, list):
#                 cmd += f" --{arg_name} {' '.join(value)}"
#             elif not isinstance(value, bool):
#                 if "regex" in arg_name:
#                     cmd += f" --{arg_name} '{value}'"
#                 else:
#                     cmd += f" --{arg_name} {value}"
#     return cmd

def add_param_to_cmd(cmd, args, aux_argset, underscore2dash=True):
    aux_args = {k: v for k, v in vars(args).items() if k in aux_argset}
    for arg, value in aux_args.items():
        if value or isinstance(value, bool):
            if underscore2dash:
                arg_name = arg.replace('_', '-')
            else:
                arg_name = arg
            if isinstance(value, bool) and value:
                cmd += f" --{arg_name}"
            elif isinstance(value, list):
                cmd += f" --{arg_name} {' '.join(map(str, value))}"
            elif not isinstance(value, bool):
                # Ensure regex patterns are properly quoted
                if "regex" in arg_name:
                    quoted_value = shlex.quote(str(value))
                    cmd += f" --{arg_name} {quoted_value}"
                else:
                    cmd += f" --{arg_name} {value}"
    return cmd

## code suggested by ChatGPT
def load_file_to_dict(file_path, file_type=None):
    """
    Load a JSON or YAML file into a dictionary.

    Parameters:
    file_path (str): Path to the JSON or YAML file.
    file_type (str, optional): The type of the file ('json' or 'yaml'). If None, the type is inferred from the file extension.

    Returns:
    dict: Dictionary containing the file's data.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    if file_type is None:
        _, file_extension = os.path.splitext(file_path)
        file_type = file_extension.lower()[1:]  # Strip the dot and use the extension

    if file_type == 'json':
        with open(file_path, 'r') as file:
            return json.load(file)
    elif file_type in ['yaml', 'yml']:
        with open(file_path, 'r') as file:
            return yaml.safe_load(file)
    else:
        raise ValueError("Unsupported file type. Please provide 'json' or 'yaml'/'yml' as file_type.")
    
def factor_id_to_name(factor_id):
    pattern = re.compile(r"(?:(t\d+)-)?(?:(f\d+)-)?(?:(p\d+)-)?(?:(a\d+)-)?(?:(r\d+))?")
    match = pattern.match(factor_id)
    if match and match.group(0) == factor_id:
        # Filter out None values and organize the result
        res = {key: value for key, value in zip(['t', 'f', 'p', 'a', 'r'], match.groups()) if value}
        if "t" in res and "f" in res:
            if "p" in res: ## projection
                if "a" in res:
                    if "r" in res:
                        return f"FICTURE {res['f']} factors - radius {res['r']}um / LDA {res['t']}-{res['p']}-{res['a']}um"
                    else:
                        return f"FICTURE {res['f']} factors / LDA {res['t']}-{res['p']}-{res['a']}um"
                else:
                    return f"FICTURE {res['f']} factors / LDA {res['t']}-{res['p']}um"
            else:
                return f"FICTURE {res['f']} factors / LDA {res['t']}um"
        else: ## not parseable
            return factor_id
    else:
        return factor_id

## transform FICTURE parameters to FACTOR assets (new standard)
def ficture_params_to_factor_assets(params, skip_raster=False):
    ## model_id
    ## proj_params -> proj_id
    ## proj_params -> decode_params -> decode_id
    suffix_factormap = "-factor-map.tsv"
    suffix_de = "-bulk-de.tsv"
    suffix_info = "-info.tsv"
    suffix_model = "-model-matrix.tsv.gz"
    suffix_post = "-posterior-counts.tsv.gz"
    suffix_rgb = "-rgb.tsv"
    suffix_hex_coarse = ".pmtiles"
    suffix_hex_fine = ".pmtiles"
    suffix_raster = "-pixel-raster.pmtiles"

    out_assets = []
    for param in params: ## train_params is a list of dictionaries
        if not "model_id" in param:
            raise ValueError(f"model_id is missing from FICTURE parameters")
        model_id = param["model_id"].replace("_", "-")

        len_proj_params = len(param["proj_params"])

        if len_proj_params == 0: ## train_param only
            out_asset = {
                "id": model_id,
                "name": factor_id_to_name(model_id),
                "model_id": model_id,
                "de": model_id + suffix_de,
                "info": model_id + suffix_info,
                "model": model_id + suffix_model,
                "post": model_id + suffix_post,
                "rgb": model_id + suffix_rgb,
                "pmtiles": {
                    "hex_coarse": model_id + suffix_hex_coarse
                }
            }
            if "factor_map" in param:
                out_asset["factor_map"] = model_id + suffix_factormap
            out_assets.append(out_asset)
        elif len_proj_params == 1:
            proj_param = param["proj_params"][0]
            if not "proj_id" in proj_param:
                raise ValueError(f"proj_id is missing from FICTURE parameter for {model_id}")
            proj_id = proj_param["proj_id"].replace("_", "-")

            len_decode_params = len(proj_param["decode_params"])

            if len_decode_params == 0:
                out_asset = {
                    "id": model_id,
                    "name": factor_id_to_name(model_id),
                    "model_id": model_id,
                    "proj_id": proj_id,
                    "de": model_id + suffix_de,
                    "info": model_id + suffix_info,
                    "model": model_id + suffix_model,
                    "post": model_id + suffix_post,
                    "rgb": model_id + suffix_rgb,
                    "pmtiles": {
                        "hex_coarse": model_id + suffix_hex_coarse,
                        "hex_fine": proj_id + suffix_hex_fine
                    }
                }
                if "factor_map" in param:
                    out_asset["factor_map"] = model_id + suffix_factormap
                out_assets.append(out_asset)
            elif len_decode_params == 1:
                decode_param = proj_param["decode_params"][0]
                if not "decode_id" in decode_param:
                    raise ValueError(f"decode_id is missing from FICTURE parameter for {model_id}-{proj_id}")
                decode_id = decode_param["decode_id"].replace("_", "-")

                out_asset = {
                    "id": model_id,
                    "name": factor_id_to_name(model_id),
                    "model_id": model_id,
                    "proj_id": proj_id,
                    "decode_id": decode_id,
                    "de": decode_id + suffix_de,
                    "info": decode_id + suffix_info,
                    "model": model_id + suffix_model,
                    "post": decode_id + suffix_post,
                    "rgb": model_id + suffix_rgb,
                    "pmtiles": {
                        "hex_coarse": model_id + suffix_hex_coarse,
                        "hex_fine": proj_id + suffix_hex_fine,
                        **({"raster": decode_id + suffix_raster} if not skip_raster else {})
                    }
                }
                if "factor_map" in param:
                    out_asset["factor_map"] = model_id + suffix_factormap
                out_assets.append(out_asset)
            else:
                for decode_param in proj_param["decode_params"]:
                    if not "decode_id" in decode_param:
                        raise ValueError(f"decode_id is missing from FICTURE parameter for {model_id}-{proj_id}")
                    decode_id = decode_param["decode_id"].replace("_", "-")

                    out_asset = {
                        "id": decode_id,
                        "name": factor_id_to_name(model_id),
                        "model_id": model_id,
                        "proj_id": proj_id,
                        "decode_id": decode_id,
                        "de": decode_id + suffix_de,
                        "info": decode_id + suffix_info,
                        "model": decode_id + suffix_model,
                        "post": decode_id + suffix_post,
                        "rgb": model_id + suffix_rgb,
                        "pmtiles": {
                            "hex_coarse": model_id + suffix_hex_coarse,
                            "hex_fine": proj_id + suffix_hex_fine,
                            **({"raster": decode_id + suffix_raster} if not skip_raster else {})
                        }
                    }
                    if "factor_map" in param:
                        out_asset["factor_map"] = model_id + suffix_factormap
                    out_assets.append(out_asset)
        else: ## multiple proj_params
            for proj_param in param["proj_params"]:
                if not "proj_id" in proj_param:
                    raise ValueError(f"proj_id is missing from FICTURE parameter for {model_id}")
                proj_id = proj_param["proj_id"].replace("_", "-")

                len_decode_params = len(proj_param["decode_params"])

                if len_decode_params == 0: ## train and projection only
                    out_asset = {
                        "id": model_id,
                        "name": factor_id_to_name(model_id),
                        "model_id": model_id,
                        "proj_id": proj_id,
                        "de": model_id + suffix_de,
                        "info": model_id + suffix_info,
                        "model": model_id + suffix_model,
                        "post": model_id + suffix_post,
                        "rgb": model_id + suffix_rgb,
                        "pmtiles": {
                            "hex_coarse": model_id + suffix_hex_coarse,
                            "hex_fine": proj_id + suffix_hex_fine
                        }
                    }
                    if "factor_map" in param:
                        out_asset["factor_map"] = model_id + suffix_factormap
                    out_assets.append(out_asset)
                else:
                    for decode_param in proj_param["decode_params"]:
                        if not "decode_id" in decode_param:
                            raise ValueError(f"decode_id is missing from FICTURE parameter for {model_id}-{proj_id}")
                        decode_id = decode_param["decode_id"].replace("_", "-")

                        out_asset = {
                            "id": decode_id,
                            "name": factor_id_to_name(model_id),
                            "model_id": model_id,
                            "proj_id": proj_id,
                            "decode_id": decode_id,
                            "de": decode_id + suffix_de,
                            "info": decode_id + suffix_info,
                            "model": decode_id + suffix_model,
                            "post": decode_id + suffix_post,
                            "rgb": model_id + suffix_rgb,
                            "pmtiles": {
                                "hex_coarse": model_id + suffix_hex_coarse,
                                "hex_fine": proj_id + suffix_hex_fine,
                                **({"raster": decode_id + suffix_raster} if not skip_raster else {})
                            }
                        }
                        if "factor_map" in param:
                            out_asset["factor_map"] = model_id + suffix_factormap
                        out_assets.append(out_asset)
    return out_assets

## code suggested by ChatGPT
def write_dict_to_file(data, file_path, file_type=None):
    """
    Write a dictionary to a JSON or YAML file.

    Parameters:
    data (dict): The dictionary to write to the file.
    file_path (str): Path to the output file.
    file_type (str, optional): The type of the file ('json' or 'yaml'). If None, the type is inferred from the file extension.

    Raises:
    ValueError: If the file type is unsupported.
    """
    if file_type is None:
        _, file_extension = os.path.splitext(file_path)
        file_type = file_extension.lower()[1:]  # Strip the dot and use the extension

    if file_type == 'json':
        with open(file_path, 'w') as file:
            json.dump(data, file, indent=4)
    elif file_type in ['yaml', 'yml']:
        with open(file_path, 'w') as file:
            yaml.safe_dump(data, file, default_flow_style=False)
    else:
        raise ValueError("Unsupported file type. Please provide 'json' or 'yaml'/'yml' as file_type.")

# def run_bash_command(command):
#     try:
#         result = subprocess.run(command, 
#                   shell=True, 
#                   stdout=subprocess.PIPE, 
#                   stderr=subprocess.PIPE, 
#                   text=True, 
#                   check=True)
#         return result.stdout
#     except subprocess.CalledProcessError as e:
#         print(f"Command failed with error:\n\t{command}\n\t{e.stderr}\n")
#         raise

def run_command(command, use_bash=False):
    executable_shell = "/bin/bash" if use_bash else "/bin/sh"
    try:
        result = subprocess.run(
            command, 
            shell=True, 
            executable=executable_shell,  # Choose the shell based on the use_bash flag
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True, 
            check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error:\n{e.stderr}")
        raise

import os

def create_symlink(A, B):
    # Purpose: Create a soft link from A to B

    # Step 0: Check if A and B are the same
    if A == B:
        print("A is the same as B. No need to create a soft link.")
        return

    # Step 1: Check if A exists
    if not os.path.exists(A):
        raise FileNotFoundError(f"The source path {A} does not exist.")
    
    # Step 2: Check if B exist
    if os.path.islink(B) and not os.path.exists(os.path.realpath(B)): # If B is a broken link, remove it
        os.remove(B)
    elif os.path.exists(B):                                           # If B is a valid file or valid link.
        real_B = os.path.realpath(B)
        real_A = os.path.realpath(A)
        
        if real_A == real_B:
            print("The source and destination files are the same. No need to create a soft link.")
            return
        else:
            if os.path.isfile(B):
                raise FileExistsError(f"The destination path {B} exists as a regular file. Please take action before creating a soft link.")
            elif os.path.islink(B):
                os.remove(B) 
        
    # Step 3: Create the soft link
    os.symlink(A, B)
    print(f"Soft link created from {A} to {B}.")