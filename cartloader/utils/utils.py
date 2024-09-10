import logging, os, shutil, sys, importlib, csv, shlex, subprocess, json, yaml

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
    if not shutil.which(app_cmd.split(" ")[0]):
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

def add_param_to_cmd(cmd, args, aux_argset):
    aux_args = {k: v for k, v in vars(args).items() if k in aux_argset}
    for arg, value in aux_args.items():
        if value or isinstance(value, bool):
            arg_name = arg.replace('_', '-')
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