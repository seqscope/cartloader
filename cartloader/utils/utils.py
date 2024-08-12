import logging, os, shutil, sys, importlib, csv, shlex, subprocess

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