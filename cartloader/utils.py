import logging, os, shutil, sys

def cmd_separator(cmds, info):
    cmds.append(rf"$(info --------------------------------------------------------------)")
    cmds.append(rf"$(info {info})")
    cmds.append(rf"$(info --------------------------------------------------------------)")
    return cmds


def scheck_app(app_cmd):
    if not shutil.which(app_cmd.split(" ")[0]):
        logging.error(f"Cannot find {app_cmd}. Please make sure that the path to --gzip is correct")
        sys.exit(1)