import os, sys, logging, importlib
import logging

def main():
    module_list = []
    scripts_dir = os.path.join(os.path.dirname(__file__), 'scripts')
    for f in os.listdir(scripts_dir):
        if f.endswith(".py") and f != "__init__.py":
            module_list.append(os.path.splitext(f)[0])

    module_map = {m:m for m in module_list}

    if len(sys.argv) < 2:
        print("Usage: cartloader <command> <args>, cartloader <command> -h to see arguments for each command")
        print("Available commands:\n\t"+"\n\t".join(sorted(list(module_map.keys()) )))
        return

    # Accept dashes in a command name (e.g. "filter-molecules") as aliases for the
    # underscore module name; Python modules cannot contain dashes, so this never
    # collides with a real command.
    function_name = sys.argv[1].replace("-", "_")
    if function_name not in module_map:
        print("Unknown command: "+sys.argv[1])
        print("Available commands:\n\t"+"\n\t".join(sorted(list(module_map.keys()) )))
        return

    module_name = "cartloader.scripts." + module_map[function_name]
    module = importlib.import_module(module_name)
    function = getattr(module, function_name)

    function(sys.argv[2:])
