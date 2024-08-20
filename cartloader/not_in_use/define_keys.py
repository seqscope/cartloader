import argparse, inspect, os, sys, gzip

def define_keys(_args):
    """
    Generate sort keys for given columns and flags.

    :param input: Path to the input file (can be .gz or plain text)
    :param columns: List of columns with flags, in the format "column,flag"
    :param sep: Separator used in the input file (default is tab)
    :return: Sort keys joined by a space
    """
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Generate sort keys for specified columns.")
    parser.add_argument("--input", help="Path to the input file (can be .gz or plain text)")
    parser.add_argument("--sort", action="store_true", default=False, help="Generate sort keys")
    parser.add_argument("--tabix", action="store_true", default=False, help="Generate tabix keys")
    parser.add_argument("--columns", nargs="*", required=True, help="If sort, define it by columns with flags in the format 'column,flag'. If tabix, define it with column.")
    parser.add_argument("--sep", default="\t", help="Separator used in the input file (default is tab)")
    parser.add_argument("--option", choices=["return","print"],default="print", help="Return or print the keys")
    args = parser.parse_args(_args)

    # must define either sort or tabix
    assert args.sort + args.tabix == 1, "Must define one of --sort or --tabix"

    # Open the input file and read the header
    if args.input.endswith(".gz"):
        with gzip.open(args.input, "rt") as f:
            header = f.readline().strip()
    else:
        with open(args.input, "r") as f:
            header = f.readline().strip()

    headers = header.split(args.sep)
    
    if args.sort:
        # Parse columns and flags
        sort_keys = []
        for col_flag in args.columns:
            column, flag = col_flag.split(",")
            assert column in headers, f"Column {column} is not found in the input file"
            assert flag in ["n", "g"], f"Flag {flag} is not valid. Use 'n' for numerical or 'g' for general sorting"
            col_idx = headers.index(column) + 1
            sort_key = f"-k{col_idx},{col_idx}{flag}"
            sort_keys.append(sort_key)
        key_param= " ".join(sort_keys)
    elif args.tabix:
        len(args.columns) == 1, "Must define one column for tabix"
        column = args.columns[0]
        assert column in headers, f"Column {column} is not found in the input file"
        col_idx = headers.index(column) + 1
        key_param = f"-b{col_idx} -e{col_idx}"
    if args.option == "return":
        return key_param
    else:
        print (key_param)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
