import pandas as pd
import argparse
import sys, os, gzip, argparse, logging, warnings, shutil
import pandas as pd

def hex_to_rgb(hex_code):
    """Convert HEX color code to normalized RGB values."""
    hex_code = hex_code.lstrip('#')
    r, g, b = [int(hex_code[i:i+2], 16) / 255.0 for i in (0, 2, 4)]
    return r, g, b

def cmap_hex2rgb(_args):
    parser = argparse.ArgumentParser(description="Convert HEX codes to RGB in a TSV file.")
    parser.add_argument("--input", help="Path to the input file")
    parser.add_argument("--col-hex", default="Hex_Code", help="Name of the column containing HEX codes")
    parser.add_argument("--col-index", default="Color_Index", help="Name of the column containing color index values")
    parser.add_argument("--output", default="output.tsv", help="Path to save the output TSV file (default: output.tsv)")
    
    args = parser.parse_args(_args)
    
    """Convert a TSV file with HEX codes to RGB and save the result."""
    # Read the input TSV file
    if args.input.endswith('.tsv'):
        df = pd.read_csv(args.input, sep='\t')
    elif args.input.endswith('.csv'):
        df = pd.read_csv(args.input, sep=',')
    else:
        raise ValueError("Input file must be in TSV or CSV format.")
    
    # Convert HEX to RGB
    df['R'], df['G'], df['B'] = zip(*df[args.col_hex].apply(hex_to_rgb))

    # color index column
    if args.col_index != 'Color_index':
        df.rename(columns={args.col_index: 'Color_index'}, inplace=True)

    if df['Color_index'].min() != 0:
        df = df.sort_values(by='Color_index').reset_index(drop=True)
        df['Color_index'] = df.index
    
    # Name
    df.insert(0, 'Name', df['Color_index'])
    df = df[['Name','Color_index', 'R', 'G', 'B']]

    df.to_csv(args.output, sep='\t', index=False, float_format='%.5f')
    print(f"Converted data saved to {args.output}")

if __name__ == "__main__":
    # Get the path to the cartloader repository
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])

