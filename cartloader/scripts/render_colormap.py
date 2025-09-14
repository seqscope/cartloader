

import argparse
import os
import sys, inspect

def parse_arguments(_args):
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Render a color map legend image from TSV.")
    parser.add_argument("--in-tsv", required=True, help="Input TSV with columns [Name, R, G, B].")
    parser.add_argument("--out-png", required=True, help="Output PNG file path.")
    parser.add_argument("--swatch-size", type=float, default=0.6, help="Size of each color swatch.")
    parser.add_argument("--font-size", type=int, default=14, help="Label font size.")
    parser.add_argument("--dpi", type=int, default=150, help="Image DPI.")
    parser.add_argument("--max-per-row", type=int, default=12, help="Maximum labels per row before wrapping to next line (<=0 for single row).")
    parser.add_argument("--label-col", default="Name", help="Label column name (default: Name).")
    parser.add_argument("--r-col", default="R", help="Red column name (default: R).")
    parser.add_argument("--g-col", default="G", help="Green column name (default: G).")
    parser.add_argument("--b-col", default="B", help="Blue column name (default: B).")
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)


def _import_libs():
    try:
        import pandas as pd  # type: ignore
        import matplotlib.pyplot as plt  # type: ignore
        from matplotlib.patches import Rectangle  # type: ignore
        return pd, plt, Rectangle
    except Exception as e:
        sys.stderr.write(
            "Missing optional dependencies. Please install pandas and matplotlib.\n"
            f"Import error: {e}\n"
        )
        sys.exit(1)


def _normalize_rgb(r, g, b):
    """Return tuple in 0-1 range. Accepts 0-1 floats or 0-255 integers/floats."""
    # If any channel > 1.0, assume 0-255 scale
    if max(r, g, b) > 1.0:
        return (float(r) / 255.0, float(g) / 255.0, float(b) / 255.0)
    return (float(r), float(g), float(b))


def render_colormap(_args):
    """
    Render a color map legend image from a TSV file.

    Input TSV format (default columns):
        Name\tColor_index\tR\tG\tB
    Where R, G, B are either 0-1 floats or 0-255 integers.

    Example:
        cartloader render_colormap \
        --in-tsv docs/tabs/colormap/seqscope_starter.t18-f12-rgb.tsv \
        --out-png docs/images/starter_vignettes/seqscope.t18-f12-rgb.png

    This utility is intentionally simple and standalone for generating legend images
    used in documentation and examples. It does not alter any pipeline outputs.
    """
    args=parse_arguments(_args)
    pd, plt, Rectangle = _import_libs()

    df = pd.read_csv(args.in_tsv, sep="\t")
    required = [args.label_col, args.r_col, args.g_col, args.b_col]
    for col in required:
        if col not in df.columns:
            raise ValueError(
                f"Column '{col}' not found in TSV. Available columns: {list(df.columns)}"
            )

    n = len(df)
    if n == 0:
        raise ValueError("No rows found in TSV; nothing to render.")

    # Layout with wrapping
    cols = n if args.max_per_row <= 0 else min(n, args.max_per_row)
    rows = (n + cols - 1) // cols

    h_spacing = args.swatch_size + 0.4  # horizontal spacing between swatches
    v_spacing = args.swatch_size + 0.5  # vertical spacing between rows
    label_offset = 0.15                 # label offset from swatch
    fig_width = max(2.0, h_spacing * cols)
    fig_height = max(args.swatch_size + 0.3, v_spacing * rows)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.axis("off")

    for i, row in df.iterrows():
        try:
            color = _normalize_rgb(row[args.r_col], row[args.g_col], row[args.b_col])
        except Exception as e:
            raise ValueError(f"Invalid RGB values in row {i}: {e}")

        r = i // cols
        c = i % cols
        # draw from top to bottom so first row appears at top
        x = c * h_spacing
        y = (rows - 1 - r) * v_spacing
        # Draw color swatch
        rect = Rectangle((x, y), args.swatch_size, args.swatch_size,
                         color=color, ec='white', linewidth=1.5)
        ax.add_patch(rect)

        # Factor label to the right of the swatch
        label = str(row[args.label_col])
        try:
            # if numeric like 0.0 → 0
            label = str(int(float(label)))
        except Exception:
            pass
        ax.text(x + args.swatch_size + label_offset, y + args.swatch_size / 2.0,
                label, va='center', ha='left', fontsize=args.font_size, color='dimgray')

    ax.set_xlim(-0.2, h_spacing * cols)
    ax.set_ylim(-0.2, v_spacing * rows)
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out_png) or ".", exist_ok=True)
    fig.savefig(args.out_png, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
