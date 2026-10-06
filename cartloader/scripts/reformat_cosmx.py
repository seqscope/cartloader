import sys, os, argparse, gzip, inspect, datetime

from cartloader.utils.feature_filter import FeatureFilter


def flexopen(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode, encoding='utf-8' if 't' in mode else None)
    return open(filename, mode, encoding='utf-8' if 't' in mode else None)


def custom_log(message):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", file=sys.stderr)


def unquote(s):
    return s.replace('"', '').replace("'", "")


def parse_arguments(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Reformat raw CosMx SMI exports (transcript / metadata / polygon CSVs) into the "
                    "generic transcript TSV + cell metadata (xy) + polygon (boundaries) CSVs consumed "
                    "by FICTURE and cartload. Global pixel coordinates are converted to microns.")
    parser.add_argument('--tx', type=str, required=True, help='Transcript file (CosMx *_tx_file.csv[.gz])')
    parser.add_argument('--poly', type=str, help='Polygon file containing cell boundaries (CosMx *-polygons.csv[.gz])')
    parser.add_argument('--meta', type=str, help='Cell metadata file (CosMx *_metadata_file.csv[.gz])')
    parser.add_argument('--out', type=str, required=True, help='Output prefix for reformatted files')
    parser.add_argument('--tx-col-fov', type=str, default='fov', help='FOV column in the transcript file')
    parser.add_argument('--tx-col-cell', type=str, default='cell_ID', help='Cell-id column in the transcript file')
    parser.add_argument('--tx-col-gene', type=str, default='target', help='Gene/target column in the transcript file')
    parser.add_argument('--tx-col-x', type=str, default='x_global_px', help='Global x (px) column in the transcript file')
    parser.add_argument('--tx-col-y', type=str, default='y_global_px', help='Global y (px) column in the transcript file')
    parser.add_argument('--tx-col-z', type=str, default='z', help='Z column in the transcript file')
    parser.add_argument('--meta-col-fov', type=str, default='fov', help='FOV column in the metadata file')
    parser.add_argument('--meta-col-cell', type=str, default='cell_ID', help='Cell-id column in the metadata file')
    parser.add_argument('--meta-col-x', type=str, default='CenterX_global_px', help='Global x (px) column in the metadata file')
    parser.add_argument('--meta-col-y', type=str, default='CenterY_global_px', help='Global y (px) column in the metadata file')
    parser.add_argument('--poly-col-fov', type=str, default='fov', help='FOV column in the polygon file')
    parser.add_argument('--poly-col-cell', type=str, default='cellID', help='Cell-id column in the polygon file')
    parser.add_argument('--poly-col-x', type=str, default='x_global_px', help='Global x (px) column in the polygon file')
    parser.add_argument('--poly-col-y', type=str, default='y_global_px', help='Global y (px) column in the polygon file')
    parser.add_argument('--um-per-px', type=float, default=0.12028, help='Microns per pixel scaling factor (default: 0.12028, i.e., 1/8.3)')
    parser.add_argument('--out-suffix-tx', type=str, default='.transcripts.tsv.gz', help='Suffix for the output transcript TSV (default: .transcripts.tsv.gz)')
    parser.add_argument('--out-suffix-poly', type=str, default='.polygons.csv.gz', help='Suffix for the output polygon CSV (default: .polygons.csv.gz)')
    parser.add_argument('--out-suffix-meta', type=str, default='.metadata.csv.gz', help='Suffix for the output metadata CSV (default: .metadata.csv.gz)')
    parser.add_argument('--offset-x', type=float, default=0.0, help='Global x offset to add to all x coordinates (in microns)')
    parser.add_argument('--offset-y', type=float, default=0.0, help='Global y offset to add to all y coordinates (in microns)')
    # Feature filtering happens as the transcript TSV is written, so a dropped feature is
    # absent from every downstream product of this ingest.
    ftr = parser.add_argument_group("Feature Filtering Parameters")
    FeatureFilter.add_arguments(ftr, "the output transcript TSV")
    return parser.parse_args(_args)


def reformat_cosmx(_args):
    args = parse_arguments(_args)
    ftr_filter = FeatureFilter.from_args(args)

    custom_log(f"Processing transcript file: {args.tx}")
    with flexopen(args.tx, 'rt') as rf, flexopen(f"{args.out}{args.out_suffix_tx}", 'wt') as wf:
        hdrs = [unquote(h) for h in rf.readline().strip().split(',')]
        icol_fov = hdrs.index(args.tx_col_fov)
        icol_cell = hdrs.index(args.tx_col_cell)
        icol_gene = hdrs.index(args.tx_col_gene)
        icol_x = hdrs.index(args.tx_col_x)
        icol_y = hdrs.index(args.tx_col_y)
        icol_z = hdrs.index(args.tx_col_z)
        wf.write("X\tY\tgene\tcount\tcell_id\tZ\n")
        nlines = 0
        for line in rf:
            toks = line.strip().split(',')
            if toks[icol_cell] == '0':
                cell_id = "UNASSIGNED"
            else:
                cell_id = f"{toks[icol_fov]}_{toks[icol_cell]}"
            gene = toks[icol_gene].replace(' ', '-').replace('"', '').replace("'", '')
            if gene.startswith('System'):
                continue
            if not ftr_filter.keep(gene):
                continue
            x_um = float(toks[icol_x]) * args.um_per_px + args.offset_x
            y_um = float(toks[icol_y]) * args.um_per_px + args.offset_y
            z = toks[icol_z]
            wf.write(f"{x_um:.3f}\t{y_um:.3f}\t{gene}\t1\t{cell_id}\t{z}\n")
            nlines += 1
            if nlines % 1000000 == 0:
                custom_log(f"  Processed {nlines} transcript lines...")
    custom_log(f"Wrote transcripts to: {args.out}{args.out_suffix_tx}")
    if ftr_filter.active:
        custom_log(f"  Dropped {ftr_filter.n_dropped} transcript rows by the feature filters")

    if args.poly is not None:
        custom_log(f"Processing polygon file: {args.poly}")
        with flexopen(args.poly, 'rt') as rf, flexopen(f"{args.out}{args.out_suffix_poly}", 'wt') as wf:
            hdrs = [unquote(h) for h in rf.readline().strip().split(',')]
            icol_fov = hdrs.index(args.poly_col_fov)
            icol_cell = hdrs.index(args.poly_col_cell)
            icol_x = hdrs.index(args.poly_col_x)
            icol_y = hdrs.index(args.poly_col_y)
            wf.write("cell_id,vertex_x,vertex_y\n")
            for line in rf:
                toks = line.strip().split(',')
                if toks[icol_cell] == '0':
                    raise ValueError("Polygon file contains cell_ID of 0, which is invalid")
                cell_id = f"{toks[icol_fov]}_{toks[icol_cell]}"
                x_um = float(toks[icol_x]) * args.um_per_px + args.offset_x
                y_um = float(toks[icol_y]) * args.um_per_px + args.offset_y
                wf.write(f"{cell_id},{x_um:.3f},{y_um:.3f}\n")
        custom_log(f"Wrote polygons to: {args.out}{args.out_suffix_poly}")

    if args.meta is not None:
        custom_log(f"Processing metadata file: {args.meta}")
        with flexopen(args.meta, 'rt') as rf, flexopen(f"{args.out}{args.out_suffix_meta}", 'wt') as wf:
            hdrs = [unquote(h) for h in rf.readline().strip().split(',')]
            icol_cell = hdrs.index(args.meta_col_cell)
            icol_fov = hdrs.index(args.meta_col_fov)
            icol_x = hdrs.index(args.meta_col_x)
            icol_y = hdrs.index(args.meta_col_y)
            wf.write("cell_id,X,Y\n")
            for line in rf:
                toks = line.strip().split(',')
                toks[0] = unquote(toks[0])
                if toks[icol_cell] == '0':
                    raise ValueError("Metadata file contains cell_ID of 0, which is invalid")
                cell_id = f"{toks[icol_fov]}_{toks[icol_cell]}"
                x_um = float(toks[icol_x]) * args.um_per_px + args.offset_x
                y_um = float(toks[icol_y]) * args.um_per_px + args.offset_y
                wf.write(f"{cell_id},{x_um:.3f},{y_um:.3f}\n")
        custom_log(f"Wrote metadata to: {args.out}{args.out_suffix_meta}")

    custom_log("reformat_cosmx finished")


if __name__ == "__main__":
    reformat_cosmx(sys.argv[1:])
