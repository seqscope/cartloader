import sys, os, argparse, gzip, inspect, datetime, subprocess, shutil


def flexopen(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode, encoding='utf-8' if 't' in mode else None)
    return open(filename, mode, encoding='utf-8' if 't' in mode else None)


def custom_log(message):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", file=sys.stderr)


def parse_arguments(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Convert a cell-by-gene matrix CSV (e.g. MERSCOPE cell_by_gene.csv) into "
                    "10x-style MEX files (barcodes/features/matrix). Each row is one cell "
                    "(first column = cell id / barcode), each remaining column a feature. "
                    "Feature columns whose name starts with the --blank-prefix (control probes) "
                    "are dropped.")
    parser.add_argument('--csv', type=str, required=True, help='Input cell-by-gene CSV (may be .gz)')
    parser.add_argument('--out-dir', type=str, required=True, help='Output MEX directory')
    parser.add_argument('--bcd', type=str, default='barcodes.tsv.gz', help='Barcode file name (default: barcodes.tsv.gz)')
    parser.add_argument('--ftr', type=str, default='features.tsv.gz', help='Feature file name (default: features.tsv.gz)')
    parser.add_argument('--mtx', type=str, default='matrix.mtx.gz', help='Matrix file name (default: matrix.mtx.gz)')
    parser.add_argument('--blank-prefix', type=str, default='Blank', help='Drop feature columns whose name starts with this prefix (default: Blank)')
    parser.add_argument('--gzip', type=str, default='gzip', help='Compression binary for the matrix file (default: gzip; use pigz for parallel)')
    parser.add_argument('--threads', type=int, default=1, help='Threads passed to pigz (only used when --gzip is pigz; default: 1)')
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)


def convert_cellxgene(_args):
    args = parse_arguments(_args)

    os.makedirs(args.out_dir, exist_ok=True)

    # Resolve the compression command for the matrix body. pigz honors -p <threads>;
    # a plain gzip does not, so only append -p when the binary looks like pigz.
    gzip_base = os.path.basename(args.gzip.split()[0])
    if shutil.which(gzip_base) is None:
        sys.exit(f"ERROR: compression binary '{gzip_base}' not found in PATH (from --gzip {args.gzip}).")
    gzip_cmd = f"{args.gzip} -p {args.threads} -c" if gzip_base == "pigz" else f"{args.gzip} -c"

    custom_log(f"Reading cell-by-gene CSV file {args.csv}....")
    with flexopen(args.csv, 'rt') as rf:
        hdrs = rf.readline().rstrip("\r\n").split(",")
        nhdrs = len(hdrs)
        # Feature columns: everything after the first (cell-id) column, excluding
        # control-probe columns (name starts with --blank-prefix).
        icols = [i for i in range(1, nhdrs) if not hdrs[i].startswith(args.blank_prefix)]
        nftrs = len(icols)
        if nftrs == 0:
            sys.exit(f"ERROR: no feature columns found in {args.csv} (after dropping "
                     f"'{args.blank_prefix}*' columns). Header had {nhdrs} columns.")

        with flexopen(os.path.join(args.out_dir, args.ftr), 'wt') as ff:
            for i in icols:
                gene_name = hdrs[i]
                ff.write(f"{gene_name}\t{gene_name}\tGene Expression\n")

        nbcds = 0        # number of cells kept (>=1 non-zero feature)
        total_nnz = 0    # total non-zero matrix entries (the MatrixMarket header count)
        nlines = 0
        mtx_path = os.path.join(args.out_dir, args.mtx)
        body_path = f"{mtx_path}.body.tmp"
        hdr_path = f"{mtx_path}.hdr.tmp"
        with flexopen(os.path.join(args.out_dir, args.bcd), 'wt') as bf, \
             open(body_path, 'wt') as mf:
            for line in rf:
                toks = line.rstrip("\r\n").split(",")
                if len(toks) != nhdrs:
                    raise ValueError(f"Expected {nhdrs} columns but found {len(toks)} at data line {nlines + 1}")
                bcd = toks[0]
                nbcds += 1
                nnz = 0
                for j, i in enumerate(icols):
                    val = int(float(toks[i]))
                    if val != 0:
                        # feature index (row) = j+1, cell index (column) = nbcds
                        mf.write(f"{j+1} {nbcds} {val}\n")
                        nnz += 1
                if nnz > 0:
                    bf.write(f"{bcd}\n")
                    total_nnz += nnz
                else:
                    nbcds -= 1  # rollback: drop empty cells from the barcode index
                nlines += 1
                if nlines % 100000 == 0:
                    custom_log(f"Processed {nlines} cells...")

        # MatrixMarket header: <#features> <#cells> <#non-zero entries>. The nnz field
        # must be the true count of matrix entries, not the number of cells.
        custom_log(f"Writing MEX files to {args.out_dir}....")
        with open(hdr_path, 'wt') as hf:
            hf.write("%%MatrixMarket matrix coordinate integer general\n")
            hf.write(f"{nftrs} {nbcds} {total_nnz}\n")

        result = subprocess.run(f"cat {hdr_path} {body_path} | {gzip_cmd} > {mtx_path}",
                                shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            sys.exit(f"ERROR during MEX matrix compression: {result.stderr.strip()}")
        os.remove(body_path)
        os.remove(hdr_path)
        custom_log(f"Finished writing MEX files: {nbcds} cells, {nftrs} features, {total_nnz} non-zero entries")
    custom_log("convert_cellxgene finished")


if __name__ == "__main__":
    convert_cellxgene(sys.argv[1:])
