import os
import pandas as pd
import re
import numpy as np

aux_sge_args = {
    "out": [
        'units_per_um', 'precision_um',
        'colname_x', 'colname_y', #'colname_count',
        'colname_feature_name' #, 'colname_feature_id'
    ],
    "inbcd": [
        'icol_bcd_barcode', 'icol_bcd_x', 'icol_bcd_y'
    ],
    "inftr": [
        'icol_ftr_id', 'icol_ftr_name'
    ],
    "inmtx": [
        'icols_mtx'
    ],
    "inpos": [
        'pos_colname_barcode', 'pos_colname_x', 'pos_colname_y', 'pos_delim'
    ],
    "incsv": [
        'csv_comment', 'csv_delim', 'csv_colname_x', 'csv_colname_y', 'csv_colname_feature_name',
        'csv_colname_count',  'csv_colnames_others', # 'csv_colname_feature_id',
        'csv_colname_phredscore', 'min_phred_score' #, 'add_molecule_id'
    ],
    "spatula": [
        'allow_duplicate_gene_names' #, 'print_feature_id'
    ],
    "ftrname": [
       # 'include_feature_list', 'exclude_feature_list','include_feature_regex', 'exclude_feature_regex'
       'include_feature_regex', 'exclude_feature_regex'
    ],
    # "ftrtype": [
    #     'include_feature_type_regex', 'csv_colname_feature_type', 'feature_type_ref', "feature_type_ref_colidx_name","feature_type_ref_colidx_type"
    # ]
}

def input_by_platform(args):
    if args.platform == "10x_visium_hd":
        in_dict = {
            "mex_bcd": os.path.join(args.in_mex, args.mex_bcd),
            "mex_ftr": os.path.join(args.in_mex, args.mex_ftr),
            "mex_mtx": os.path.join(args.in_mex, args.mex_mtx),
            "pos_parquet": args.pos_parquet,
        }
    elif args.platform == "seqscope":
        in_dict = {
            "mex_bcd": os.path.join(args.in_mex, args.mex_bcd),
            "mex_ftr": os.path.join(args.in_mex, args.mex_ftr),
            "mex_mtx": os.path.join(args.in_mex, args.mex_mtx),
        }
    elif args.platform in ["10x_xenium", "cosmx_smi", "bgi_stereoseq", "vizgen_merscope", "pixel_seq", "nova_st", "generic"]:
        if args.in_csv is not None:
            in_dict={
                "in_csv": args.in_csv
            }
        elif args.in_parquet is not None:
            in_dict={
                "in_parquet": args.in_parquet
            }
        else:
            raise ValueError(f"Provide --in-csv or --in-parquet for {args.platform}")
    else:
        raise ValueError(f"Unsupported platform: {args.platform}")

    for key, value in in_dict.items():
        assert value is not None, f"Missing input file for --{key}: None provided"
        assert os.path.exists(value), f"Input file --{key} does not exist: {value}"
    
    input_files = list(in_dict.values())

    return input_files


def create_minmax(cmds, args):
    minmax_cmd = f"""{args.gzip} -cd {args.out_dir}/{args.out_transcript} | awk 'BEGIN{{FS=OFS="\\t"}} NR==1{{for(i=1;i<=NF;i++){{if($i=="X")x=i;if($i=="Y")y=i}}print $x,$y;next}}{{print $x,$y}}' | awk -F'\\t' ' BEGIN {{ min1 = "undef"; max1 = "undef"; min2 = "undef"; max2 = "undef"; }} {{ if (NR == 2 || $1 < min1) min1 = $1; if (NR == 2 || $1 > max1) max1 = $1; if (NR == 2 || $2 < min2) min2 = $2; if (NR == 2 || $2 > max2) max2 = $2; }} END {{ print "xmin\\t", min1; print "xmax\\t", max1; print "ymin\\t", min2; print "ymax\\t", max2; }}' > {args.out_dir}/{args.out_minmax}"""
    cmds.append(minmax_cmd)
    return cmds

def update_csvformat_by_platform(args):
    platform_mappings = {
        "10x_xenium": {
            "x": "x_location",
            "y": "y_location",
            "feature_name": "feature_name",
            "count": None,
            "delim": ",",
            "comment": False
        },
        "bgi_stereoseq": {
            "x": "x",
            "y": "y",
            "feature_name": "geneID",
            "count": "MIDCounts",
            "delim": None,
            "comment": False
        },
        "cosmx_smi": {
            "x": "x_global_px",
            "y": "y_global_px",
            "feature_name": "target",
            "count": None,
            "delim": ",",
            "comment": False
        },
        "vizgen_merscope": {
            "x": "global_x",
            "y": "global_y",
            "feature_name": "gene",
            "count": None,
            "delim": ",",
            "comment": False
        },
        "pixel_seq": {
            "x": "xcoord",
            "y": "ycoord",
            "feature_name": "geneName",
            "count": None,
            "delim": None,
            "comment": False
        },
        "nova_st": {
            "x": "x",
            "y": "y",
            "feature_name": "geneID",
            "count": "MIDCount",
            "delim": None,
            "comment": True
        },
        "generic": {
            "x": "X",
            "y": "Y",
            "feature_name": "gene",
            "count": "count",
            "delim": None, # none for using default tab
            "comment": False
        }
    }
    # Update arguments based on platform
    platform_settings = platform_mappings[args.platform]
    if args.csv_colname_x is None:
        args.csv_colname_x = platform_settings["x"]
    if args.csv_colname_y is None:
        args.csv_colname_y = platform_settings["y"]
    if args.csv_colname_feature_name is None:
        args.csv_colname_feature_name = platform_settings["feature_name"]
    if args.csv_colname_count is None:
        args.csv_colname_count = platform_settings["count"]
    if args.csv_delim is None:
        args.csv_delim = platform_settings["delim"]
    if args.csv_comment is False:
        args.csv_comment = platform_settings["comment"]
    return args

# def read_minmax(minmax_f):
#     minmax = {}
#     with open(minmax_f, "r") as f:
#         for line in f:
#             key, value = line.strip().split("\t")
#             minmax[key] = float(value)
#     return minmax