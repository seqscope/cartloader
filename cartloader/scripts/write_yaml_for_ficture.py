import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, yaml, glob
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                     description="""
                                     Summarize the ficture input/output files into a yaml file.
                                     Two options are available:
                                     (1) Specifying the FICTURE parameters (train-width, n-factor, fit-width, anchor-res) to summarize the output files the out-yaml file. For example:
                                        cartloader write_yaml_for_ficture --out-dir /path/to/outdir --in-transcript /path/to/transcripts.unsorted.tsv.gz --in-cstranscript /path/to/transcripts.sorted.tsv.gz --in-minmax /path/to/coordinate_minmax.tsv --in-feature /path/to/feature.clean.tsv.gz --data-id data_id --platform platform --train-width 12 --n-factor 6 --fit-width 12 --anchor-res 4 --out-yaml /path/to/ficture.yaml
                                     (2) Use --all-ficture to summarize all available output files from the out-dir in one yaml file. For example:
                                        cartloader write_yaml_for_ficture --out-dir /path/to/outdir --in-transcript /path/to/transcripts.unsorted.tsv.gz --in-cstranscript /path/to/transcripts.sorted.tsv.gz --in-minmax /path/to/coordinate_minmax.tsv --in-feature /path/to/feature.clean.tsv.gz --data-id data_id --platform platform --all-ficture --out-yaml /path/to/ficture.all.yaml
                                     """)

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/Output parameters")
    inout_params.add_argument('--out-dir', required= True, type=str, help='Output directory')
    inout_params.add_argument('--out-yaml', type=str, default=None, help='The output YAML file summarizing the output files')
    inout_params.add_argument('--in-transcript', type=str, default=None, help='Input unsorted transcript-indexed SGE file in TSV format, e.g., transcripts.unsorted.tsv.gz')
    inout_params.add_argument('--in-cstranscript', type=str, default=None, help='Input transcript-indexed SGE file that sorted by x and y coordinates in TSV format, e.g., transcripts.sorted.tsv.gz')
    inout_params.add_argument('--in-minmax', type=str, default=None, help='Input coordinate minmax TSV file, e.g., coordinate_minmax.tsv')
    inout_params.add_argument('--in-feature', type=str, default=None,  help='Input TSV file that specify which genes to use as input, e.g., feature.clean.tsv.gz')

    #params setting
    key_params = parser.add_argument_group("Key Parameters", "Key parameters that requires user's attention")
    key_params.add_argument('--data-id', type=str, default=None, help='Data ID or Human-readable Description')
    key_params.add_argument('--platform', type=str, default=None, help='Platform')
    key_params.add_argument('--all-ficture', action='store_true', default=False, help='Use glob to find summarize all available output files from the out-dir in one yaml file. When applied, no need to provide --train-width,--n-factor, --fit-width, --anchor-res. Default: False')
    key_params.add_argument('--train-width', type=int, default=None, help='Hexagon flat-to-flat width (in um) during training. Use comma to specify multiple values')
    key_params.add_argument('--n-factor', type=int, default=None, help='Number of factors to train. Use comma to specify multiple values')
    key_params.add_argument('--fit-width', type=int, default=None, help='Hexagon flat-to-flat width (in um) during model fitting (default: same to train-width)')
    key_params.add_argument('--anchor-res', type=int, default=4, help='Anchor resolution for decoding')

    #action if the out yaml file exists (overwrite or add to the existing file)
    parser.add_argument('--append', action='store_true', default=False, help='Append to the existing yaml file if it exists. If not, overwrite the existing file')
    parser.add_argument('--field-type', type=str, default="dict", choices=["list","dict"], help='How to organize the output structure in the YAML file. If --all-ficture is used, it has to be defined as list. By default, when --all-ficture is not used, it is defined as dict. Options: list, dict')
  
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    
    return parser.parse_args(_args) 

def update_yaml_for_input(yaml_content, args):
  input_field={
    'tsvsge_unsorted': args.in_transcript,
    'tsvsge_sorted': args.in_cstranscript,
    'feature': args.in_feature,
    'xy_range': args.in_minmax,
  }
  yaml_content['input'].update(input_field)

def ini_yaml(args):
  yaml_content = {
  'data_id': args.data_id,
  'platform': args.platform,
  'input': {},
  'output': {},
  }
  return yaml_content

def valid_yaml_field(key, config, format):
    if key not in config:
        if format == 'list':
          config[key] = []
        elif format == 'dict':
          config[key] = {}

# for specific chars
def files_existence_by_fn(fnlist, out_dir):
  for fn in fnlist:
    if not os.path.exists(os.path.join(out_dir, fn)):
      print(f"File {fn} does not exist in {out_dir}")
      return False
  return True

def update_yaml_for_hexagon(yaml_content, out_dir, train_width, field_type):
  hexagon_fn= f"hexagon.d_{train_width}.tsv.gz"
  if files_existence_by_fn([hexagon_fn], out_dir):
    valid_yaml_field('hexagon', yaml_content['output'], field_type)
    hexagon_entry= {
      'hexagon_width': train_width,
      'hexagon_sge': hexagon_fn,
    }
    if field_type == "list":
      yaml_content['output']['hexagon'].append(hexagon_entry)
    elif field_type == "dict":
      yaml_content['output']['hexagon'] = hexagon_entry

def update_yaml_for_lda(yaml_content, out_dir, train_width, n_factor, field_type):
  # prefix
  prefix = f"nF{n_factor}.d_{train_width}"
  color_prefix = prefix
  # fn
  fitres_fn = f"{prefix}.fit_result.tsv.gz"
  de_fn     = f"{prefix}.bulk_chisq.tsv"
  coarse_fn = f"{prefix}.coarse.png"
  top_fn    = f"{prefix}.coarse.top.png"
  rgb_fn    = f"{color_prefix}.rgb.tsv"
  cbar_fn   = f"{color_prefix}.cbar.png"
  rep_fn    = f"{prefix}.factor.info.html"
  if files_existence_by_fn([fitres_fn, rgb_fn, de_fn, coarse_fn, top_fn, cbar_fn, rep_fn], out_dir):
    valid_yaml_field('lda', yaml_content['output'], field_type)
    lda_entry= {
      # params
      'train_width': train_width,
      'n_factor': n_factor,
      # the deployed results in novaflow
      'fit_results': fitres_fn, 
      'cmap': rgb_fn, 
      'de': de_fn,
      # other results may be useful
      'coarse': coarse_fn,
      'top': top_fn,
      'cbar': cbar_fn,
      'report': rep_fn,
    }
    if field_type == "list":
      yaml_content['output']['lda'].append(lda_entry)
    elif field_type == "dict":
      yaml_content['output']['lda'] = lda_entry

def update_yaml_for_transform(yaml_content, out_dir, train_width, n_factor, fit_width, anchor_res, field_type):
  # prefix
  prefix      = f"nF{n_factor}.d_{train_width}.prj_{fit_width}.r_{anchor_res}"
  color_prefix  = f"nF{n_factor}.d_{train_width}"
  # fn
  fitres_fn   = f"{prefix}.fit_result.tsv.gz"
  de_fn       = f"{prefix}.bulk_chisq.tsv"
  coarse_fn   = f"{prefix}.coarse.png"
  top_fn      = f"{prefix}.coarse.top.png"
  rgb_fn      = f"{color_prefix}.rgb.tsv"
  cbar_fn     = f"{color_prefix}.cbar.png"
  rep_fn      = f"{prefix}.factor.info.html"
  if files_existence_by_fn([fitres_fn, rgb_fn, de_fn, coarse_fn, top_fn, cbar_fn, rep_fn], out_dir):
    valid_yaml_field('transform', yaml_content['output'], field_type)
    transform_entry= {
      # params
      'train_width': train_width,
      'n_factor': n_factor,
      'fit_width': fit_width,
      'anchor_res': anchor_res,
      # the deployed results in novaflow
      'fit_results': fitres_fn, 
      'cmap': rgb_fn, 
      'de': de_fn,
      # other results may be useful
      'coarse': coarse_fn,
      'top': top_fn,
      'cbar': cbar_fn,
      'report': rep_fn,
    }
    if field_type == "list":
      yaml_content['output']['transform'].append(transform_entry)
    elif field_type == "dict":
      yaml_content['output']['transform'] = transform_entry

def update_yaml_for_decode(yaml_content, out_dir, train_width, n_factor, fit_width, anchor_res, radius, field_type):
  # prefix
  prefix        = f"nF{n_factor}.d_{train_width}.decode.prj_{fit_width}.r_{anchor_res}_{radius}"
  color_prefix  = f"nF{n_factor}.d_{train_width}"
  # fn
  pixel_fn    = f"{prefix}.pixel.sorted.tsv.gz"
  de_fn       = f"{prefix}.bulk_chisq.tsv"
  pixel_fn    = f"{prefix}.pixel.png"
  rgb_fn      = f"{color_prefix}.rgb.tsv"
  cbar_fn     = f"{color_prefix}.cbar.png"
  rep_fn      = f"{prefix}.factor.info.html"
  if files_existence_by_fn([pixel_fn, rgb_fn, de_fn], out_dir):
    valid_yaml_field('decode', yaml_content['output'], field_type)
    decode_entry= {
      # params
      'train_width': train_width,
      'n_factor': n_factor,
      'fit_width': fit_width,
      'anchor_res': anchor_res,
      'radius': radius,
      # the deployed results in novaflow
      'pixel_sorted': pixel_fn, 
      'cmap': rgb_fn, 
      'de': de_fn,
      # other results may be useful
      'pixel': pixel_fn,
      'cbar': cbar_fn,
      'report': rep_fn,
    }
    if field_type == "list":
      yaml_content['output']['decode'].append(decode_entry)
    elif field_type == "dict":
      yaml_content['output']['decode'] = decode_entry

# for existing files
def extract_chars_from_existing(out_dir, regex, fn_pattern):
    # Extract all number groups from files matching the pattern
    charlist = [
        tuple(int(group) for group in match.groups())
        for file in glob.glob(os.path.join(out_dir, fn_pattern))
        if (match := regex.match(os.path.basename(file)))
    ]
    return charlist

def update_yaml_for_existing(yaml_content, out_dir, field_type="list"):
  lda_charlist = extract_chars_from_existing(out_dir, regex=re.compile(r"nF(\d+)\.d_(\d+)\.fit_result\.tsv\.gz"), fn_pattern="nF*.d_*.fit_result.tsv.gz",)
  if len(lda_charlist) > 0:
    for train_width, n_factor in lda_charlist:
      update_yaml_for_hexagon(yaml_content, out_dir, train_width, field_type)
      update_yaml_for_lda(yaml_content, out_dir, train_width, n_factor, field_type)
  tsf_charlist = extract_chars_from_existing(out_dir, regex=re.compile(r"nF(\d+)\.d_(\d+)\.prj_(\d+)\.r_(\d+)\.fit_result\.tsv\.gz"), fn_pattern="nF*.d_*.prj_*.r_*.fit_result.tsv.gz")
  if len(tsf_charlist) > 0:
    for train_width, n_factor, fit_width, anchor_res in tsf_charlist:
      update_yaml_for_transform(yaml_content, out_dir, train_width, n_factor, fit_width, anchor_res, field_type)
  dc_charlist = extract_chars_from_existing(out_dir, regex=re.compile(r"nF(\d+)\.d_(\d+)\.decode\.prj_(\d+)\.r_(\d+)_(\d+)\.fit_result\.tsv\.gz"), fn_pattern="nF*.d_*.decode.prj_*.r_*.fit_result.tsv.gz")
  if len(dc_charlist) > 0:
    for train_width, n_factor, fit_width, anchor_res, radius in dc_charlist:
      update_yaml_for_decode(yaml_content, out_dir, train_width, n_factor, fit_width, anchor_res, radius, field_type)

def write_yaml_for_ficture(_args):
    args = parse_arguments(_args)
    # check args for all_ficture
    if not args.all_ficture:
      assert args.train_width is not None, "Provide the train width for write_yaml_for_ficture"
      assert args.n_factor is not None, "Provide the number of factors used in LDA train for write_yaml_for_ficture"
      if args.fit_width is None:
        args.fit_width = args.train_width
      radius = args.anchor_res + 1
    # out_yaml
    if args.out_yaml is None:
      if args.all_ficture:
        args.out_yaml= os.path.join(args.out_dir, "ficture.all.yaml")
      else:
        args.out_yaml = os.path.join(args.out_dir, f"ficture.nF{args.n_factor}.d_{args.train_width}.decode.prj_{args.fit_width}.r_{args.anchor_res}_{radius}.yaml") 
    # init yaml by --append or not
    if os.path.exists(args.out_yaml) and args.append:
      with open(args.out_yaml, 'r') as rf:
        yaml_content = yaml.load(rf, Loader=yaml.FullLoader)
    else:
      warnings.warn(f"File {args.out_yaml} exists. Overwriting the file...")
      yaml_content = ini_yaml(args)
    # yaml content
    update_yaml_for_input(yaml_content, args)
    if not args.all_ficture:
      # update the yaml output field
      update_yaml_for_hexagon(yaml_content, args.out_dir, args.train_width, args.field_type)
      update_yaml_for_lda(yaml_content, args.out_dir, args.train_width, args.n_factor, args.field_type)
      update_yaml_for_transform(yaml_content, args.out_dir, args.train_width, args.n_factor, args.fit_width, args.anchor_res, args.field_type)
      update_yaml_for_decode(yaml_content, args.out_dir, args.train_width, args.n_factor, args.fit_width, args.anchor_res, radius, args.field_type)
    else:
      print("Given --all-ficture, define the output structure as list")
      update_yaml_for_existing(yaml_content, args.out_dir, "list")
    # write yaml
    yaml_content_str = yaml.dump(yaml_content, sort_keys=False, default_flow_style=False)
    with open(args.out_yaml, 'w') as file:
        file.write(yaml_content_str)
    # with open(args.out_yaml, 'w') as wf:
    #   yaml.dump(yaml_content_str, wf, default_flow_style=False, sort_keys=False)

if __name__ == "__main__":
    # get the cartloader path
    global cartloader_repo
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
