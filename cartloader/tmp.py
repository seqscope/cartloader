    # notes:
    #   - visium_hd: units_per_um: um_x = ssr.cur_sbcd.px / units_per_um;
    #   - xenium: skip. The input x y are in um, no need to convert.
    #   - stereoseq: $F[{args.col_x}]/{args.scale_factor},$F[{args.col_y}]/{args.scale_factor}


# def clean_feature_by_geneinfo(cmds, args):
#     # feature
#     cmds = cmd_separator(cmds, f"Filtering feature by reference {args.geneinfo}...")
#     cleanftr_cmd = clean_feature_by_ref(args)
#     cmds.append(cleanftr_cmd)

def create_cleanftr(cmds, args):
    # kept_gene_type=$(echo "{params.kept_gene_type}" | sed 's/,/\|/')
    # rm_gene_regex=$(echo "{params.rm_gene_regex}" | sed 's/\^/\\t/g')
    # echo -e "gene_id\tgene\tgn\tgt\tspl\tunspl\tambig" > {sdgeAR_ftr_tabqc_unzip}
    # awk 'BEGIN{{FS=OFS="\t"}} NR==FNR{{ft[$1]=$1; next}} ($1 in ft && $4 + 0 > 50){{print $0}}' \
    #     <(zcat {geneinfo} | grep -P "${{kept_gene_type}}" | cut -f 4 ) \
    #     <(zcat {output.sdgeAR_ftr_tab})| \
    #     grep -vP "${{rm_gene_regex}}" >> {sdgeAR_ftr_tabqc_unzip}
    # gzip -f {sdgeAR_ftr_tabqc_unzip}
    script_path = f"{args.out_dir}/write_cleanftr.sh"
    with open(script_path, "w")  as f:
        f.write(r"""#!/bin/bash
input=$1
output=$2
geneinfo=$3
kept_gene_type=$4
rm_gene_regex=$5

# remove the ".gz" in the output file
output_unzip=${output%.gz}                

echo -e "gene_id\tgene\tgn\tgt\tspl\tunspl\tambig" > ${output_unzip}
kept_gene_type=$(echo "${kept_gene_type}" | sed 's/,/\|/')
rm_gene_regex=$(echo "${rm_gene_regex}" | sed 's/\^/\\t/g')
awk 'BEGIN{FS=OFS="\t"} NR==FNR{ft[$1]=$1; next} ($1 in ft && $4 + 0 > 50){print $0}' \
    <(zcat ${geneinfo} | grep -P "${kept_gene_type}" | cut -f 4 ) \
    <(zcat ${input})| \
    grep -vP "${rm_gene_regex}" >> ${output_unzip}
gzip -f ${output_unzip}
rm output_unzip
""") 
    cmds=cmd_separator(cmds, f"Creating clean feature file to {args.out_feature}...")
    cmds.append(f"bash {script_path} {args.out_transcript} {args.out_feature} {args.kept_gene_type} {args.rm_gene_regex} {args.geneinfo}")
    return cmds

