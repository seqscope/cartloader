suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(argparse)
  library(data.table)
  library(tidyr)
  library(patchwork)
})

log_message <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

parser <- ArgumentParser(description='Perform heatmap analysis between LDA results and clustering results')
parser$add_argument('--results', type="character", required=TRUE, default='Result file in the LDA output')
parser$add_argument('--clust', type="character", required=TRUE, default='Cluster assignments')
parser$add_argument('--de-results', type="character", required=TRUE, default='Result file containing chi-squared values for top genes')
parser$add_argument('--de-clust', type="character", required=TRUE, default='DE results for clusters')
parser$add_argument('--out', type="character", required=TRUE, default='Prefix of the output files')
parser$add_argument('--offset-data', type='integer', default=2, help='Offset for the results file')
parser$add_argument('--cell-clust', type='character', nargs='+',default="cell_id", help='Column name containing the cell IDs in the clustering file')
parser$add_argument('--colname-clust', type='character', required=TRUE, help='Column name containing the cluster names in the clustering file')
parser$add_argument('--top-k', type='integer', default=5, help='Number of top genes to annotate each factor')
parser$add_argument('--min-count', type='integer', default=100, help='Prune factors with fewer than this many cells')
parser$add_argument('--draw', action='store_true',help='Draw the UMAP plot based on the clustering results')
parser$add_argument('--width', type='integer', help='Width of the heatmap PDF')
parser$add_argument('--height', type='integer', help='Height of the heatmap PDF')

anno_top_k_genes <- function(chisqf, top_k) {
    df_chisq = read.table(chisqf, sep="\t", header=TRUE) %>% arrange(desc(Chi2))
    ## rename Feature to gene, if it exists
    if ("Feature" %in% colnames(df_chisq)) {
        df_chisq = df_chisq %>% rename(gene = Feature)
    }
    if (!"gene" %in% colnames(df_chisq)) {
        stop("The chisqf file must contain a 'gene' or 'Feature' column.")
    }
    df_tmp = df_chisq %>% 
        group_by(factor) %>% 
        slice_head(n = top_k)
    unique_factors = unique(df_tmp$factor)
    sprintf_str = "%02d"
    if ( length(unique_factors) >= 100 ) {
        sprintf_str = "%03d"
    } else if ( length(unique_factors) >= 1000 ) {
        sprintf_str = "%04d"
    }
    return(df_tmp %>% 
        summarise(gene_string = paste0(sprintf(sprintf_str, unique(factor)), "_", paste(gene, collapse = "_")), .groups = "drop") %>% 
        pull(gene_string))
}

load_lda_results <- function(tsvf, chisqf, offset_data, top_k, min_count, colnames_cell_clust) {
    ## read the results TSV file
    df_tsv = fread(tsvf, header = TRUE, sep = '\t') 
    ## remove the leading '#' in the first column name if exists
    setnames(df_tsv, 1, sub("^#", "", names(df_tsv)[1]))
    df_tsv = df_tsv %>% select(-any_of(c("random_key", "topK", "topP")))

    ibeg = offset_data

    print(paste("Number of factors:", ncol(df_tsv)-ibeg+1))

    factor_alias = anno_top_k_genes(chisqf, top_k)
    
    mat = as.matrix(df_tsv[, ibeg:(ncol(df_tsv))])
    mat_colsums = colSums(mat)

    print(paste("Number of factors with >100 posterior counts:", sum(mat_colsums>100)))

    #print(head(df_tsv))

    df_tsv = df_tsv[, .SD, .SDcols = c(rep(TRUE,length(colnames_cell_clust)), mat_colsums>=min_count)]  # Select 1st and 3rd columns


    colnames(df_tsv)[1:(ibeg-1)] = colnames_cell_clust
    colnames(df_tsv)[ibeg:(ncol(df_tsv))] = factor_alias[1:(ncol(df_tsv) - ibeg + 1)]
    iend = ncol(df_tsv)

    return(df_tsv)
}

lda_celltype_heatmaps = function(df_tsv, df_anno, celltype_name, offset_data, colnames_cell_clust) {
    ibeg = offset_data
    iend = ncol(df_tsv)
    df_merged = df_tsv %>% inner_join(df_anno, by = all_of(colnames_cell_clust))

    ## compute anno_by_factor
    df_merged %>% select(all_of(ibeg:iend), cluster) %>%
        group_by(cluster) %>%
        summarise(across(everything(), sum), .groups = "drop") -> anno_by_factor

    anno_by_factor_row <- anno_by_factor %>%
        rowwise() %>%
        mutate(across(-cluster, ~ . / sum(c_across(-cluster)))) %>%
        ungroup()

    return( list(count=anno_by_factor, norm_frac=anno_by_factor_row) )
}

args <- parser$parse_args()

log_message("Analysis started")

print(args$cell_clust)

## make tsv file
log_message("Loading LDA results")
df_tsv = load_lda_results(args$results, args$de_results, args$offset_data, args$top_k, args$min_count, args$cell_clust)

## load the clustering results
log_message("Loading clustering results")
df_clust = fread(args$clust)
## remove the leading '#' in the first column name if exists
setnames(df_clust, 1, sub("^#", "", names(df_clust)[1]))

#print(head(df_clust)) 
df_clust = df_clust %>% select(all_of(args$cell_clust), cluster = !!args$colname_clust)
print(head(df_clust)) 


## create heatmaps
log_message("Creating heatmaps")
heatmaps = lda_celltype_heatmaps(df_tsv, df_clust, args$colname_clust, args$offset_data, args$cell_clust)

if ( length(args$de_clust) > 0 ) {
    log_message("Loading DE results for clusters")
    clust_alias = anno_top_k_genes(args$de_clust, args$top_k)

    #print(head(heatmaps))
    heatmaps$count$cluster = clust_alias
    heatmaps$norm_frac$cluster = clust_alias
}

log_message("Saving heatmaps")
write.table(heatmaps$count, file=paste0(args$out, ".counts.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(heatmaps$norm_frac, file=paste0(args$out, ".normfrac.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

if ( args$draw ) {
    if ( is.null(args$width) || is.null(args$height) ) {
        ## automatically calculate the width and height
        max_length_clust_alias = max(nchar(heatmaps$count$cluster))
        max_length_factor_alias = max(nchar(colnames(heatmaps$count)[-1]))
        n_clust = nrow(heatmaps$norm_frac)
        n_factor = ncol(heatmaps$norm_frac) - 1
        args$width = 2 + max_length_factor_alias / 10 + n_clust / 10
        args$height = 2 + max_length_clust_alias / 10 + n_factor / 10 
        print(paste0("Auto-calculated width:", args$width, " height:", args$height))        
    }
    pdf(paste0(args$out, ".pdf"), width = args$width, height = args$height)
    df_melt <- tidyr::pivot_longer(heatmaps$norm_frac, cols = -cluster, names_to = "Row", values_to = "Value") %>% rename(Column = cluster)
    p2 <- ggplot(df_melt, aes(x = as.factor(Column), y = as.factor(Row), fill = log10(Value+1e-4))) +
        geom_tile(color = "white") +
        scale_fill_distiller(palette = "Spectral", direction = -1, name = "log10(frac)") +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),  # rotate x-axis labels
            axis.text.y = element_text(size = 10),              # adjust row label size
            axis.title = element_blank()
        )
    print(p2)
    dev.off()
}