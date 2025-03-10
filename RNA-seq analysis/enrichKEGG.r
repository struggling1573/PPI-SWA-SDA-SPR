print_help <- function() {
    cat("
This script is used to perform KEGG enrichment analysis was performed using clusterProfiler.

Usage: Rscript enrichKEGG.r <gene2kegg_lib> <genelist> <ouput> 

Arguments:
  <gene2kegg_lib>  The relationship between GeneID, KEGG pathway ID and KEGG pathway description, etc.
  <genelist>       List of genes to be annotated.
  <output>         The annotion results.

Author: Hangyu Wang
Date: Dec 29 2024
Unit: Southwest University
Contact: wanghyx666@163.com
")
}

args <- commandArgs()
if ("-h" %in% args) {
    print_help()
    q(save = "no")
}

if (length(commandArgs(trailingOnly = TRUE)) < 3) {
    stop("Please provide both input and output file paths. Use -h for help.")
}

gene2kegg_lib <- commandArgs(trailingOnly = TRUE)[1]
genelist <- commandArgs(trailingOnly = TRUE)[2]
output <- commandArgs(trailingOnly = TRUE)[3]

library (clusterProfiler)
library (ggplot2)
library (dplyr)
library (ggrepel)

kegg_anno <- read.table(gene2kegg_lib,sep="\t",header = T)
gene_select <- read.delim(file = genelist, stringsAsFactors = FALSE,header = F)$V1
kegg_rich <- enricher(gene = gene_select,
                       TERM2GENE = kegg_anno[c('KO', 'GID')], 
                       TERM2NAME = kegg_anno[c('KO', 'KO_des')], 
                       pvalueCutoff = 1, 
                       pAdjustMethod = 'BH', 
                       qvalueCutoff = 1)
df_kegg <- kegg_rich@result
df_kegg <- df_kegg[order(df_kegg$pvalue), ]
head(df_kegg)
write.table (df_kegg, output, quote = FALSE, row.names = F, sep = "\t")
