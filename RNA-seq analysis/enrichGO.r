print_help <- function() {
    cat("
This script is used to perform GO enrichment analysis was performed using clusterProfiler.

Usage: Rscript enrichGO.r <gene2go_lib> <genelist> <ouput> 

Arguments:
  <gene2go_lib>  The relationship between GeneID, GO pathway ID, GO pathway description and Ontology.
  <genelist>     List of genes to be annotated.
  <output>       The annotion results.

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

gene2go_lib <- commandArgs(trailingOnly = TRUE)[1]
genelist <- commandArgs(trailingOnly = TRUE)[2]
output <- commandArgs(trailingOnly = TRUE)[3]

library (clusterProfiler)
library (ggplot2)
library (dplyr)
library (ggrepel)

go_anno <- read.table(gene2go_lib,sep="\t",header = T)
gene_select <- read.delim(file = genelist, stringsAsFactors = FALSE,header = F)$V1
go_rich <- enricher(gene = gene_select,
                     TERM2GENE = go_anno[c('GO', 'GID')], 
                     TERM2NAME = go_anno[c('GO', 'Description')], 
                     pvalueCutoff = 1, 
                     pAdjustMethod = 'BH', 
                     qvalueCutoff = 1)
go2name <- go_anno %>% select(GO, Ontology)
go2name = unique(go2name)
colnames(go2name) <- c ('ID', 'Ontology')
df_go <- merge (go_rich, go2name[c ('ID', 'Ontology')], by = 'ID')
df_go <- df_go %>%
  relocate(., all_of(colnames(df_go)[13]),.after = 1)
df_go <- df_go[order(df_go$pvalue), ]
write.table(df_go, output, sep = '\t', row.names = FALSE, quote = FALSE)