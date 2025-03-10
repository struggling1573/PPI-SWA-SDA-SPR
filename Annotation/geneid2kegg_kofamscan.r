print_help <- function() {
    cat("
This script is used to mapping the GeneID to KEGG pathway ID to KEGG pathway description with the annotion results from kofamscan software.

Usage: Rscript geneid2kegg_kofamscan.r <gene_eggnog> <kegg_lib> <ouput> 

Arguments:
  <gene_kofa>    The results from kofamscan software.
  <kegg_lib>     The KEGG file.
  <output>       The gene2kegg annotion results.

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

gene_kofa <- commandArgs(trailingOnly = TRUE)[1]
kegg_lib <- commandArgs(trailingOnly = TRUE)[2]
gene_kegg <- commandArgs(trailingOnly = TRUE)[3]

library(dplyr)
library(stringr)
options(stringsAsFactors = F)

gene2ko<-read.delim(gene_kofa,header=F,sep="\t")
gene2ko[gene2ko==""] <- NA
colnames(gene2ko) <- c("GID", "KO")
pathway2name <- read.delim (kegg_lib, header = T)
kegg_anno<- merge (gene2ko, pathway2name,by = 'KO')

write.table (kegg_anno, gene_kegg, quote = FALSE, row.names = F, sep = "\t")
