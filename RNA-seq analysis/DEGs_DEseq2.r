print_help <- function() {
    cat("
This script is used to perform differential gene analysis with DEseq2.

Usage: Rscript DEG_class.r <express_matrix> <sample_group> <output_file> <fc_cutoff> <pvalue>

Arguments:
  <express_matrix>  The gene expression matrix, sep by tab. Like this:
    #gene_id              ck_1 ck_2 ck_3 ck_4
    #ckhap1_chr01G000330  153  369  523  427
    #ckhap1_chr01G000380  113  0    330  0
  <sample_group>    The sample and it's group, sep by tab. Like this:
    #sample    group
    #ck_1  female
    #ck_2  female
    #ck_3  male
    #ck_4  male
  <output_file>     The output file path where the analysis results will be saved in TXT format.
  <fc_cutoff>       The cutoff value of logFC which will to be used to class the gene to 'up', 'down' or 'normal'.
  <pvalue>          The cutoff value of pvalue (no special requirements, the recommended consistency of 0.05).

Author: Hangyu Wang
Date: Dec 28 2024
Unit: Southwest University
Contact: wanghyx666@163.com
")
}

args <- commandArgs()
if ("-h" %in% args) {
    print_help()
    q(save = "no")
}

if (length(commandArgs(trailingOnly = TRUE)) < 5) {
    stop("Please provide both input and output file paths. Use -h for help.")
}

express_matrix <- commandArgs(trailingOnly = TRUE)[1]
sample_group <- commandArgs(trailingOnly = TRUE)[2]
output <- commandArgs(trailingOnly = TRUE)[3]
fc_cutoff <- commandArgs(trailingOnly = TRUE)[4]
pvalue <- commandArgs(trailingOnly = TRUE)[5]

library(ggplot2)
library(DESeq2)

exprSet <- read.table(express_matrix, header = T)
rownames(exprSet) <- exprSet[, 1]
exprSet <- exprSet[, -1]
exprSet <- round(exprSet)

metadata <- read.table(sample_group, header = T)
rownames(metadata) <- metadata[, 1]

dds <-DESeqDataSetFromMatrix(countData=exprSet, 
                             colData=metadata, 
                             design=~group)
dds <- DESeq(dds)
contrast=c("group", unique(metadata$group))
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
dd1Ordered = dd1[order(dd1$pvalue),]
DEG_dd1 = as.data.frame(dd1Ordered)
DEG_dd1$GeneID <- rownames(DEG_dd1)
DEG_dd1 <- DEG_dd1[, c(ncol(DEG_dd1), 1:(ncol(DEG_dd1)-1))]
fc_cutoff <- as.numeric(fc_cutoff)
pvalue <- as.numeric(pvalue)
DEG_dd1$regulated <- "normal"
loc_up <- intersect(which(DEG_dd1$log2FoldChange>log2(fc_cutoff)),
                    which(DEG_dd1$pvalue<pvalue))
loc_down <- intersect(which(DEG_dd1$log2FoldChange < (-log2(fc_cutoff))),
                      which(DEG_dd1$pvalue<pvalue))
DEG_dd1$regulated[loc_up] <- "up"
DEG_dd1$regulated[loc_down] <- "down" 

write.table(DEG_dd1, output, quote = FALSE, row.names = F, sep = "\t")
