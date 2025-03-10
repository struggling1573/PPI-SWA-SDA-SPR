print_help <- function() {
    cat("
This script is used to perform differential gene analysis without biological replicates with edgeR.
We set the bcv value to 0.1 by default, if you have other special requirements, please modify the script.

Usage: Rscript script.r <input_file> <output_file> <fc_cutoff> <pvalue>

Arguments:
  <input_file>  The input file path, which should be a tab-delimited text file with gene expression data.
                The first row should contain column names and the first column should contain gene names.
  <output_file> The output file path where the analysis results will be saved in TXT format.
  <fc_cutoff>   The cutoff value of logFC which will to be used to class the gene to 'up', 'down' or 'normal'.
  <pvalue>      The cutoff value of pvalue.

Author: Hangyu Wang
Date: Dec 26 2024
Unit: Southwest University
Contact: wanghyx666@163.com
")
}

args <- commandArgs()
if ("-h" %in% args) {
    print_help()
    q(save = "no")
}

if (length(commandArgs(trailingOnly = TRUE)) < 4) {
    stop("Please provide both input and output file paths. Use -h for help.")
}

input <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]
fc_cutoff <- commandArgs(trailingOnly = TRUE)[3]
pvalue <- commandArgs(trailingOnly = TRUE)[4]

library(edgeR)

exprSet <- read.table(file = input, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
group_list <- factor(c(rep("Contral", 1), rep("Treat", 1)))
exprSet <- DGEList(counts = exprSet, group = group_list)
bcv = 0.1
et <- exactTest(exprSet, dispersion = bcv^2)
DEG_edgeR=as.data.frame(topTags(et, n = nrow(exprSet$counts)))
DEG_edgeR$GeneID <- rownames(DEG_edgeR)
DEG_edgeR <- DEG_edgeR[, c(ncol(DEG_edgeR), 1:(ncol(DEG_edgeR)-1))]
fc_cutoff <- as.numeric(fc_cutoff)
pvalue <- as.numeric(pvalue)
DEG_edgeR$regulated <- "normal"
loc_up <- intersect(which(DEG_edgeR$logFC>log2(fc_cutoff)),
                    which(DEG_edgeR$PValue<pvalue))
loc_down <- intersect(which(DEG_edgeR$logFC < (-log2(fc_cutoff))),
                      which(DEG_edgeR$PValue<pvalue))
DEG_edgeR$regulated[loc_up] <- "up"
DEG_edgeR$regulated[loc_down] <- "down" 

write.table(DEG_edgeR, output, quote = FALSE, sep = "\t")
