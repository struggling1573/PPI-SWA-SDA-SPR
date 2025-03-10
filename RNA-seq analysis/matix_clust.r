print_help <- function() {
    cat("
This script is used to clust sample by gene expression matrix.

Usage: Rscript matix_clust.r <express_matrix>
Arguments:
  <express_matrix>  The gene expression matrix, sep by tab. Like this:
    #gene_id              ck_1 ck_2 ck_3 ck_4
    #ckhap1_chr01G000330  153  369  523  427
    #ckhap1_chr01G000380  113  0    330  0

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

if (length(commandArgs(trailingOnly = TRUE)) < 1) {
    stop("Please provide both input and output file paths. Use -h for help.")
}

express_matrix <- commandArgs(trailingOnly = TRUE)[1]

library(ggplot2)
library(ggdendro)

exprSet <- read.table(express_matrix, header = T)
rownames(exprSet) <- exprSet[, 1]
exprSet <- exprSet[, -1]
exprSet <- round(exprSet)
sampleDists <- dist(t(exprSet))
hc <- hclust(sampleDists, method = "ward.D2") 
ggdendrogram(hc, rotate = FALSE, size = 1, leaf_labels = FALSE)
ggsave("matrix_clust.png",width = 4,height = 3,dpi=300)

