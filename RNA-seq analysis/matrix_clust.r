print_help <- function() {
    cat("
This script is used to clust sample by gene expression matrix.

Usage: Rscript matix_clust.r <express_matrix> [<samples>]
Arguments:
  <express_matrix>  The gene expression matrix, sep by tab. Like this:
    #gene_id              ck_1 ck_2 ck_3 ck_4
    #ckhap1_chr01G000330  153  369  523  427
    #ckhap1_chr01G000380  113  0    330  0
  <samples>         Comma-separated list of samples to include in clustering. (Optional)

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
    stop("Please provide the expression matrix file path. Use -h for help.")
}

# 获取输入文件路径
express_matrix <- commandArgs(trailingOnly = TRUE)[1]

# 检查是否提供了第二个参数
if (length(commandArgs(trailingOnly = TRUE)) > 1) {
    samples <- strsplit(commandArgs(trailingOnly = TRUE)[2], ",")[[1]]
    output_file <- "matrix_clust_severalSamples.png"
} else {
    samples <- NULL
    output_file <- "matrix_clust_allSamples.png"
}

library(ggplot2)

exprSet <- read.table(express_matrix, header = T)
rownames(exprSet) <- exprSet[, 1]
exprSet <- exprSet[, -1]
exprSet <- round(exprSet)
print(colnames(exprSet))

colnames(exprSet) <- trimws(tolower(colnames(exprSet)))
if (!is.null(samples)) {
    samples <- trimws(tolower(samples))
    samples <- gsub("-", ".", samples)
    exprSet_selected <- exprSet[, samples, drop = FALSE]
    missing_samples <- setdiff(samples, colnames(exprSet))
    if (length(missing_samples) > 0) {
        stop(paste("The following samples are not found in the expression matrix:", paste(missing_samples, collapse = ", ")))
    }
} else {
    exprSet_selected <- exprSet
}

sampleDists <- dist(t(exprSet_selected))
hc <- hclust(sampleDists, method = "ward.D2")

png(output_file, width = 4, height = 3, units = "in", res = 300)
plot(hc, main = "Sample Clustering Dendrogram", xlab = "", ylab = "Distance", sub = "")
dev.off()
