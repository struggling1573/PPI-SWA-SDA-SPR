print_help <- function() {
    cat("
This script is used to perform differential gene analysis with DESeq2.

Usage: Rscript DEGs_DEseq2.r <express_matrix> <sample_group> <output_file> <fc_cutoff> <pvalue> <library_size_estimation> <low_expression_filter>

Arguments:
  <express_matrix>  The gene expression matrix, separated by tab. For example:
    #gene_id              ck_1 ck_2 ck_3 ck_4
    #ckhap1_chr01G000330  153  369  523  427
    #ckhap1_chr01G000380  113  0    330  0
  <sample_group>    The sample and its group, separated by tab. For example:
    #sample    group
    #ck_1  female
    #ck_2  female
    #ck_3  male
    #ck_4  male
  <output_file>     The output file path where the analysis results will be saved in TXT format.
  <fc_cutoff>       The cutoff value of logFC which will be used to classify the gene as 'up', 'down' or 'normal'.
  <pvalue>          The cutoff value of p - value (no special requirements, the recommended value is 0.05).
  <library_size_estimation> y/n(default: y) Whether to estimate the sequencing library size.
  <low_expression_filter> y/n(default: y) Whether to filter low - expression genes. Low - expressed genes are those whose expression counts add up to less than 10 in all samples.

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
library_size_estimation <- ifelse(length(commandArgs(trailingOnly = TRUE)) >= 6, commandArgs(trailingOnly = TRUE)[6], "y")
low_expression_filter <- ifelse(length(commandArgs(trailingOnly = TRUE)) >= 7, commandArgs(trailingOnly = TRUE)[7], "y")

library(ggplot2)
library(DESeq2)

exprSet <- read.table(express_matrix, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
GeneID <- exprSet[, 1]
exprSet <- exprSet[, -1]
rownames(exprSet) <- GeneID
colnames(exprSet) <- gsub("[-.]", "_", colnames(exprSet))

cat("Processed sample names in the gene expression matrix:\n")
print(colnames(exprSet))

metadata <- read.table(sample_group, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata$sample <- gsub("[-]", "_", metadata$sample)
rownames(metadata) <- metadata$sample
cat("Processed sample names in the sample grouping information:\n")
print(rownames(metadata))
available_samples <- intersect(rownames(metadata), colnames(exprSet))
if (length(available_samples) == 0) {
    stop("No matching sample names. Please check the sample name format in the input files.")
}
exprSet <- exprSet[, available_samples, drop = FALSE]
metadata <- metadata[available_samples, , drop = FALSE]
exprSet <- as.matrix(exprSet)
storage.mode(exprSet) <- "integer" 
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = metadata,
                              design = ~ group)

if (library_size_estimation == "y") {
    dds <- estimateSizeFactors(dds)
    cat("estimateSizeFactors:\n")
    print(sizeFactors(dds))
}

if (low_expression_filter == "y") {
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
}

dds <- DESeq(dds)

contrast_groups <- c("group", unique(metadata$group))
results <- results(dds, contrast = contrast_groups, alpha = 0.05)

ordered_results <- results[order(results$pvalue), ]
DEG_results <- as.data.frame(ordered_results)
DEG_results$GeneID <- rownames(DEG_results)
DEG_results <- DEG_results[, c("GeneID", setdiff(names(DEG_results), "GeneID"))]

fc_cutoff <- as.numeric(fc_cutoff)
pvalue_cutoff <- as.numeric(pvalue)
DEG_results$regulated <- "normal"

up_genes <- which(
  DEG_results$log2FoldChange > log2(fc_cutoff) & 
  DEG_results$pvalue < pvalue_cutoff
)

down_genes <- which(
  DEG_results$log2FoldChange < -log2(fc_cutoff) & 
  DEG_results$pvalue < pvalue_cutoff
)

DEG_results$regulated[up_genes] <- "up"
DEG_results$regulated[down_genes] <- "down"

write.table(DEG_results, output, sep = "\t", quote = FALSE, row.names = FALSE)
