print_help <- function() {
    cat("
This script is used to classify the user-provided gene lists into tissue-enriched or tissue-specific transcripts with the TissueEnrich R package.

Usage: Rscript TissueEnrich.r <express_matrix> <sample_group> <output_file>

Arguments:
  <express_matrix>   The gene expression matrix, separated by tab. For example:
    #gene_id              ck_1 ck_2 ck_3 ck_4
    #ckhap1_chr01G000330  153  369  523  427
    #ckhap1_chr01G000380  113  0    330  0
  <sample_group>     The sample and its group, separated by tab. For example:
    #sample    group
    #ck_1  female
    #ck_2  female
    #ck_3  male
    #ck_4  male
  <output_file>      The output file path where the analysis results will be saved in xls format.

Author: Hangyu Wang
Date: Jan 05 2025
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

express_matrix <- commandArgs(trailingOnly = TRUE)[1]
sample_group <- commandArgs(trailingOnly = TRUE)[2]
output <- commandArgs(trailingOnly = TRUE)[3]

library(TissueEnrich)
library(dplyr)
library(tidyr)
library(tibble)

exp_data <- read.table(express_matrix, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
GeneID <- exp_data[, 1]
exp_data <- exp_data[, -1]
rownames(exp_data) <- GeneID
colnames(exp_data) <- gsub("[-.]", "_", colnames(exp_data))

cat("Sample names of the processed gene expression matrix:\n")
print(colnames(exp_data))

classify <- read.table(sample_group, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(classify) <- c("Sample", "Group")
classify$Sample <- gsub("[-]", "_", classify$Sample)
rownames(classify) <- classify$Sample

cat("Sample names of the processed sample group information:\n")
print(rownames(classify))

available_samples <- intersect(rownames(classify), colnames(exp_data))
if (length(available_samples) == 0) {
    stop("No matching sample names. Please check the sample name format in the input files.")
}

exp_data <- exp_data[, available_samples, drop = FALSE]
classify <- classify[available_samples, , drop = FALSE]

cat("First few columns of the processed gene expression matrix:\n")
print(head(exp_data))
cat("Processed sample group information (show all):\n")
print(classify)

exp_data_temp1 <- exp_data %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression")
merged_data <- left_join(exp_data_temp1, classify, by = "Sample")
average_expression <- merged_data %>%
  group_by(Gene, Group) %>%
  summarize(Average_Expression = mean(Expression))
exp_data_temp2 <- average_expression %>%
  pivot_wider(names_from = Group, values_from = Average_Expression)
average_expression_matrix <- as.matrix(exp_data_temp2[, -1]) 
average_expression_matrix <- round(average_expression_matrix, 2)
rownames(average_expression_matrix) <- exp_data_temp2$Gene
Summarized_data <- SummarizedExperiment(assays = SimpleList(average_expression_matrix),
                                        rowData = row.names(average_expression_matrix),
                                        colData = colnames(average_expression_matrix))
num_tissues <- ncol(average_expression_matrix)
print(num_tissues)
maxNumberOfTissues <- min(num_tissues - 1, 7) 
GeneRetrieval_output <- teGeneRetrieval(Summarized_data,
                                        foldChangeThreshold = 5,
                                        maxNumberOfTissues = maxNumberOfTissues,
                                        expressedGeneThreshold = 1)
result <- as.data.frame(assay(GeneRetrieval_output))
write.table(result, output, quote = F, sep = '\t', row.names = F)
