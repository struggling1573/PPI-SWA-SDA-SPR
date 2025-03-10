print_help <- function() {
    cat("
This script is used to convert the gene expression raw count matrix into TPM and FPKM matrices

Usage: Rscript reads2TPM.r <gtf> <raw_count> <prefix>

Arguments:
  <gtf>       Genome annotation information in gff format, which needs to be consistent with the documents used in previous analyses.
  <raw_count> The gene raw count matrix.
  <prefix>    Prefix of output.

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

if (length(commandArgs(trailingOnly = TRUE)) < 3) {
    stop("Please provide both input and output file paths. Use -h for help.")
}

library(GenomicFeatures)

gtffile <- commandArgs(trailingOnly = TRUE)[1]
readsfile <- commandArgs(trailingOnly = TRUE)[2]
out_name <- commandArgs(trailingOnly = TRUE)[3]
Counts <- read.table(file = readsfile, header = TRUE, row.names = 1)
txdb <- makeTxDbFromGFF(gtffile, format="gtf")
ebg <- exonsBy(txdb, by="gene")
ebgList <- sum(width(reduce(ebg)))
genes <- intersect(rownames(Counts), names(ebgList))
Length <- as.vector(ebgList[genes])
TPM <- t(t(Counts / t(Length)) * 1e6 / colSums(Counts / t(Length)))
FPKM <- t(t(Counts / t(Length)) * 1e9 / colSums(Counts))
TPM <- as.data.frame(TPM)
TPM$GeneID <- rownames(TPM)
TPM <- TPM[, c(ncol(TPM), 1:(ncol(TPM)-1))]
FPKM <- as.data.frame(FPKM)
FPKM$GeneID <- rownames(FPKM)
FPKM <- FPKM[, c(ncol(FPKM), 1:(ncol(FPKM)-1))]
outTPM <- paste(out_name, "_TPM.xls", sep = "") 
outFPKM <- paste(out_name, "_FPKM.xls", sep = "") 

write.table(TPM, file=outTPM, quote = FALSE, row.names = F, sep = "\t")
write.table(FPKM, file=outFPKM , quote = FALSE, row.names = F, sep = "\t")
