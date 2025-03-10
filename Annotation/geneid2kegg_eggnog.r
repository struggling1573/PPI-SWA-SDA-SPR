print_help <- function() {
    cat("
This script is used to mapping the GeneID to KEGG pathway ID to KEGG pathway description with the annotion results from EggNOG software.

Usage: Rscript geneid2kegg_eggnog.r <gene_eggnog> <kegg_lib> <ouput> 

Arguments:
  <gene_eggnog>  The results from EggNOG software.
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

gene_eggnog <- commandArgs(trailingOnly = TRUE)[1]
kegg_lib <- commandArgs(trailingOnly = TRUE)[2]
gene_kegg <- commandArgs(trailingOnly = TRUE)[3]

library(dplyr)
library(stringr)
options(stringsAsFactors = F)

egg<-read.delim(gene_eggnog,header=T,sep="\t")
egg[egg==""] <- NA
koterms <- egg %>%dplyr::select(query, KEGG_ko) %>% na.omit() %>% filter(str_detect(KEGG_ko,"ko"))
all_ko_list <- str_split(koterms$KEGG_ko,",")
gene2ko <- data.frame(GID = rep(koterms$query, times = sapply(all_ko_list, length)), KEGG_ko = unlist(all_ko_list)) %>% filter(str_detect(KEGG_ko,"ko")) 
gene2ko$KO <- str_replace (gene2ko$KEGG_ko, "ko:","")
pathway2name <- read.delim (kegg_lib, header = T)
kegg_anno<- merge (gene2ko, pathway2name,by = 'KO')

write.table (kegg_anno, gene_kegg, quote = FALSE, row.names = F, sep = "\t")
