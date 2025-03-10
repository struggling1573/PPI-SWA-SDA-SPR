print_help <- function() {
    cat("
This script is used to mapping the GeneID to GO pathway ID to GO pathway description with the annotion results from EggNOG software.

Usage: Rscript geneid2go_eggnog.r <gene_eggnog> <go_lib> <ouput> 

Arguments:
  <gene_eggnog>  The results from EggNOG software.
  <go_lib>       The GO file.
  <output>       The gene2go annotion results.

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
go_lib <- commandArgs(trailingOnly = TRUE)[2]
gene_go <- commandArgs(trailingOnly = TRUE)[3]

library(dplyr)
library(stringr)
options(stringsAsFactors = F)

egg<-read.delim(gene_eggnog,header=T,sep="\t")
egg[egg==""] <- NA
colnames(egg)
goterms <- egg %>%dplyr::select(query, GOs) %>% na.omit() %>% filter(str_detect(GOs,"GO")) 
all_go_list <- str_split(goterms$GOs,",") 
gene2go <- data.frame(GID = rep(goterms$query, times = sapply(all_go_list, length)), GO = unlist(all_go_list)) %>% filter(str_detect(GO,"GO")) 
go2name <- read.delim (go_lib, header = FALSE, stringsAsFactors = FALSE)
names (go2name) <- c ('GO', 'Description', 'Ontology')
go_anno <- merge (gene2go, go2name, by = 'GO', all.x = TRUE)

# save data
write.table (go_anno, gene_go, quote = FALSE, row.names = F, sep = "\t")

