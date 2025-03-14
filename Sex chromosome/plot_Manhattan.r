print_help <- function() {
    cat("
This script is used to draw a Manhattan plot.

Usage: Rscript script_name.r <input_file> <output_file>

Arguments:
  <input_file>   The input file containing GWAS results. It should have columns 'SNP', 'CHR', 'BP', 'P'.
  <output_file>  The output file path where the Manhattan plot will be saved in PNG format.

Author: Hangyu Wang
Date: Feb 04 2025
Unit: Southwest University
Contact: wanghyx666@163.com
")
}

args <- commandArgs()
if ("-h" %in% args) {
    print_help()
    q(save = "no")
}

if (length(commandArgs(trailingOnly = TRUE)) < 2) {
    stop("Please provide both input and output file paths. Use -h for help.")
}

input_file <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]

library(tidyverse)
library(ggplot2)

gwasResults <- read.table(input_file, header = T, na.strings = c("NA"))
colnames(gwasResults) <- c("SNP", "CHR", "BP", "P")
chr_len <- gwasResults %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP))
chr_pos <- chr_len  %>%
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)
Snp_pos <- chr_pos %>%
  left_join(gwasResults, ., by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + total)
X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center = ( max(BPcum) + min(BPcum) ) / 2 )
ggplot(Snp_pos, aes(x = BPcum, y = P)) +
  geom_point( aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("#030708", "#bebebe"), 25 )) +
  scale_x_continuous( label = X_axis$CHR, breaks = X_axis$center ) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  labs(title = "", x = "", y = "female/male read depth ratio") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    strip.background = element_rect(fill = "grey90", colour = "black", size = 1),
    strip.text = element_text(size = 20, face = "plain"),
    axis.title = element_text(face = "plain", size = 16, colour = 'black'), 
    axis.text = element_text(face = "plain", size = 14, colour = 'black'), 
    axis.line = element_line(size = 0.5, colour = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
ggsave(filename = output_file, width = 9, height = 3)
