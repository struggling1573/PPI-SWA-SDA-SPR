print_help <- function() {
    cat("
This script is used to plot a percentage stacked bar chart based on allele count data.

Usage: Rscript plot.r <input> <output> <sample> 

Arguments:
  <input>    The input file containing allele count data with a header.
  <output>   The output file path for saving the generated plot.
  <sample>   The sample name used in the plot title.

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
    stop("Please provide both input and output file paths and sample name. Use -h for help.")
}

input <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]
sample <- commandArgs(trailingOnly = TRUE)[3]

library(ggplot2)
library(dplyr)
library(tidyr)

data <- read.table(input, header = T)
data_long <- data %>%
  gather(key = "Type", value = "Value", -CHR)
data_filtered <- data_long %>% filter(.[[2]]!= "AlleNumber" & .[[2]]!= "Orthe")
data_filtered <- data_filtered %>%
  group_by(CHR) %>%
  mutate(Percentage = Value / sum(Value) * 100)
plot_title <- paste("Genotypic statistics in ", sample, sep = "")
ggplot(data_filtered, aes(x=CHR, y=Percentage, fill=Type)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("Hexasomic"="#289584", "Tetrasomic"="#ccb280", "Disomic"="#3d3a6b", "NoAlle"="#65b961", "Others"="#faebd7")) +
  labs(y="Proportion(%)", x="", title = plot_title, fill="Type") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill="grey90", colour="black", size=1),
    strip.text = element_text(size=12, face="bold"),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(output, width = 5, height = 3)    
