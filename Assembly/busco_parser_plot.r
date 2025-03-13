print_help <- function() {
    cat("
This script generates stacked horizontal bar charts for BUSCO assessment results
while strictly preserving the original species order from input data.

Usage: Rscript busco_plot.R <input_tsv> <output_plot> [options]

Arguments:
  <input_tsv>    BUSCO summary table (TSV format)
  <output_plot>  Output plot file (PDF/PNG/SVG supported)

Options:
  -w <number>    Figure width in inches (default: 10)
  -H <number>    Figure height in inches (default: auto)
  -h             Show this help message

Author: Haoyu Wang
Date: Dec 30  2024
Affiliation: Southwest University
Contact: wanghyx666@163.com
")
}

args <- commandArgs(trailingOnly = TRUE)

if ("-h" %in% args) {
    print_help()
    quit(status = 0)
}

fig_width <- 10
fig_height <- NULL
input_tsv <- NULL
output_file <- NULL

i <- 1
while(i <= length(args)) {
    if(args[i] == "-w") {
        fig_width <- as.numeric(args[i+1])
        i <- i + 2
    } else if(args[i] == "-H") {
        fig_height <- as.numeric(args[i+1])
        i <- i + 2
    } else if(grepl("^-", args[i])) {
        stop("Unknown option: ", args[i])
    } else {
        if(is.null(input_tsv)) {
            input_tsv <- args[i]
        } else if(is.null(output_file)) {
            output_file <- args[i]
        } else {
            stop("Unexpected argument: ", args[i])
        }
        i <- i + 1
    }
}

if(is.null(input_tsv) || is.null(output_file)) {
    stop("Missing required arguments. Usage:\n  Rscript busco_plot.R <input> <output> [-w width] [-H height]")
}

if(!file.exists(input_tsv)) {
    stop("Input file not found: ", input_tsv)
}

required_packages <- c("ggplot2", "tidyr", "dplyr", "forcats", "scales")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if(length(missing_packages) > 0) {
    stop("Missing packages: ", paste(missing_packages, collapse = ", "),
         "\nInstall with: install.packages(c('", paste(missing_packages, collapse = "','"), "'))")
}

library(ggplot2)
library(tidyr)
library(dplyr)
library(forcats)
library(scales)

df <- read.delim(input_tsv, stringsAsFactors = FALSE) %>%
    mutate(species = factor(species, levels = unique(species)))

required_cols <- c("species", "percent_S", "percent_D", "percent_F", "percent_M")
missing_cols <- setdiff(required_cols, colnames(df))
if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

plot_data <- df %>%
    mutate(across(starts_with("percent_"), ~ as.numeric(sub("%", "", .)))) %>%
    select(species, 
           "Single-copy" = percent_S,
           "Duplicated" = percent_D,
           "Fragmented" = percent_F,
           "Missing" = percent_M) %>%
    pivot_longer(
        cols = -species,
        names_to = "Category",
        values_to = "Percentage"
    ) %>%
    mutate(
        Category = factor(Category, 
                         levels = c("Missing", "Fragmented", "Duplicated", "Single-copy"))
    )

n_species <- n_distinct(plot_data$species)
base_height <- 8
per_species_height <- 0.3
auto_height <- base_height + n_species * per_species_height

if(is.null(fig_height)) {
    fig_height <- auto_height
    message("Auto-set height to ", round(fig_height, 1), " inches for ", n_species, " species")
}

color_palette <- c(
    "Single-copy" = "#4E79A7",
    "Duplicated" = "#F28E2B",
    "Fragmented" = "#59A14F",
    "Missing" = "#E15759"
)

smart_wrap <- function(x, width = 25) {
    sapply(x, function(s) {
        if(nchar(s) > width) {
            paste(strwrap(s, width), collapse = "\n")
        } else {
            s
        }
    })
}

p <- ggplot(plot_data, aes(x = Percentage, y = species, fill = Category)) +
    geom_col(position = "stack", width = 0.85) +
    scale_fill_manual(values = color_palette) +
    labs(
        title = "BUSCO Assessment Results",
        x = "Percentage (%)",
        y = NULL,
        fill = "Category"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.y = element_text(size = rel(0.9), margin = margin(r = 10), hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 15)),
        legend.position = "right",
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90"),
        plot.margin = margin(1, 2, 1, 1.5, "cm")
    ) +
    scale_x_continuous(labels = label_percent(scale = 1)) +
    scale_y_discrete(labels = smart_wrap, limits = rev) +
    coord_cartesian(clip = "off")

ggsave(
    filename = output_file,
    plot = p,
    width = fig_width,
    height = fig_height,
    dpi = 300
)

message("\nSuccessfully generated: ", output_file)
message("Dimensions: ", fig_width, "\" x ", fig_height, "\"")
