print_help <- function() {
    cat("
This script is used to get intersection between several files.

Usage: Rscript venn_intersect.r <file1> <file2> ...

Arguments:
  <file1>    file1
  <file2>    file2
  ....

Author: Hangyu Wang
Date: Jan 06 2024
Unit: Southwest University
Contact: wanghyx666@163.com
")
}

args <- commandArgs()
if ("-h" %in% args) {
    print_help()
    q(save = "no")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("At least two files are needed for comparison.\nPlease provide both input and output file paths. Use -h for help.", call. = FALSE)
}
read_file_to_vector <- function(file_path) {
  read.csv(file_path, header = FALSE, stringsAsFactors = FALSE)$V1
}
files_content <- lapply(args, read_file_to_vector)
common_to_all <- Reduce(intersect, files_content)
unique_results <- list()
pairwise_common_results <- list()
for (i in seq_along(files_content)) {
  current_file <- files_content[[i]]
  other_files <- files_content[-i]
  unique_to_current <- setdiff(current_file, unlist(other_files))
  unique_results[[args[i]]] <- unique_to_current
  for (j in seq_along(other_files)) {
    if (i != j + i - 1) {
      common_with_other <- setdiff(intersect(current_file, other_files[[j]]), common_to_all)
      pairwise_common_key <- paste0("common_", args[i], "_", args[j + i - 1])
      if (!exists(pairwise_common_key, where = pairwise_common_results)) {
        pairwise_common_results[[pairwise_common_key]] <- common_with_other
      } else {
        existing_common <- pairwise_common_results[[pairwise_common_key]]
        combined_common <- union(existing_common, common_with_other)
        pairwise_common_results[[pairwise_common_key]] <- combined_common
      }
    }
  }
}

for (file_name in names(unique_results)) {
  write.table(unique_results[[file_name]], file = paste0("vennOut_unique_to_", file_name, ".csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}
for (common_key in names(pairwise_common_results)) {
  write.table(pairwise_common_results[[common_key]], file = paste0("vennOut_", common_key, ".csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}
write.table(common_to_all, file = "vennOut_common_to_all.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

