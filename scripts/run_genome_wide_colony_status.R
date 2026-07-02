#!/usr/bin/env Rscript

## Assign intended and non-target editing status to genome-wide colonies.

args <- commandArgs(trailingOnly = TRUE)

read_args <- function(x) {
  values <- list()
  for (i in seq(1, length(x), by = 2)) {
    if (i == length(x) || !startsWith(x[i], "--")) {
      stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
    }
    name <- gsub("-", "_", sub("^--", "", x[i]))
    values[[name]] <- x[i + 1]
  }
  values
}

if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
  cat(paste(
    "Usage:",
    "  Rscript scripts/run_genome_wide_colony_status.R \\",
    "    --picked-colonies-round-1 picked_round1.csv \\",
    "    --picked-colonies-round-2 picked_round2.csv \\",
    "    --target-coverage-round-1 target_coverage_round1.tsv \\",
    "    --target-coverage-round-2 target_coverage_round2.tsv \\",
    "    --mean-coverage-round-1 mean_coverage_round1.csv \\",
    "    --mean-coverage-round-2 mean_coverage_round2.csv \\",
    "    --normalized-variants-round-1 variants_round1.csv \\",
    "    --normalized-variants-round-2 variants_round2.csv \\",
    "    --design-annotations design_annotations.tsv \\",
    "    --deepcpf1-scores deepcpf1_scores.csv \\",
    "    --dr-scores direct_repeat_scores.tsv \\",
    "    --normalized-indels normalized_indel_designs.vcf \\",
    "    --output-dir results/genome_wide_colony_status",
    "",
    "Intended edits require AO >= 2.",
    "Non-target variants require AF >= 50 and DP >= 4.",
    "QUAL is retained but is not used as a filter.",
    sep = "\n"
  ))
  quit(status = 0)
}

opts <- read_args(args)
input_names <- c(
  "picked_colonies_round_1",
  "picked_colonies_round_2",
  "target_coverage_round_1",
  "target_coverage_round_2",
  "mean_coverage_round_1",
  "mean_coverage_round_2",
  "normalized_variants_round_1",
  "normalized_variants_round_2",
  "design_annotations",
  "deepcpf1_scores",
  "dr_scores",
  "normalized_indels"
)
missing <- setdiff(c(input_names, "output_dir"), names(opts))
if (length(missing) > 0) {
  stop("Missing arguments: ", paste(missing, collapse = ", "), call. = FALSE)
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", file_arg[1])
project_root <- dirname(dirname(normalizePath(script_file)))
source(file.path(project_root, "R", "genome_wide_colony_status.R"))

input_paths <- opts[input_names]
outputs <- run_genome_wide_colony_status(input_paths, opts$output_dir)

message("Colony rows: ", nrow(outputs$colony_status))
message("Assayed colonies: ", sum(outputs$colony_status$assayed_colony, na.rm = TRUE))
message("Edited colonies: ", sum(outputs$colony_status$edited_colony, na.rm = TRUE))
message("Design rows: ", nrow(outputs$design_status))
