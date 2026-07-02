#!/usr/bin/env Rscript

## Calculate HDR efficiencies for the donor-position and shorter-guide assays.

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
    "  Rscript scripts/run_donor_guide_variant_efficiency.R \\",
    "    --fn-donor-variants path/to/fn_donor_variants.csv \\",
    "    --enas-donor-variants path/to/enas_donor_variants.csv \\",
    "    --short-guide-variants path/to/short_guide_variants.csv \\",
    "    --output-dir results/donor_guide_variant_efficiency",
    sep = "\n"
  ))
  quit(status = 0)
}

opts <- read_args(args)
required <- c(
  "fn_donor_variants",
  "enas_donor_variants",
  "short_guide_variants",
  "output_dir"
)
missing <- setdiff(required, names(opts))
if (length(missing) > 0) {
  stop("Missing arguments: ", paste(missing, collapse = ", "), call. = FALSE)
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", file_arg[1])
project_root <- dirname(dirname(normalizePath(script_file)))
source(file.path(project_root, "R", "amplicon_variant_efficiency.R"))

input_paths <- list(
  fn_donor_variants = opts$fn_donor_variants,
  enas_donor_variants = opts$enas_donor_variants,
  short_guide_variants = opts$short_guide_variants
)

outputs <- prepare_donor_and_guide_length_outputs(input_paths, opts$output_dir)

message("Fn donor rows: ", nrow(outputs$fn_donor))
message("enAs donor rows: ", nrow(outputs$enas_donor))
message("Short-guide rows: ", nrow(outputs$short_guide))
