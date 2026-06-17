#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) {
    y
  } else {
    x
  }
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/run_donor_guide_variant_efficiency.R --output-dir results/donor_guide_variant_efficiency",
      "",
      "Optional path overrides:",
      "  --fn-donor-variants path/to/20250312_ys88_donor_mut/data/raw/master_df_filtered.csv",
      "  --enas-donor-variants path/to/20250611_amplicons_ys85_misc/data/raw/master_df_filtered.csv",
      "  --short-guide-variants path/to/20250611_amplicons_shorterguide/data/raw/master_df_filtered.csv",
      "",
      "Outputs:",
      "  fn_donor_position_hdr_from_variants.csv",
      "  enas_donor_position_hdr_from_variants.csv",
      "  donor_position_hdr_from_variants.csv",
      "  short_guide_hdr_from_variants.csv",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  out <- list()
  i <- 1

  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }

    key <- sub("^--", "", key)
    value <- args[[i + 1]]
    if (is.na(value) || startsWith(value, "--")) {
      stop("Missing value for --", key, call. = FALSE)
    }

    out[[gsub("-", "_", key)]] <- value
    i <- i + 2
  }

  out
}

if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
  usage()
  quit(status = 0)
}

opts <- parse_args(args)

if (is.null(opts$output_dir)) {
  usage()
  stop("Missing required argument: --output-dir", call. = FALSE)
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[1] %||% "scripts/run_donor_guide_variant_efficiency.R"), mustWork = FALSE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(project_root, "config", "paths.R"))
source(file.path(project_root, "R", "amplicon_variant_efficiency.R"))

paths <- get_analysis_paths()

default_inputs <- list(
  fn_donor_variants = file.path(paths$presentations_amplicon_root, "20250312_ys88_donor_mut", "data", "raw", "master_df_filtered.csv"),
  enas_donor_variants = file.path(paths$presentations_amplicon_root, "20250611_amplicons_ys85_misc", "data", "raw", "master_df_filtered.csv"),
  short_guide_variants = file.path(paths$presentations_amplicon_root, "20250611_amplicons_shorterguide", "data", "raw", "master_df_filtered.csv")
)

input_paths <- list(
  fn_donor_variants = opts$fn_donor_variants %||% default_inputs$fn_donor_variants,
  enas_donor_variants = opts$enas_donor_variants %||% default_inputs$enas_donor_variants,
  short_guide_variants = opts$short_guide_variants %||% default_inputs$short_guide_variants
)

outputs <- prepare_donor_and_guide_length_outputs(input_paths, opts$output_dir)

message("Fn donor rows: ", nrow(outputs$fn_donor))
message("enAs donor rows: ", nrow(outputs$enas_donor))
message("Short-guide rows: ", nrow(outputs$short_guide))
message("Wrote donor/guide outputs to: ", opts$output_dir)
