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
      "  Rscript scripts/run_amplicon_variant_efficiency.R \\",
      "    --variants path/to/genomic_ampli_enas_fn_harmonized_20260318.csv \\",
      "    --output results/amplicon_panel_efficiency.csv \\",
      "    [--guide-features results/panel_sequence_manifest/combined_sequence_manifest.csv]",
      "",
      "Input must be a long amplicon variant table with at least:",
      "  MUTATION, frc_alt, frc_ref",
      "",
      "Output:",
      "  A long efficiency table with HDR, reference, and non-HDR percentages.",
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

required <- c("variants", "output")
missing <- setdiff(required, names(opts))
if (length(missing) > 0) {
  usage()
  stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[1] %||% "scripts/run_amplicon_variant_efficiency.R"), mustWork = FALSE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(project_root, "R", "amplicon_variant_efficiency.R"))

efficiency <- amplicon_efficiency_from_long_table(
  opts$variants,
  guide_features_path = opts$guide_features %||% NULL
)

write_csv_mkdir(efficiency, opts$output)

message("Rows written: ", nrow(efficiency))
message("Output written to: ", opts$output)
