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
      "  Rscript scripts/run_fcs_efficiency.R \\",
      "    --fcs-dir data/raw/FACS_DATA \\",
      "    --platemap data/raw/platemap.csv \\",
      "    --output results/fcs_efficiency \\",
      "    [--metadata data/raw/metadata.csv] \\",
      "    [--guide-features data/raw/guides_donor_selected_info.csv] \\",
      "    [--colony-counts data/raw/colony_counts.csv] \\",
      "    [--assay-col assay] [--min-events 500]",
      "",
      "Outputs:",
      "  pop_stats.csv       Raw GFP-gate population frequencies",
      "  event_counts.csv    Event counts used for the minimum-event filter",
      "  fcs_efficiency.csv  Annotated GFP-positive and editing-efficiency table",
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

required <- c("fcs_dir", "platemap", "output")
missing <- setdiff(required, names(opts))
if (length(missing) > 0) {
  usage()
  stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[1] %||% "scripts/run_fcs_efficiency.R"), mustWork = FALSE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(project_root, "R", "fcs_efficiency.R"))

run_fcs_efficiency(
  fcs_dir = opts$fcs_dir,
  platemap_path = opts$platemap,
  output_dir = opts$output,
  metadata_path = opts$metadata %||% NULL,
  guide_features_path = opts$guide_features %||% NULL,
  colony_counts_path = opts$colony_counts %||% NULL,
  min_events = as.integer(opts$min_events %||% 500),
  assay_col = opts$assay_col %||% "assay",
  control_well = opts$control_well %||% "A09"
)

message("Wrote FCS efficiency outputs to: ", opts$output)
