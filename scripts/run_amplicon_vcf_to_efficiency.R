#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) y else x
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/run_endogenous_amplicon_vcf_to_efficiency.R \\",
      "    --variants path/to/concatenated_freebayes_variants.csv \\",
      "    --designs path/to/intended_loci.csv \\",
      "    [--guide-features path/to/endogenous_guide_features.csv] \\",
      "    --output-dir results/amplicon_efficiency",
      "",
      "Compatibility alias:",
      "  Rscript scripts/run_amplicon_vcf_to_efficiency.R ...",
      "",
      "Required long-table columns (common aliases are accepted):",
      "  Sample, guide/crRNA ID, CHROM, POS, REF, ALT, TYPE, AO, RO, DP",
      "  frc_alt and frc_ref must not be precomputed; they are calculated here",
      "",
      "Design assignment columns:",
      "  crRNA_id/Guide plus either match/intended_match, or exact",
      "  intended_chromosome, intended_position, and alternate_allele",
      "",
      "Rules:",
      "  intended edit detected: AO >= 2",
      "  editing-window frequencies: count / (RO + sum(AO at the intended locus))",
      "  QUAL is retained when present but is never used as a filter",
      sep = "\n"
    )
  )
}

parse_args <- function(values) {
  out <- list()
  i <- 1

  while (i <= length(values)) {
    key <- values[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    if (i == length(values) || startsWith(values[[i + 1]], "--")) {
      stop("Missing value for ", key, call. = FALSE)
    }

    out[[gsub("-", "_", sub("^--", "", key))]] <- values[[i + 1]]
    i <- i + 2
  }

  out
}

if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
  usage()
  quit(status = 0)
}

opts <- parse_args(args)
required <- c("variants", "designs", "output_dir")
missing <- setdiff(required, names(opts))
if (length(missing) > 0) {
  usage()
  stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", file_arg[1] %||% "scripts/run_amplicon_vcf_to_efficiency.R"),
  mustWork = FALSE
)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
source(file.path(project_root, "R", "amplicon_variant_efficiency.R"))

intended_ao_min <- as.numeric(opts$intended_ao_min %||% 2)
min_dp <- as.numeric(opts$min_dp %||% 0)

designs <- read_variant_table(opts$designs)
if (!is.null(opts$guide_features)) {
  designs <- merge_endogenous_design_features(
    designs,
    read_variant_table(opts$guide_features)
  )
}

annotated <- annotate_endogenous_amplicon_variants(
  long_variant_table = opts$variants,
  design_table = designs,
  intended_ao_min = intended_ao_min,
  min_dp = min_dp
)

editing_window <- build_endogenous_editing_window(annotated)
efficiency <- endogenous_amplicon_efficiency_table(editing_window)
hdr_efficiency <- endogenous_hdr_efficiency_table(annotated)

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
annotated_path <- file.path(
  opts$output_dir,
  "amplicon_variants_annotated.csv"
)
efficiency_path <- file.path(
  opts$output_dir,
  "amplicon_efficiency.csv"
)
editing_window_path <- file.path(
  opts$output_dir,
  "amplicon_editing_window.csv"
)
hdr_efficiency_path <- file.path(
  opts$output_dir,
  "amplicon_hdr_efficiency.csv"
)

readr::write_csv(annotated, annotated_path)
readr::write_csv(editing_window, editing_window_path)
readr::write_csv(efficiency, efficiency_path)
readr::write_csv(hdr_efficiency, hdr_efficiency_path)

message("Annotated variant rows written: ", nrow(annotated))
message("Annotated variants: ", annotated_path)
message("Editing window: ", editing_window_path)
message("Efficiency table: ", efficiency_path)
message("Sample-level HDR efficiency: ", hdr_efficiency_path)
