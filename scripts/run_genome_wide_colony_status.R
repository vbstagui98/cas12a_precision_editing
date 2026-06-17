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
      "  Rscript scripts/run_genome_wide_colony_status.R \\",
      "    --output-dir results/genome_wide_colony_status",
      "",
      "Optional path overrides:",
      "  --picked-colonies-round-1 path/to/picked_colonies_round_1.csv",
      "  --picked-colonies-round-2 path/to/picked_colonies_round_2.csv",
      "  --target-coverage-round-1 path/to/target_coverage_round_1.tsv",
      "  --target-coverage-round-2 path/to/target_coverage_round_2.tsv",
      "  --mean-coverage-round-1 path/to/mean_coverage_round_1.csv",
      "  --mean-coverage-round-2 path/to/mean_coverage_round_2.csv",
      "  --normalized-variants-round-1 path/to/variants_norm_round_1.csv",
      "  --normalized-variants-round-2 path/to/variants_norm_round_2.csv",
      "  --design-annotations path/to/designs_812_sv_my_variants_annotated.tsv",
      "  --deepcpf1-scores path/to/input_deepcpf1_20260319.scored.csv",
      "  --dr-scores path/to/guides_dr_scores_20260325.tsv",
      "  --normalized-indels path/to/normalized_indel_designs.vcf",
      "",
      "Rules:",
      "  Intended editing: intended HDR allele has AO >= 2.",
      "  Non-target variants: AF >= 50 and DP >= 4.",
      "  No variant quality threshold is applied.",
      "",
      "Outputs:",
      "  genome_wide_colony_editing_status.csv",
      "  genome_wide_design_editing_status.csv",
      "  genome_wide_non_target_af50_dp4_variants.csv",
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
script_path <- normalizePath(sub("^--file=", "", file_arg[1] %||% "scripts/run_genome_wide_colony_status.R"), mustWork = FALSE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(project_root, "config", "paths.R"))
source(file.path(project_root, "R", "genome_wide_colony_status.R"))

paths <- get_analysis_paths()

default_inputs <- list(
  picked_colonies_round_1 = file.path(paths$genome_wide_root, "REDI", "picked_colonies_20250613_annotated.csv"),
  picked_colonies_round_2 = file.path(
    paths$genome_wide_root,
    "REDI",
    "second_REDI_20250917",
    "scritps_20251216",
    "results",
    "colonies_with_BC1_20260212.csv"
  ),
  target_coverage_round_1 = file.path(paths$genome_wide_root, "Tn5", "Analysis", "stats", "target_basecov_20260222.tsv"),
  target_coverage_round_2 = file.path(
    paths$genome_wide_root,
    "Tn5",
    "second_round_20260130",
    "read_stats",
    "target_coverage_20260213.tsv"
  ),
  mean_coverage_round_1 = file.path(
    paths$genome_wide_root,
    "Tn5",
    "second_round_20260130",
    "read_stats",
    "mean_cov_20260223.csv"
  ),
  mean_coverage_round_2 = file.path(
    paths$genome_wide_root,
    "Tn5",
    "second_round_20260130",
    "read_stats",
    "mean_cov_20260223_second.csv"
  ),
  normalized_variants_round_1 = file.path(
    paths$genome_wide_root,
    "Tn5",
    "Analysis",
    "results",
    "variants_norm_WGS_20260222.csv"
  ),
  normalized_variants_round_2 = file.path(
    paths$genome_wide_root,
    "Tn5",
    "second_round_20260130",
    "results",
    "variants_norm_WGS_20260219.csv"
  ),
  design_annotations = file.path(project_root, "scripts", "designs_812_sv_my_variants_annotated.tsv"),
  deepcpf1_scores = file.path(project_root, "input_deepcpf1_20260319.scored.csv"),
  dr_scores = file.path(
    paths$genome_wide_root,
    "Tn5",
    "Analysis",
    "annotations_guides",
    "guides_dr_scores_20260325.tsv"
  ),
  normalized_indels = file.path(
    paths$genome_wide_root,
    "Tn5",
    "second_round_20260130",
    "annotations_ref",
    "normalized_indel_designs_20260223.vcf"
  )
)

input_paths <- list(
  picked_colonies_round_1 = opts$picked_colonies_round_1 %||% default_inputs$picked_colonies_round_1,
  picked_colonies_round_2 = opts$picked_colonies_round_2 %||% default_inputs$picked_colonies_round_2,
  target_coverage_round_1 = opts$target_coverage_round_1 %||% default_inputs$target_coverage_round_1,
  target_coverage_round_2 = opts$target_coverage_round_2 %||% default_inputs$target_coverage_round_2,
  mean_coverage_round_1 = opts$mean_coverage_round_1 %||% default_inputs$mean_coverage_round_1,
  mean_coverage_round_2 = opts$mean_coverage_round_2 %||% default_inputs$mean_coverage_round_2,
  normalized_variants_round_1 = opts$normalized_variants_round_1 %||% default_inputs$normalized_variants_round_1,
  normalized_variants_round_2 = opts$normalized_variants_round_2 %||% default_inputs$normalized_variants_round_2,
  design_annotations = opts$design_annotations %||% default_inputs$design_annotations,
  deepcpf1_scores = opts$deepcpf1_scores %||% default_inputs$deepcpf1_scores,
  dr_scores = opts$dr_scores %||% default_inputs$dr_scores,
  normalized_indels = opts$normalized_indels %||% default_inputs$normalized_indels
)

outputs <- run_genome_wide_colony_status(input_paths, opts$output_dir)

colony_status <- outputs$colony_status
design_status <- outputs$design_status

message("Colony rows: ", nrow(colony_status))
message("Assayed colonies: ", sum(colony_status$assayed_colony, na.rm = TRUE))
message("Edited colonies: ", sum(colony_status$edited_colony, na.rm = TRUE))
message("Design rows: ", nrow(design_status))
message("Assayed designs: ", sum(design_status$assayed_design, na.rm = TRUE))
message("Edited designs: ", sum(design_status$edited_design, na.rm = TRUE))
message("Wrote genome-wide outputs to: ", opts$output_dir)
