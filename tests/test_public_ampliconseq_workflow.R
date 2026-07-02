#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", file_arg[1])
project_root <- dirname(dirname(normalizePath(script_file)))

source(file.path(project_root, "R", "amplicon_variant_efficiency.R"))

paths <- list(
  fn_donor_variants = file.path(
    project_root,
    "data/publication_inputs/ampliconseq_fn_donor_variants.csv"
  ),
  enas_donor_variants = file.path(
    project_root,
    "data/publication_inputs/ampliconseq_enas_donor_short_guide_variants.csv"
  ),
  short_guide_variants = file.path(
    project_root,
    "data/publication_inputs/ampliconseq_enas_donor_short_guide_variants.csv"
  )
)

output_dir <- tempfile("ampliconseq_test_")
results <- prepare_donor_and_guide_length_outputs(paths, output_dir)

stopifnot(
  nrow(results$fn_donor) == 182,
  nrow(results$enas_donor) == 180,
  nrow(results$short_guide) == 124,
  nrow(results$figure_3) == 124,
  isTRUE(all.equal(sum(results$fn_donor$HDR_pct), 11031.18423882289)),
  isTRUE(all.equal(sum(results$enas_donor$HDR_pct), 11952.914159520275)),
  isTRUE(all.equal(sum(results$short_guide$HDR_pct), 7095.541466815113)),
  isTRUE(all.equal(sum(results$figure_3$efficiency_pct), 9128.668892476493)),
  all(file.exists(file.path(
    output_dir,
    c(
      "fn_donor_position_hdr_from_variants.csv",
      "enas_donor_position_hdr_from_variants.csv",
      "donor_position_hdr_from_variants.csv",
      "short_guide_hdr_from_variants.csv",
      "figure_3_ampliconseq_efficiency.csv"
    )
  )))
)

message("Public ampliconseq workflow regression test passed.")
