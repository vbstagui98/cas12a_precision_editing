#!/usr/bin/env Rscript

source("R/amplicon_variant_efficiency.R")

designs <- tibble::tibble(
  promoter = "RPR1",
  Guide = "G03",
  match = "RPR1_G03_chrXV_565436_A"
)

variants <- tibble::tibble(
  promoter = rep("RPR1", 5),
  Guide = rep("G03", 5),
  Sample = paste0("sample_", seq_len(5)),
  CHROM = rep("chrXV", 5),
  POS = c(565436, 565436, 565500, 565501, 565436),
  REF = c("T", "T", "G", "G", "T"),
  ALT = c("A", "A", "C", "C", "T"),
  TYPE = c("snp", "snp", "snp", "snp", "REF"),
  pos_mismatch = c(565436, 565436, 565500, 565501, 565436),
  mismatches = c("A", "A", "C", "C", "A"),
  AO = c(2, 1, 2, 2, 2),
  RO = c(2, 3, 2, 1, 2),
  DP = c(4, 4, 4, 3, 4)
)

annotated <- annotate_endogenous_amplicon_variants(
  variants,
  designs,
  min_dp = 0
)

stopifnot(
  annotated$variant_annotation[[1]] == "HDR",
  annotated$variant_annotation[[2]] == "HDR",
  annotated$variant_annotation[[3]] == "other_variant",
  annotated$variant_annotation[[4]] == "other_variant",
  annotated$variant_annotation[[5]] == "REF",
  annotated$MUTATION[[1]] == "HDR",
  annotated$MUTATION[[2]] == "HDR",
  annotated$MUTATION[[3]] == "snp",
  annotated$MUTATION[[5]] == "REF",
  !any(annotated$variant_annotation == "unintended_AF50_DP4", na.rm = TRUE),
  !any(annotated$variant_annotation == "HDR_below_read_threshold", na.rm = TRUE)
)

raw_counts <- tibble::tibble(
  promoter = rep("RPR1", 4),
  Guide = rep("G03", 4),
  Sample = c("edited", "edited", "edited", "not_detected"),
  CHROM = rep("chrXV", 4),
  POS = c(565436, 565436, 565500, 565500),
  REF = c("T", "T", "G", "G"),
  ALT = c("A", "G", "C", "C"),
  TYPE = rep("snp", 4),
  AO = c(20, 10, 80, 10),
  RO = c(70, 70, 20, 90),
  DP = rep(100, 4),
  frc_alt = rep(999, 4),
  frc_ref = rep(999, 4)
)

raw_annotated <- annotate_endogenous_amplicon_variants(
  raw_counts,
  designs,
  min_dp = 0
)
editing_window <- build_endogenous_editing_window(raw_annotated)
raw_efficiency <- endogenous_amplicon_efficiency_table(editing_window)

edited_hdr <- raw_efficiency[
  raw_efficiency$sample == "edited" & raw_efficiency$editing_class == "HDR",
]
edited_ref <- raw_efficiency[
  raw_efficiency$sample == "edited" & raw_efficiency$editing_class == "REF",
]
missing_hdr <- raw_efficiency[
  raw_efficiency$sample == "not_detected" & raw_efficiency$editing_class == "HDR",
]
missing_ref <- raw_efficiency[
  raw_efficiency$sample == "not_detected" & raw_efficiency$editing_class == "REF",
]

stopifnot(
  edited_hdr$efficiency_pct[[1]] == 20,
  edited_ref$efficiency_pct[[1]] == 70,
  edited_hdr$non_hdr_pct[[1]] == 10,
  missing_hdr$efficiency_pct[[1]] == 0,
  missing_ref$efficiency_pct[[1]] == 100,
  sum(editing_window$Sample == "edited" & editing_window$editing_class == "HDR") == 1,
  sum(editing_window$Sample == "edited" & editing_window$editing_class == "snp") == 1,
  sum(editing_window$Sample == "edited" & editing_window$editing_class == "REF") == 1,
  all(editing_window$editing_window_pct[editing_window$Sample == "edited"] != 999)
)

cat("amplicon variant efficiency tests passed\n")
