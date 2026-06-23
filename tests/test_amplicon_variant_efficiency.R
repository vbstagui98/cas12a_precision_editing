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

annotated <- annotate_endogenous_amplicon_variants(variants, designs)

stopifnot(
  annotated$variant_annotation[[1]] == "HDR",
  annotated$variant_annotation[[2]] == "HDR_below_read_threshold",
  annotated$variant_annotation[[3]] == "other_variant",
  annotated$variant_annotation[[4]] == "other_variant",
  annotated$variant_annotation[[5]] == "REF"
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

raw_annotated <- annotate_endogenous_amplicon_variants(raw_counts, designs)
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
  missing_ref$efficiency_pct[[1]] == 100
)

publication_table <- tibble::tibble(
  cas_variant = "FnCas12a",
  crRNA_promoter = "RPR1",
  direct_repeat = "As_DR",
  donor_recruitment = "LexA-FHA",
  crRNA_id = "G03",
  replicate = "D1",
  timepoint = "T1",
  generations = 0,
  editing_class = "HDR",
  efficiency_pct = 15,
  reference_pct = 80
)

features <- tibble::tibble(
  crRNA_id = "G03",
  target = "ADE2_3_donor_2",
  variant_position_annotation = "chrXV:565450",
  crRNA_sequence = "CCGGTTGTGGTATATTTGGTGTG",
  donor_sequence = "ACGT",
  mutation_type = "1 nt substitution",
  alternate_allele = "A",
  pam = "TTTC",
  pre_pam_nt = "T",
  distance_from_pam = 14,
  strand = "+",
  crRNA_gc_pct = 47.83,
  deepcpf1_score = 0.4609433
)

features_path <- tempfile(fileext = ".csv")
readr::write_csv(features, features_path)
efficiency <- amplicon_efficiency_from_long_table(
  publication_table,
  guide_features_path = features_path
)
unlink(features_path)

stopifnot(
  efficiency$guide_id[[1]] == "G03",
  efficiency$mutation_class[[1]] == "HDR",
  efficiency$allele_frequency_pct[[1]] == 15,
  efficiency$ref_pct[[1]] == 80,
  efficiency$non_hdr_pct[[1]] == 5,
  efficiency$deepcpf1_score[[1]] == 0.4609433
)

cat("amplicon variant efficiency tests passed\n")
