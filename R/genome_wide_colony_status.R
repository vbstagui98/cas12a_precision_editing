suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) {
    y
  } else {
    x
  }
}

clean_table_names <- function(x) {
  x <- gsub("([a-z0-9])([A-Z])", "\1_\2", x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("%", "pct", x, fixed = TRUE)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)

  empty_idx <- which(x == "")
  if (length(empty_idx) > 0) {
    x[empty_idx] <- paste0("col_", empty_idx)
  }

  make.unique(x, sep = "_")
}

read_csv_clean <- function(path, ...) {
  data <- readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal", ...)
  names(data) <- clean_table_names(names(data))
  tibble::as_tibble(data)
}

read_tsv_clean <- function(path, ...) {
  data <- readr::read_tsv(path, show_col_types = FALSE, name_repair = "minimal", ...)
  names(data) <- clean_table_names(names(data))
  tibble::as_tibble(data)
}

write_csv_mkdir <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(data, path)
  invisible(path)
}

check_required_files <- function(paths) {
  missing_paths <- names(paths)[!file.exists(unlist(paths, use.names = FALSE))]

  if (length(missing_paths) > 0) {
    stop(
      "Missing required input files: ",
      paste(missing_paths, collapse = ", "),
      call. = FALSE
    )
  }
}

assert_columns <- function(data, required_columns, data_name) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      data_name,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

validate_unique_keys <- function(data, key_columns, data_name) {
  duplicated_keys <- data %>%
    dplyr::count(dplyr::across(dplyr::all_of(key_columns)), name = "n") %>%
    dplyr::filter(.data$n > 1)

  if (nrow(duplicated_keys) > 0) {
    stop(
      data_name,
      " contains duplicated keys for: ",
      paste(key_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(data)
}

parse_number_safe <- function(x) {
  readr::parse_number(as.character(x), locale = readr::locale(decimal_mark = "."))
}

as_logical_flag <- function(x) {
  value <- tolower(as.character(x))
  dplyr::case_when(
    value %in% c("true", "t", "yes", "y", "1") ~ TRUE,
    value %in% c("false", "f", "no", "n", "0") ~ FALSE,
    TRUE ~ NA
  )
}

gc_fraction <- function(sequence) {
  sequence <- toupper(as.character(sequence))
  base_count <- nchar(sequence)
  gc_count <- stringr::str_count(sequence, "[GC]")

  dplyr::if_else(base_count > 0, gc_count / base_count, NA_real_)
}

gc_fraction_bin <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x < 0.30 ~ "<30%",
    x < 0.40 ~ "30-39%",
    x < 0.50 ~ "40-49%",
    x < 0.60 ~ "50-59%",
    TRUE ~ ">=60%"
  )
}

read_genome_wide_variant_calls <- function(path) {
  read_csv_clean(path) %>%
    dplyr::select(-dplyr::any_of(c("x", "x_1", "col_1")))
}

read_normalized_indel_designs <- function(path) {
  data <- readr::read_tsv(
    path,
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE,
    name_repair = "unique_quiet"
  )

  names(data) <- paste0("x", seq_len(ncol(data)))

  data %>%
    dplyr::transmute(
      uid = .data$x3,
      chr = .data$x1,
      pos = as.integer(.data$x2),
      ref = .data$x4,
      alt = .data$x5
    )
}

read_dr_scores <- function(path) {
  dr_scores <- read_tsv_clean(path) %>%
    dplyr::transmute(
      uid = .data$guide_id,
      dr_guide_seq_input = .data$guide_seq_input,
      dr_guide_seq_rna = .data$guide_seq_rna,
      dr_mfe_kcal_mol = as.numeric(.data$mfe_kcal_mol),
      dr_structure = .data$structure,
      dr_len = as.integer(.data$dr_len),
      dr_nt_paired_to_dr = as.integer(.data$dr_nt_paired_to_dr),
      dr_nt_paired_to_spacer = as.integer(.data$dr_nt_paired_to_spacer),
      dr_nt_unpaired = as.integer(.data$dr_nt_unpaired),
      dr_perturbed = as_logical_flag(.data$dr_perturbed),
      dr_score = as.numeric(.data$score)
    )

  validate_unique_keys(dr_scores, "uid", "guide direct-repeat scores")
  dr_scores
}

load_picked_colonies <- function(picked_colonies_round_1, picked_colonies_round_2) {
  round_1 <- read_csv_clean(picked_colonies_round_1) %>%
    dplyr::transmute(
      uid = .data$design,
      bc1 = .data$bc1,
      colony_id = .data$colony_id,
      round = "round_1"
    )

  round_2 <- read_csv_clean(picked_colonies_round_2) %>%
    dplyr::transmute(
      uid = .data$design,
      bc1 = .data$bc1,
      colony_id = .data$colony_id,
      round = "round_2"
    )

  picked <- dplyr::bind_rows(round_1, round_2)
  validate_unique_keys(picked, c("round", "colony_id"), "picked colonies")
  picked
}

load_target_coverage <- function(target_coverage_round_1, target_coverage_round_2) {
  round_1 <- read_tsv_clean(target_coverage_round_1) %>%
    dplyr::mutate(
      colony_id = stringr::str_replace(.data$sample, ".*VB_(.*?)_S.*", "\1"),
      round = "round_1"
    )

  round_2 <- read_tsv_clean(target_coverage_round_2) %>%
    dplyr::mutate(
      colony_id = stringr::str_replace(.data$sample, ".*VB_(.*?)_S.*", "\1"),
      round = "round_2"
    )

  coverage <- dplyr::bind_rows(round_1, round_2) %>%
    dplyr::transmute(
      round = .data$round,
      colony_id = .data$colony_id,
      sample = .data$sample,
      target_basecov = as.numeric(.data$coverage),
      coverage_chr = .data$chrom,
      coverage_pos = as.integer(.data$pos)
    )

  validate_unique_keys(coverage, c("round", "colony_id"), "target coverage")
  coverage
}

load_mean_coverage <- function(mean_coverage_round_1, mean_coverage_round_2) {
  mean_cov <- dplyr::bind_rows(
    read_csv_clean(mean_coverage_round_1),
    read_csv_clean(mean_coverage_round_2)
  ) %>%
    dplyr::mutate(sample = stringr::str_remove(.data$file, "\.bincov\.txt$")) %>%
    dplyr::transmute(
      sample = .data$sample,
      mean_cov = as.numeric(.data$mean_cov)
    )

  validate_unique_keys(mean_cov, "sample", "mean coverage")
  mean_cov
}

load_design_features <- function(design_annotations, deepcpf1_scores, dr_scores) {
  designs <- read_tsv_clean(design_annotations)
  assert_columns(
    designs,
    c(
      "uid", "guide_id", "guide_seq", "guide_nopam", "donor_wt",
      "pam", "gc_content", "max_a", "max_t", "num_t",
      "pam_dist_var_start", "pam_dist_var_end", "vardiff",
      "chr", "pos", "ref", "alt", "in_sv_hotspot",
      "exact_high_score", "hotspot_high_score"
    ),
    "design_annotations"
  )

  designs <- designs %>%
    dplyr::transmute(
      uid = .data$uid,
      guide_id = .data$guide_id,
      guide_seq = .data$guide_seq,
      guide_nopam = .data$guide_nopam,
      donor_wt = .data$donor_wt,
      pam = .data$pam,
      gc_content = as.numeric(.data$gc_content),
      max_a = as.integer(.data$max_a),
      max_t = as.integer(.data$max_t),
      num_t = as.integer(.data$num_t),
      pam_dist_var_start = as.integer(.data$pam_dist_var_start),
      pam_dist_var_end = as.integer(.data$pam_dist_var_end),
      vardiff = as.integer(.data$vardiff),
      chr = .data$chr,
      pos = as.integer(.data$pos),
      ref = .data$ref,
      alt = .data$alt,
      in_sv_hotspot = as_logical_flag(.data$in_sv_hotspot),
      exact_high_score = as_logical_flag(.data$exact_high_score),
      hotspot_high_score = as_logical_flag(.data$hotspot_high_score)
    )

  validate_unique_keys(designs, "uid", "design annotations")

  deepcpf1 <- read_csv_clean(deepcpf1_scores) %>%
    dplyr::transmute(
      guide_id = .data$guide_id,
      deepcpf1_target_seq = .data$target_seq,
      deepcpf1_score = as.numeric(.data$y_pred)
    )

  validate_unique_keys(deepcpf1, "guide_id", "DeepCpf1 scores")

  dr <- read_dr_scores(dr_scores)

  designs %>%
    dplyr::left_join(deepcpf1, by = "guide_id") %>%
    dplyr::left_join(dr, by = "uid")
}

load_variant_support <- function(normalized_variants_round_1, normalized_variants_round_2) {
  variants <- dplyr::bind_rows(
    read_genome_wide_variant_calls(normalized_variants_round_1),
    read_genome_wide_variant_calls(normalized_variants_round_2)
  )

  assert_columns(
    variants,
    c("sample", "chrom", "pos", "ref", "alt", "ao"),
    "normalized variant table"
  )

  variants %>%
    dplyr::transmute(
      sample = .data$sample,
      chr = .data$chrom,
      pos = as.integer(.data$pos),
      ref = .data$ref,
      alt = .data$alt,
      sample_variant_id = paste(.data$sample, .data$chrom, .data$pos, .data$ref, .data$alt, sep = "__"),
      supporting_reads = parse_number_safe(.data$ao)
    )
}

load_variant_details <- function(normalized_variants_round_1, normalized_variants_round_2) {
  variants <- dplyr::bind_rows(
    read_genome_wide_variant_calls(normalized_variants_round_1),
    read_genome_wide_variant_calls(normalized_variants_round_2)
  )

  assert_columns(
    variants,
    c("sample", "chrom", "pos", "ref", "alt", "ao", "ro", "dp", "frc_alt", "frc_ref", "type"),
    "normalized variant table"
  )

  variants %>%
    dplyr::transmute(
      sample = .data$sample,
      chr = .data$chrom,
      pos = as.integer(.data$pos),
      ref = .data$ref,
      alt = .data$alt,
      ao = parse_number_safe(.data$ao),
      ro = parse_number_safe(.data$ro),
      dp = parse_number_safe(.data$dp),
      af = parse_number_safe(.data$frc_alt),
      ref_af = parse_number_safe(.data$frc_ref),
      variant_type = .data$type,
      sample_variant_id = paste(.data$sample, .data$chrom, .data$pos, .data$ref, .data$alt, sep = "__")
    )
}

make_colony_assay_table <- function(input_paths, target_basecov_min = 2, intended_ao_min = 2) {
  picked <- load_picked_colonies(
    input_paths$picked_colonies_round_1,
    input_paths$picked_colonies_round_2
  )
  target_coverage <- load_target_coverage(
    input_paths$target_coverage_round_1,
    input_paths$target_coverage_round_2
  )
  mean_coverage <- load_mean_coverage(
    input_paths$mean_coverage_round_1,
    input_paths$mean_coverage_round_2
  )
  designs <- load_design_features(
    input_paths$design_annotations,
    input_paths$deepcpf1_scores,
    input_paths$dr_scores
  )
  normalized_indels <- read_normalized_indel_designs(input_paths$normalized_indels)
  variant_support <- load_variant_support(
    input_paths$normalized_variants_round_1,
    input_paths$normalized_variants_round_2
  )

  validate_unique_keys(normalized_indels, "uid", "normalized indel design loci")

  missing_designs <- dplyr::anti_join(picked, designs %>% dplyr::select(uid), by = "uid")
  if (nrow(missing_designs) > 0) {
    stop(
      "Picked colonies reference design IDs missing from design annotations: ",
      paste(utils::head(sort(unique(missing_designs$uid)), 10), collapse = ", "),
      call. = FALSE
    )
  }

  colony_annotations <- picked %>%
    dplyr::inner_join(designs, by = "uid") %>%
    dplyr::left_join(target_coverage, by = c("round", "colony_id")) %>%
    dplyr::left_join(mean_coverage, by = "sample")

  if (any(is.na(colony_annotations$deepcpf1_score))) {
    missing_scores <- unique(colony_annotations$guide_id[is.na(colony_annotations$deepcpf1_score)])
    stop(
      "Missing DeepCpf1 scores for guide IDs: ",
      paste(utils::head(sort(missing_scores), 10), collapse = ", "),
      call. = FALSE
    )
  }

  if (any(is.na(colony_annotations$dr_score))) {
    missing_dr_scores <- unique(colony_annotations$uid[is.na(colony_annotations$dr_score)])
    stop(
      "Missing direct-repeat scores for design IDs: ",
      paste(utils::head(sort(missing_dr_scores), 10), collapse = ", "),
      call. = FALSE
    )
  }

  snv_assays <- colony_annotations %>%
    dplyr::filter(.data$vardiff == 0) %>%
    dplyr::mutate(
      target_chr = .data$chr,
      target_pos = as.integer(.data$pos),
      target_ref = .data$ref,
      target_alt = .data$alt
    )

  indel_assays <- colony_annotations %>%
    dplyr::filter(.data$vardiff != 0) %>%
    dplyr::select(-dplyr::all_of(c("chr", "pos", "ref", "alt"))) %>%
    dplyr::left_join(normalized_indels, by = "uid") %>%
    dplyr::rename(
      target_chr = chr,
      target_pos = pos,
      target_ref = ref,
      target_alt = alt
    )

  missing_indels <- indel_assays %>%
    dplyr::filter(is.na(.data$target_chr) | is.na(.data$target_pos))
  if (nrow(missing_indels) > 0) {
    stop(
      "Missing normalized indel target loci for design IDs: ",
      paste(utils::head(sort(unique(missing_indels$uid)), 10), collapse = ", "),
      call. = FALSE
    )
  }

  intended_support <- variant_support %>%
    dplyr::group_by(.data$sample_variant_id) %>%
    dplyr::summarise(
      intended_supporting_reads = max(.data$supporting_reads, na.rm = TRUE),
      .groups = "drop"
    )

  colony_assays <- dplyr::bind_rows(snv_assays, indel_assays) %>%
    dplyr::mutate(
      sample_variant_id = paste(
        .data$sample,
        .data$target_chr,
        .data$target_pos,
        .data$target_ref,
        .data$target_alt,
        sep = "__"
      ),
      assayed_colony = !is.na(.data$target_basecov) & .data$target_basecov >= target_basecov_min,
      coverage_gap = .data$mean_cov - .data$target_basecov,
      pam_dist = pmax(-3, pmin(.data$pam_dist_var_start, .data$pam_dist_var_end)),
      pam_dist_bin = dplyr::case_when(
        .data$pam_dist <= 0 ~ "PAM",
        .data$pam_dist <= 6 ~ "1-6",
        .data$pam_dist <= 11 ~ "7-11",
        .data$pam_dist <= 17 ~ "12-17",
        .data$pam_dist <= 23 ~ "18-23",
        TRUE ~ ">23"
      ),
      variant_class = dplyr::if_else(.data$vardiff == 0, "SNV", "Indel"),
      deepcpf1_guide_match = substr(.data$deepcpf1_target_seq, 4, 30) == .data$guide_seq,
      ntttv_pattern = substr(.data$deepcpf1_target_seq, 3, 7),
      pam_upstream_nt = substr(.data$deepcpf1_target_seq, 3, 3),
      full_pam_5mer = paste0(.data$pam_upstream_nt, .data$pam),
      guide_gc_23 = .data$gc_content,
      guide_gc_20 = gc_fraction(substr(.data$guide_nopam, 1, 20)),
      guide_gc_18 = gc_fraction(substr(.data$guide_nopam, 1, 18)),
      donor_wt_gc = gc_fraction(.data$donor_wt),
      guide_gc_23_bin = gc_fraction_bin(.data$guide_gc_23),
      guide_gc_20_bin = gc_fraction_bin(.data$guide_gc_20),
      guide_gc_18_bin = gc_fraction_bin(.data$guide_gc_18),
      donor_wt_gc_bin = gc_fraction_bin(.data$donor_wt_gc),
      dr_canonical_structure = .data$dr_structure == ".....(((((....)))))........................",
      dr_spacer_pair_any = .data$dr_nt_paired_to_spacer > 0
    ) %>%
    dplyr::left_join(intended_support, by = "sample_variant_id") %>%
    dplyr::mutate(
      intended_supporting_reads = tidyr::replace_na(.data$intended_supporting_reads, 0),
      edited_colony = .data$assayed_colony & .data$intended_supporting_reads >= intended_ao_min
    )

  if (any(!colony_assays$deepcpf1_guide_match, na.rm = TRUE)) {
    mismatched_guides <- unique(colony_assays$uid[!colony_assays$deepcpf1_guide_match])
    stop(
      "DeepCpf1 target contexts do not match the guide sequence for design IDs: ",
      paste(utils::head(sort(mismatched_guides), 10), collapse = ", "),
      call. = FALSE
    )
  }

  deepcpf1_design_annotations <- colony_assays %>%
    dplyr::distinct(.data$uid, .data$deepcpf1_score) %>%
    dplyr::mutate(
      deepcpf1_quartile = dplyr::ntile(.data$deepcpf1_score, 4),
      deepcpf1_quartile = factor(
        .data$deepcpf1_quartile,
        levels = 1:4,
        labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")
      ),
      high_deepcpf1 = .data$deepcpf1_score >= stats::median(.data$deepcpf1_score, na.rm = TRUE)
    )

  colony_assays %>%
    dplyr::left_join(deepcpf1_design_annotations, by = c("uid", "deepcpf1_score")) %>%
    dplyr::mutate(
      preferred_pam = .data$pam %in% c("TTTA", "TTTC"),
      preferred_distance = .data$pam_dist_bin %in% c("1-6", "7-11", "12-17"),
      preferred_rule_count = .data$preferred_pam + .data$preferred_distance + .data$high_deepcpf1,
      recommended_design = .data$preferred_rule_count == 3
    )
}

make_non_target_tables <- function(colony_assays, input_paths, non_target_af_min = 50, non_target_dp_min = 4) {
  variant_details <- load_variant_details(
    input_paths$normalized_variants_round_1,
    input_paths$normalized_variants_round_2
  )

  intended_sites <- colony_assays %>%
    dplyr::filter(.data$assayed_colony, !is.na(.data$sample)) %>%
    dplyr::distinct(
      sample,
      intended_uid = .data$uid,
      intended_sample_variant_id = .data$sample_variant_id
    )

  variants_annotated <- variant_details %>%
    dplyr::inner_join(intended_sites, by = "sample") %>%
    dplyr::mutate(
      is_intended_variant = .data$sample_variant_id == .data$intended_sample_variant_id,
      non_target_af50_dp4 = !.data$is_intended_variant &
        .data$af >= non_target_af_min &
        .data$dp >= non_target_dp_min,
      variant_id = paste(.data$chr, .data$pos, .data$ref, .data$alt, sep = ":")
    )

  non_target_variants <- variants_annotated %>%
    dplyr::filter(.data$non_target_af50_dp4) %>%
    dplyr::select(
      sample,
      intended_uid,
      chr,
      pos,
      ref,
      alt,
      variant_type,
      ao,
      ro,
      dp,
      af,
      ref_af,
      variant_id
    )

  colony_summary <- non_target_variants %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::summarise(
      non_target_af50_dp4_count = dplyr::n(),
      max_non_target_af = max(.data$af, na.rm = TRUE),
      max_non_target_dp = max(.data$dp, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    non_target_variants = non_target_variants,
    colony_summary = colony_summary
  )
}

build_genome_wide_colony_status <- function(
  input_paths,
  target_basecov_min = 2,
  intended_ao_min = 2,
  non_target_af_min = 50,
  non_target_dp_min = 4
) {
  check_required_files(input_paths)

  colony_assays <- make_colony_assay_table(
    input_paths,
    target_basecov_min = target_basecov_min,
    intended_ao_min = intended_ao_min
  )

  non_target <- make_non_target_tables(
    colony_assays,
    input_paths,
    non_target_af_min = non_target_af_min,
    non_target_dp_min = non_target_dp_min
  )

  colony_status <- colony_assays %>%
    dplyr::left_join(non_target$colony_summary, by = "sample") %>%
    dplyr::mutate(
      non_target_af50_dp4_count = tidyr::replace_na(.data$non_target_af50_dp4_count, 0L),
      has_non_target_af50_dp4 = .data$non_target_af50_dp4_count > 0,
      colony_edit_status = dplyr::case_when(
        !.data$assayed_colony ~ "Not assayed",
        .data$edited_colony & .data$has_non_target_af50_dp4 ~ "Intended edit plus non-target variant",
        .data$edited_colony ~ "Intended edit only",
        .data$has_non_target_af50_dp4 ~ "Non-target variant only",
        TRUE ~ "Unedited"
      )
    )

  design_status <- colony_status %>%
    dplyr::group_by(.data$uid) %>%
    dplyr::summarise(
      round = dplyr::first(.data$round),
      guide_id = dplyr::first(.data$guide_id),
      pam = dplyr::first(.data$pam),
      pam_dist = dplyr::first(.data$pam_dist),
      pam_dist_bin = dplyr::first(.data$pam_dist_bin),
      deepcpf1_score = dplyr::first(.data$deepcpf1_score),
      deepcpf1_quartile = dplyr::first(as.character(.data$deepcpf1_quartile)),
      pam_upstream_nt = dplyr::first(.data$pam_upstream_nt),
      full_pam_5mer = dplyr::first(.data$full_pam_5mer),
      ntttv_pattern = dplyr::first(.data$ntttv_pattern),
      guide_gc_23 = dplyr::first(.data$guide_gc_23),
      guide_gc_20 = dplyr::first(.data$guide_gc_20),
      guide_gc_18 = dplyr::first(.data$guide_gc_18),
      donor_wt_gc = dplyr::first(.data$donor_wt_gc),
      guide_gc_23_bin = dplyr::first(.data$guide_gc_23_bin),
      guide_gc_20_bin = dplyr::first(.data$guide_gc_20_bin),
      guide_gc_18_bin = dplyr::first(.data$guide_gc_18_bin),
      donor_wt_gc_bin = dplyr::first(.data$donor_wt_gc_bin),
      max_a = dplyr::first(.data$max_a),
      max_t = dplyr::first(.data$max_t),
      num_t = dplyr::first(.data$num_t),
      variant_class = dplyr::first(.data$variant_class),
      target_chr = dplyr::first(.data$target_chr),
      target_pos = dplyr::first(.data$target_pos),
      target_ref = dplyr::first(.data$target_ref),
      target_alt = dplyr::first(.data$target_alt),
      in_sv_hotspot = dplyr::first(.data$in_sv_hotspot),
      exact_high_score = dplyr::first(.data$exact_high_score),
      hotspot_high_score = dplyr::first(.data$hotspot_high_score),
      dr_mfe_kcal_mol = dplyr::first(.data$dr_mfe_kcal_mol),
      dr_nt_paired_to_dr = dplyr::first(.data$dr_nt_paired_to_dr),
      dr_nt_paired_to_spacer = dplyr::first(.data$dr_nt_paired_to_spacer),
      dr_nt_unpaired = dplyr::first(.data$dr_nt_unpaired),
      dr_perturbed = dplyr::first(.data$dr_perturbed),
      dr_score = dplyr::first(.data$dr_score),
      dr_canonical_structure = dplyr::first(.data$dr_canonical_structure),
      dr_spacer_pair_any = dplyr::first(.data$dr_spacer_pair_any),
      preferred_pam = dplyr::first(.data$preferred_pam),
      preferred_distance = dplyr::first(.data$preferred_distance),
      high_deepcpf1 = dplyr::first(.data$high_deepcpf1),
      preferred_rule_count = dplyr::first(.data$preferred_rule_count),
      recommended_design = dplyr::first(.data$recommended_design),
      picked_colonies = dplyr::n(),
      assayed_colonies = sum(.data$assayed_colony, na.rm = TRUE),
      edited_colonies = sum(.data$edited_colony, na.rm = TRUE),
      max_intended_supporting_reads = max(.data$intended_supporting_reads, na.rm = TRUE),
      assayed_colonies_with_non_target_af50_dp4 = sum(.data$has_non_target_af50_dp4, na.rm = TRUE),
      non_target_af50_dp4_variants = sum(.data$non_target_af50_dp4_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      assayed_design = .data$assayed_colonies > 0,
      edited_design = .data$edited_colonies > 0,
      other_variant_design = .data$non_target_af50_dp4_variants > 0,
      colony_edit_rate = dplyr::if_else(
        .data$assayed_colonies > 0,
        .data$edited_colonies / .data$assayed_colonies,
        NA_real_
      ),
      design_edit_status = dplyr::case_when(
        .data$other_variant_design ~ "Other variants observed",
        .data$edited_design ~ "Intended edit only",
        TRUE ~ "Unedited"
      )
    )

  colony_status <- colony_status %>%
    dplyr::select(
      uid,
      round,
      sample,
      colony_id,
      guide_id,
      target_chr,
      target_pos,
      target_ref,
      target_alt,
      pam,
      pam_upstream_nt,
      full_pam_5mer,
      pam_dist,
      pam_dist_bin,
      deepcpf1_score,
      deepcpf1_quartile,
      guide_gc_23,
      guide_gc_20,
      guide_gc_18,
      donor_wt_gc,
      max_a,
      max_t,
      num_t,
      variant_class,
      in_sv_hotspot,
      exact_high_score,
      hotspot_high_score,
      dr_mfe_kcal_mol,
      dr_nt_paired_to_dr,
      dr_nt_paired_to_spacer,
      dr_nt_unpaired,
      dr_perturbed,
      dr_score,
      dr_canonical_structure,
      dr_spacer_pair_any,
      target_basecov,
      mean_cov,
      assayed_colony,
      intended_supporting_reads,
      edited_colony,
      non_target_af50_dp4_count,
      max_non_target_af,
      max_non_target_dp,
      has_non_target_af50_dp4,
      colony_edit_status
    )

  list(
    colony_status = colony_status,
    design_status = design_status,
    non_target_variants = non_target$non_target_variants
  )
}

run_genome_wide_colony_status <- function(input_paths, output_dir) {
  tables <- build_genome_wide_colony_status(input_paths)

  output_paths <- list(
    colony_status = file.path(output_dir, "genome_wide_colony_editing_status.csv"),
    design_status = file.path(output_dir, "genome_wide_design_editing_status.csv"),
    non_target_variants = file.path(output_dir, "genome_wide_non_target_af50_dp4_variants.csv")
  )

  write_csv_mkdir(tables$colony_status, output_paths$colony_status)
  write_csv_mkdir(tables$design_status, output_paths$design_status)
  write_csv_mkdir(tables$non_target_variants, output_paths$non_target_variants)

  c(tables, list(paths = output_paths))
}
