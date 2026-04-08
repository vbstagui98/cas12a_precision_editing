validate_unique_keys <- function(data, key_columns, data_name) {
  duplicated_keys <- data %>%
    dplyr::count(dplyr::across(dplyr::all_of(key_columns)), name = "n") %>%
    dplyr::filter(.data$n > 1)

  if (nrow(duplicated_keys) > 0) {
    stop(
      sprintf(
        "%s contains duplicated keys for: %s",
        data_name,
        paste(key_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(data)
}

genome_wide_dr_sequence <- function() {
  "TAATTTCTACTCTTGTAGAT"
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

genome_wide_design_rule_files <- function(project_root, paths) {
  list(
    design_annotations = require_file(
      file.path(project_root, "scripts", "designs_812_sv_my_variants_annotated.tsv"),
      "workspace design annotations"
    ),
    deepcpf1_scores = require_file(
      file.path(project_root, "input_deepcpf1_20260319.scored.csv"),
      "workspace DeepCpf1 scores"
    ),
    variant_annotations = require_file(
      file.path(project_root, "results", "variants", "vars_20241106.annotated.tsv"),
      "workspace variant annotations"
    ),
    dr_scores = require_file(
      file.path(
        paths$genome_wide_root,
        "Tn5",
        "Analysis",
        "annotations_guides",
        "guides_dr_scores_20260325.tsv"
      ),
      "guide direct-repeat folding scores"
    ),
    score_sv_sites = require_file(
      file.path(paths$downloads_root, "s288c_NGG_PAM_sites_editing_SCORE.tsv"),
      "SCORE SV site table"
    ),
    picked_colonies_round_1 = require_file(
      file.path(paths$genome_wide_root, "REDI", "picked_colonies_20250613_annotated.csv"),
      "genome-wide picked colonies round 1"
    ),
    picked_colonies_round_2 = require_file(
      file.path(
        paths$genome_wide_root,
        "REDI",
        "second_REDI_20250917",
        "scritps_20251216",
        "results",
        "colonies_with_BC1_20260212.csv"
      ),
      "genome-wide picked colonies round 2"
    ),
    target_coverage_round_1 = require_file(
      file.path(paths$genome_wide_root, "Tn5", "Analysis", "stats", "target_basecov_20260222.tsv"),
      "target coverage round 1"
    ),
    target_coverage_round_2 = require_file(
      file.path(
        paths$genome_wide_root,
        "Tn5",
        "second_round_20260130",
        "read_stats",
        "target_coverage_20260213.tsv"
      ),
      "target coverage round 2"
    ),
    mean_coverage_round_1 = require_file(
      file.path(
        paths$genome_wide_root,
        "Tn5",
        "second_round_20260130",
        "read_stats",
        "mean_cov_20260223.csv"
      ),
      "mean coverage round 1"
    ),
    mean_coverage_round_2 = require_file(
      file.path(
        paths$genome_wide_root,
        "Tn5",
        "second_round_20260130",
        "read_stats",
        "mean_cov_20260223_second.csv"
      ),
      "mean coverage round 2"
    ),
    normalized_variants_round_1 = require_file(
      file.path(paths$genome_wide_root, "Tn5", "Analysis", "results", "variants_norm_WGS_20260222.csv"),
      "normalized variants round 1"
    ),
    normalized_variants_round_2 = require_file(
      file.path(
        paths$genome_wide_root,
        "Tn5",
        "second_round_20260130",
        "results",
        "variants_norm_WGS_20260219.csv"
      ),
      "normalized variants round 2"
    ),
    normalized_indels = require_file(
      file.path(
        paths$genome_wide_root,
        "Tn5",
        "second_round_20260130",
        "annotations_ref",
        "normalized_indel_designs_20260223.vcf"
      ),
      "normalized indel design loci"
    )
  )
}

read_genome_wide_variant_calls <- function(path) {
  read_csv_clean(path) %>%
    dplyr::select(-dplyr::any_of(c("x", "...1")))
}

read_genome_wide_normalized_indels <- function(path) {
  data <- readr::read_tsv(
    require_file(path),
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE
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

read_genome_wide_dr_scores <- function(path) {
  dr_scores <- read_delim_clean(path, delim = "\t") %>%
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
      dr_perturbed = tolower(as.character(.data$dr_perturbed)) %in% c("yes", "true", "1"),
      dr_score = as.numeric(.data$score)
    )

  validate_unique_keys(dr_scores, "uid", "guide direct-repeat folding scores")
  dr_scores
}

read_score_sv_sites <- function(path) {
  score_sv <- read_delim_clean(path, delim = "\t") %>%
    dplyr::transmute(
      chr = paste0("chr", .data$chr),
      pos = as.integer(.data$pos),
      score_sv = as.numeric(.data$score_sv_minmax),
      score_guide = .data$guide,
      score_site_id = .data$id
    )

  validate_unique_keys(score_sv, c("chr", "pos"), "SCORE SV sites")
  score_sv
}

annotate_nearest_score_sv <- function(data, chr_col, pos_col, score_sv_sites) {
  chr_col <- rlang::ensym(chr_col)
  pos_col <- rlang::ensym(pos_col)
  chr_name <- rlang::as_name(chr_col)
  pos_name <- rlang::as_name(pos_col)

  split_scores <- split(score_sv_sites, score_sv_sites$chr)

  nearest_match <- function(chr_value, pos_value) {
    chr_scores <- split_scores[[chr_value]]

    if (is.null(chr_scores) || is.na(pos_value) || nrow(chr_scores) == 0) {
      return(
        data.frame(
          nearest_score_sv = NA_real_,
          nearest_score_sv_pos = NA_integer_,
          nearest_score_sv_guide = NA_character_,
          nearest_score_sv_site_id = NA_character_,
          stringsAsFactors = FALSE
        )
      )
    }

    nearest_idx <- which.min(abs(chr_scores$pos - pos_value))
    data.frame(
      nearest_score_sv = chr_scores$score_sv[[nearest_idx]],
      nearest_score_sv_pos = chr_scores$pos[[nearest_idx]],
      nearest_score_sv_guide = chr_scores$score_guide[[nearest_idx]],
      nearest_score_sv_site_id = chr_scores$score_site_id[[nearest_idx]],
      stringsAsFactors = FALSE
    )
  }

  nearest_matrix <- purrr::map2(
    data[[chr_name]],
    data[[pos_name]],
    nearest_match
  ) %>%
    purrr::list_rbind() %>%
    tibble::as_tibble(.name_repair = "minimal")

  names(nearest_matrix) <- c(
    "nearest_score_sv",
    "nearest_score_sv_pos",
    "nearest_score_sv_guide",
    "nearest_score_sv_site_id"
  )

  nearest_matrix <- nearest_matrix %>%
    dplyr::mutate(
      nearest_score_sv = as.numeric(.data$nearest_score_sv),
      nearest_score_sv_pos = as.integer(.data$nearest_score_sv_pos)
    )

  dplyr::bind_cols(data, nearest_matrix)
}

load_genome_wide_picked_colonies <- function(file_map) {
  round_1 <- read_csv_clean(file_map$picked_colonies_round_1) %>%
    dplyr::transmute(
      uid = .data$design,
      bc1 = .data$bc1,
      colony_id = .data$colony_id,
      round = "round_1"
    )

  round_2 <- read_csv_clean(file_map$picked_colonies_round_2) %>%
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

load_genome_wide_target_coverage <- function(file_map) {
  round_1 <- read_delim_clean(file_map$target_coverage_round_1, delim = "\t") %>%
    dplyr::mutate(
      colony_id = stringr::str_replace(.data$sample, ".*VB_(.*?)_S.*", "\\1"),
      round = "round_1"
    )

  round_2 <- read_delim_clean(file_map$target_coverage_round_2, delim = "\t") %>%
    dplyr::mutate(
      colony_id = stringr::str_replace(.data$sample, ".*VB_(.*?)_S.*", "\\1"),
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

load_genome_wide_mean_coverage <- function(file_map) {
  mean_cov <- dplyr::bind_rows(
    read_csv_clean(file_map$mean_coverage_round_1),
    read_csv_clean(file_map$mean_coverage_round_2)
  ) %>%
    dplyr::mutate(
      sample = stringr::str_remove(.data$file, "\\.bincov\\.txt$")
    ) %>%
    dplyr::transmute(
      sample = .data$sample,
      mean_cov = as.numeric(.data$mean_cov)
    )

  validate_unique_keys(mean_cov, "sample", "mean coverage")
  mean_cov
}

load_workspace_design_annotations <- function(file_map) {
  design_annotations <- read_delim_clean(file_map$design_annotations, delim = "\t")
  assert_columns(
    design_annotations,
    c(
      "uid", "guide_id", "guide_seq", "guide_nopam", "donor_wt",
      "pam", "gc_content", "max_a", "max_t", "num_t",
      "pam_dist_var_start", "pam_dist_var_end", "vardiff",
      "chr", "pos", "ref", "alt", "in_sv_hotspot",
      "exact_high_score", "hotspot_high_score"
    ),
    "design_annotations"
  )

  validate_unique_keys(design_annotations, "uid", "design annotations")

  deepcpf1 <- read_csv_clean(file_map$deepcpf1_scores) %>%
    dplyr::transmute(
      guide_id = .data$guide_id,
      deepcpf1_target_seq = .data$target_seq,
      deepcpf1_score = as.numeric(.data$y_pred)
    )

  validate_unique_keys(deepcpf1, "guide_id", "DeepCpf1 scores")

  variant_annotations <- read_delim_clean(file_map$variant_annotations, delim = "\t") %>%
    dplyr::transmute(
      chr = .data$chr,
      pos = as.integer(.data$pos),
      ref = .data$ref,
      alt = .data$alt,
      region_class = .data$region_class,
      snpeff_primary_impact = .data$snpeff_primary_impact,
      snpeff_primary_effect = .data$snpeff_primary_effect,
      spcas9_targetable = .data$spcas9_targetable
    )

  validate_unique_keys(
    variant_annotations,
    c("chr", "pos", "ref", "alt"),
    "variant annotations"
  )

  list(
    design_annotations = design_annotations,
    deepcpf1 = deepcpf1,
    variant_annotations = variant_annotations
  )
}

load_genome_wide_normalized_variants <- function(file_map) {
  variants <- dplyr::bind_rows(
    read_genome_wide_variant_calls(file_map$normalized_variants_round_1),
    read_genome_wide_variant_calls(file_map$normalized_variants_round_2)
  ) %>%
    dplyr::transmute(
      sample = .data$sample,
      variant_chr = .data$chrom,
      variant_pos = as.integer(.data$pos),
      variant_ref = .data$ref,
      variant_alt = .data$alt,
      sample_variant_id = paste(
        .data$sample,
        .data$chrom,
        .data$pos,
        .data$ref,
        .data$alt,
        sep = "__"
      )
    ) %>%
    dplyr::distinct()

  variants
}

load_genome_wide_variant_call_details <- function(file_map) {
  dplyr::bind_rows(
    read_genome_wide_variant_calls(file_map$normalized_variants_round_1),
    read_genome_wide_variant_calls(file_map$normalized_variants_round_2)
  ) %>%
    dplyr::transmute(
      sample = .data$sample,
      chr = .data$chrom,
      pos = as.integer(.data$pos),
      ref = .data$ref,
      alt = .data$alt,
      qual = as.numeric(.data$qual),
      dp = as.numeric(.data$dp),
      af = as.numeric(.data$frc_alt),
      variant_type = .data$type,
      sample_variant_id = paste(
        .data$sample,
        .data$chrom,
        .data$pos,
        .data$ref,
        .data$alt,
        sep = "__"
      )
    )
}

build_genome_wide_design_rule_tables <- function(project_root, paths, target_basecov_min = 2) {
  file_map <- genome_wide_design_rule_files(project_root, paths)
  annotations <- load_workspace_design_annotations(file_map)
  dr_scores <- read_genome_wide_dr_scores(file_map$dr_scores)
  score_sv_sites <- read_score_sv_sites(file_map$score_sv_sites)
  picked <- load_genome_wide_picked_colonies(file_map)
  target_coverage <- load_genome_wide_target_coverage(file_map)
  mean_coverage <- load_genome_wide_mean_coverage(file_map)
  normalized_indels <- read_genome_wide_normalized_indels(file_map$normalized_indels)
  normalized_variants <- load_genome_wide_normalized_variants(file_map)

  missing_designs <- dplyr::anti_join(
    picked,
    annotations$design_annotations %>% dplyr::select(uid),
    by = "uid"
  )

  if (nrow(missing_designs) > 0) {
    stop(
      sprintf(
        "Picked colonies reference %d design IDs that are missing from the workspace annotations.",
        nrow(missing_designs)
      ),
      call. = FALSE
    )
  }

  colony_annotations <- picked %>%
    dplyr::inner_join(annotations$design_annotations, by = "uid") %>%
    dplyr::left_join(annotations$deepcpf1, by = "guide_id") %>%
    dplyr::left_join(dr_scores, by = "uid") %>%
    dplyr::left_join(target_coverage, by = c("round", "colony_id")) %>%
    dplyr::left_join(mean_coverage, by = "sample")

  if (any(is.na(colony_annotations$deepcpf1_score))) {
    missing_scores <- unique(colony_annotations$guide_id[is.na(colony_annotations$deepcpf1_score)])
    stop(
      sprintf(
        "Missing DeepCpf1 scores for guide IDs: %s",
        paste(utils::head(sort(missing_scores), 10), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (any(is.na(colony_annotations$dr_score))) {
    missing_dr_scores <- unique(colony_annotations$uid[is.na(colony_annotations$dr_score)])
    stop(
      sprintf(
        "Missing DR-folding scores for design IDs: %s",
        paste(utils::head(sort(missing_dr_scores), 10), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  validate_unique_keys(normalized_indels, "uid", "normalized indel loci")

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

  colony_assays <- dplyr::bind_rows(snv_assays, indel_assays) %>%
    dplyr::left_join(
      annotations$variant_annotations,
      by = c(
        "target_chr" = "chr",
        "target_pos" = "pos",
        "target_ref" = "ref",
        "target_alt" = "alt"
      )
    ) %>%
    annotate_nearest_score_sv(target_chr, target_pos, score_sv_sites = score_sv_sites) %>%
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
      edited_colony = .data$assayed_colony & .data$sample_variant_id %in% normalized_variants$sample_variant_id,
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
      pam_dist_bin = factor(.data$pam_dist_bin, levels = c("1-6", "7-11", "12-17", "18-23", "PAM", ">23")),
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
    )

  if (any(!colony_assays$deepcpf1_guide_match, na.rm = TRUE)) {
    mismatched_guides <- unique(colony_assays$uid[!colony_assays$deepcpf1_guide_match])
    stop(
      sprintf(
        "DeepCpf1 target contexts do not match the guide sequence for design IDs: %s",
        paste(utils::head(sort(mismatched_guides), 10), collapse = ", ")
      ),
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

  colony_assays <- colony_assays %>%
    dplyr::left_join(
      deepcpf1_design_annotations,
      by = c("uid", "deepcpf1_score")
    ) %>%
    dplyr::mutate(
      preferred_pam = .data$pam %in% c("TTTA", "TTTC"),
      preferred_distance = .data$pam_dist_bin %in% c("1-6", "7-11", "12-17"),
      preferred_rule_count = .data$preferred_pam + .data$preferred_distance + .data$high_deepcpf1,
      recommended_design = .data$preferred_rule_count == 3
    )

  design_outcomes <- colony_assays %>%
    dplyr::group_by(.data$uid) %>%
    dplyr::summarise(
      round = dplyr::first(.data$round),
      guide_id = dplyr::first(.data$guide_id),
      pam = dplyr::first(.data$pam),
      pam_dist = dplyr::first(.data$pam_dist),
      pam_dist_bin = dplyr::first(.data$pam_dist_bin),
      deepcpf1_score = dplyr::first(.data$deepcpf1_score),
      deepcpf1_quartile = dplyr::first(.data$deepcpf1_quartile),
      pam_upstream_nt = dplyr::first(.data$pam_upstream_nt),
      full_pam_5mer = dplyr::first(.data$full_pam_5mer),
      ntttv_pattern = dplyr::first(.data$ntttv_pattern),
      gc_content = dplyr::first(.data$gc_content),
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
      vardiff = dplyr::first(.data$vardiff),
      variant_class = dplyr::first(.data$variant_class),
      region_class = dplyr::first(.data$region_class),
      snpeff_primary_impact = dplyr::first(.data$snpeff_primary_impact),
      snpeff_primary_effect = dplyr::first(.data$snpeff_primary_effect),
      spcas9_targetable = dplyr::first(.data$spcas9_targetable),
      in_sv_hotspot = dplyr::first(.data$in_sv_hotspot),
      exact_high_score = dplyr::first(.data$exact_high_score),
      hotspot_high_score = dplyr::first(.data$hotspot_high_score),
      nearest_score_sv = dplyr::first(.data$nearest_score_sv),
      nearest_score_sv_pos = dplyr::first(.data$nearest_score_sv_pos),
      nearest_score_sv_guide = dplyr::first(.data$nearest_score_sv_guide),
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
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      assayed_design = .data$assayed_colonies > 0,
      edited_design = .data$edited_colonies > 0,
      colony_edit_rate = dplyr::if_else(
        .data$assayed_colonies > 0,
        .data$edited_colonies / .data$assayed_colonies,
        NA_real_
      )
    )

  list(
    files = file_map,
    colony_assays = colony_assays,
    design_outcomes = design_outcomes
  )
}

binom_wilson_interval <- function(successes, trials, conf_level = 0.95) {
  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)
  trials_safe <- ifelse(trials > 0, trials, NA_real_)
  rate <- successes / trials_safe
  denominator <- 1 + (z_value^2 / trials_safe)
  center <- (rate + (z_value^2 / (2 * trials_safe))) / denominator
  half_width <- (
    z_value *
      sqrt((rate * (1 - rate) + (z_value^2 / (4 * trials_safe))) / trials_safe)
  ) / denominator

  tibble::tibble(
    rate = rate,
    lower = pmax(0, center - half_width),
    upper = pmin(1, center + half_width)
  )
}

build_genome_wide_stage_summary <- function(design_outcomes) {
  summarise_stage <- function(data, label) {
    tibble::tibble(
      round = label,
      picked_designs = nrow(data),
      assayed_designs = sum(data$assayed_design),
      edited_designs = sum(data$edited_design),
      picked_colonies = sum(data$picked_colonies),
      assayed_colonies = sum(data$assayed_colonies),
      edited_colonies = sum(data$edited_colonies)
    )
  }

  dplyr::bind_rows(
    summarise_stage(design_outcomes, "overall"),
    design_outcomes %>%
      dplyr::group_by(.data$round) %>%
      dplyr::summarise(
        picked_designs = dplyr::n(),
        assayed_designs = sum(.data$assayed_design),
        edited_designs = sum(.data$edited_design),
        picked_colonies = sum(.data$picked_colonies),
        assayed_colonies = sum(.data$assayed_colonies),
        edited_colonies = sum(.data$edited_colonies),
        .groups = "drop"
      )
  )
}

summarise_edit_rates <- function(design_outcomes, group_vars) {
  summary_data <- design_outcomes %>%
    dplyr::filter(.data$assayed_design) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      total_guides = dplyr::n(),
      successful_guides = sum(.data$edited_design),
      total_assayed_colonies = sum(.data$assayed_colonies),
      edited_colonies = sum(.data$edited_colonies),
      .groups = "drop"
    )

  design_ci <- binom_wilson_interval(summary_data$successful_guides, summary_data$total_guides)
  colony_ci <- binom_wilson_interval(summary_data$edited_colonies, summary_data$total_assayed_colonies)

  dplyr::bind_cols(
    summary_data,
    design_ci %>%
      dplyr::rename(
        design_success_rate = rate,
        design_success_lower = lower,
        design_success_upper = upper
      ),
    colony_ci %>%
      dplyr::rename(
        colony_edit_rate = rate,
        colony_edit_lower = lower,
        colony_edit_upper = upper
      )
  )
}

build_multi_colony_consistency_summary <- function(design_outcomes, min_assayed_colonies = 2) {
  design_outcomes %>%
    dplyr::filter(.data$assayed_design, .data$assayed_colonies >= min_assayed_colonies) %>%
    dplyr::mutate(
      consistency = dplyr::case_when(
        .data$edited_colonies == 0 ~ "None edited",
        .data$edited_colonies == .data$assayed_colonies ~ "All edited",
        TRUE ~ "Mixed"
      ),
      consistency = factor(.data$consistency, levels = c("All edited", "Mixed", "None edited"))
    ) %>%
    dplyr::count(.data$assayed_colonies, .data$consistency, name = "design_count") %>%
    dplyr::group_by(.data$assayed_colonies) %>%
    dplyr::mutate(
      total_designs = sum(.data$design_count),
      fraction = .data$design_count / .data$total_designs
    ) %>%
    dplyr::ungroup()
}

build_genome_wide_background_variant_tables <- function(
  file_map,
  colony_assays,
  score_sv_sites = read_score_sv_sites(file_map$score_sv_sites),
  qual_min = 20,
  af_min = 50,
  dp_min = 4
) {
  variant_calls <- load_genome_wide_variant_call_details(file_map)

  intended_sites <- colony_assays %>%
    dplyr::transmute(
      sample = .data$sample,
      intended_sample_variant_id = .data$sample_variant_id,
      assayed_colony = .data$assayed_colony,
      edited_colony = .data$edited_colony
    ) %>%
    dplyr::distinct()

  background_variants <- variant_calls %>%
    dplyr::left_join(
      intended_sites %>%
        dplyr::select(sample, intended_sample_variant_id),
      by = "sample"
    ) %>%
    dplyr::mutate(
      is_intended_variant = .data$sample_variant_id == .data$intended_sample_variant_id
    ) %>%
    dplyr::select(-.data$intended_sample_variant_id) %>%
    annotate_nearest_score_sv(chr, pos, score_sv_sites = score_sv_sites) %>%
    dplyr::mutate(
      high_conf_non_target = !.data$is_intended_variant &
        .data$qual >= qual_min &
        .data$af >= af_min &
        .data$dp >= dp_min,
      high_score_sv_site = .data$nearest_score_sv >= 0.9,
      variant_id = paste(.data$chr, .data$pos, .data$ref, .data$alt, sep = ":")
    )

  high_conf_non_target_variants <- background_variants %>%
    dplyr::filter(.data$high_conf_non_target)

  colony_background_summary <- intended_sites %>%
    dplyr::left_join(
      high_conf_non_target_variants %>%
        dplyr::group_by(.data$sample) %>%
        dplyr::summarise(
          high_conf_non_target_count = dplyr::n(),
          max_non_target_score_sv = max(.data$nearest_score_sv, na.rm = TRUE),
          mean_non_target_score_sv = mean(.data$nearest_score_sv, na.rm = TRUE),
          frac_non_target_high_score = mean(.data$high_score_sv_site, na.rm = TRUE),
          max_non_target_qual = max(.data$qual, na.rm = TRUE),
          max_non_target_af = max(.data$af, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "sample"
    ) %>%
    dplyr::mutate(
      high_conf_non_target_count = dplyr::coalesce(.data$high_conf_non_target_count, 0L)
    )

  recurrent_high_conf_variants <- high_conf_non_target_variants %>%
    dplyr::count(
      .data$variant_id,
      .data$chr,
      .data$pos,
      .data$ref,
      .data$alt,
      .data$nearest_score_sv,
      name = "sample_count",
      sort = TRUE
    )

  list(
    all_variants = background_variants,
    high_conf_non_target_variants = high_conf_non_target_variants,
    colony_background_summary = colony_background_summary,
    recurrent_high_conf_variants = recurrent_high_conf_variants
  )
}

summarise_rule_comparison <- function(
  preferred_data,
  other_data,
  priority,
  rule,
  preferred_group,
  other_group,
  note = NULL
) {
  safe_fraction <- function(numerator, denominator) {
    if (denominator == 0) {
      return(NA_real_)
    }
    numerator / denominator
  }

  preferred_design_success <- mean(preferred_data$edited_design)
  other_design_success <- mean(other_data$edited_design)
  preferred_colony_rate <- safe_fraction(
    sum(preferred_data$edited_colonies),
    sum(preferred_data$assayed_colonies)
  )
  other_colony_rate <- safe_fraction(
    sum(other_data$edited_colonies),
    sum(other_data$assayed_colonies)
  )

  tibble::tibble(
    priority = priority,
    rule = rule,
    preferred_group = preferred_group,
    other_group = other_group,
    preferred_designs = nrow(preferred_data),
    other_designs = nrow(other_data),
    preferred_design_success_rate = preferred_design_success,
    other_design_success_rate = other_design_success,
    preferred_colony_edit_rate = preferred_colony_rate,
    other_colony_edit_rate = other_colony_rate,
    note = note %||% ""
  )
}

build_genome_wide_design_rule_summary <- function(design_outcomes) {
  assayed <- design_outcomes %>% dplyr::filter(.data$assayed_design)

  preferred_pam <- assayed %>% dplyr::filter(.data$preferred_pam)
  other_pam <- assayed %>% dplyr::filter(!.data$preferred_pam)

  preferred_distance <- assayed %>% dplyr::filter(.data$preferred_distance)
  other_distance <- assayed %>% dplyr::filter(!.data$preferred_distance)

  deep_top <- assayed %>% dplyr::filter(.data$deepcpf1_quartile == "Q4 highest")
  deep_bottom <- assayed %>% dplyr::filter(.data$deepcpf1_quartile == "Q1 lowest")

  preferred_bundle <- assayed %>% dplyr::filter(.data$recommended_design)
  other_bundle <- assayed %>% dplyr::filter(!.data$recommended_design)

  high_poly_a <- assayed %>% dplyr::filter(.data$max_a >= 5)
  lower_poly_a <- assayed %>% dplyr::filter(.data$max_a < 5)

  hotspot_guides <- assayed %>% dplyr::filter(.data$in_sv_hotspot)
  non_hotspot_guides <- assayed %>% dplyr::filter(!.data$in_sv_hotspot)

  dplyr::bind_rows(
    summarise_rule_comparison(
      preferred_pam,
      other_pam,
      priority = "Primary",
      rule = "Prefer TTTA or TTTC PAMs over TTTG",
      preferred_group = "TTTA or TTTC",
      other_group = "TTTG",
      note = "Consistent in both marginal summaries and the weighted model."
    ),
    summarise_rule_comparison(
      preferred_distance,
      other_distance,
      priority = "Primary",
      rule = "Place the edit 1-17 nt from the PAM",
      preferred_group = "1-17 nt from PAM",
      other_group = "PAM-overlapping or 18-23 nt",
      note = "PAM-overlapping and 18-23 nt designs were the main distance failures."
    ),
    summarise_rule_comparison(
      deep_top,
      deep_bottom,
      priority = "Primary",
      rule = "Prioritize guides with higher DeepCpf1 scores",
      preferred_group = "DeepCpf1 Q4",
      other_group = "DeepCpf1 Q1",
      note = "Treat the model score as a prioritization signal, not an absolute cutoff."
    ),
    summarise_rule_comparison(
      preferred_bundle,
      other_bundle,
      priority = "Primary",
      rule = "Best simple bundle: preferred PAM + 1-17 nt distance + DeepCpf1 above the median",
      preferred_group = "All 3 preferred rules",
      other_group = "Anything else",
      note = "This is the cleanest simple design bundle in the dataset."
    ),
    summarise_rule_comparison(
      lower_poly_a,
      high_poly_a,
      priority = "Secondary",
      rule = "Use long A-runs only as a tie-breaker filter",
      preferred_group = "max_A < 5",
      other_group = "max_A >= 5",
      note = "The drop for long A-runs was visible only in a small subset."
    ),
    summarise_rule_comparison(
      non_hotspot_guides,
      hotspot_guides,
      priority = "Context",
      rule = "Avoid SV hotspot loci if you have equivalent alternatives",
      preferred_group = "Not in SV hotspot",
      other_group = "SV hotspot",
      note = "Small sample size; useful as a contextual filter, not a primary rule."
    )
  )
}

fit_genome_wide_weighted_model <- function(design_outcomes) {
  model_data <- design_outcomes %>%
    dplyr::filter(.data$assayed_design) %>%
    dplyr::mutate(
      round = factor(.data$round, levels = c("round_1", "round_2")),
      pam = factor(.data$pam, levels = c("TTTA", "TTTC", "TTTG")),
      pam_dist_bin = factor(
        .data$pam_dist_bin,
        levels = c("1-6", "7-11", "12-17", "18-23", "PAM", ">23")
      )
    )

  stats::glm(
    cbind(edited_colonies, assayed_colonies - edited_colonies) ~
      round + pam + pam_dist_bin + deepcpf1_score + max_a + max_t + vardiff,
    family = stats::binomial(),
    data = model_data
  )
}

tidy_genome_wide_weighted_model <- function(model) {
  broom::tidy(model, conf.int = TRUE) %>%
    dplyr::mutate(
      odds_ratio = exp(.data$estimate),
      odds_ratio_low = exp(.data$conf.low),
      odds_ratio_high = exp(.data$conf.high),
      feature_group = dplyr::case_when(
        stringr::str_starts(.data$term, "pam") ~ "PAM",
        stringr::str_starts(.data$term, "pam_dist_bin") ~ "Edit position",
        .data$term == "deepcpf1_score" ~ "Guide score",
        .data$term %in% c("max_a", "max_t") ~ "Homopolymers",
        .data$term == "vardiff" ~ "Variant size",
        stringr::str_starts(.data$term, "round") ~ "Round",
        TRUE ~ "Other"
      ),
      term_label = dplyr::case_when(
        .data$term == "pamTTTC" ~ "TTTC vs TTTA",
        .data$term == "pamTTTG" ~ "TTTG vs TTTA",
        .data$term == "pam_dist_bin7-11" ~ "7-11 nt vs 1-6 nt",
        .data$term == "pam_dist_bin12-17" ~ "12-17 nt vs 1-6 nt",
        .data$term == "pam_dist_bin18-23" ~ "18-23 nt vs 1-6 nt",
        .data$term == "pam_dist_binPAM" ~ "PAM-overlap vs 1-6 nt",
        .data$term == "pam_dist_bin>23" ~ ">23 nt vs 1-6 nt",
        .data$term == "deepcpf1_score" ~ "DeepCpf1 score",
        .data$term == "max_a" ~ "max A run",
        .data$term == "max_t" ~ "max T run",
        .data$term == "vardiff" ~ "Variant length difference",
        .data$term == "roundround_2" ~ "Round 2 vs Round 1",
        TRUE ~ .data$term
      )
    ) %>%
    dplyr::filter(.data$term != "(Intercept)")
}

fit_genome_wide_extended_weighted_model <- function(design_outcomes) {
  model_data <- design_outcomes %>%
    dplyr::filter(.data$assayed_design) %>%
    dplyr::mutate(
      pam = factor(.data$pam, levels = c("TTTA", "TTTC", "TTTG")),
      pam_dist_bin = factor(
        .data$pam_dist_bin,
        levels = c("1-6", "7-11", "12-17", "18-23", "PAM", ">23")
      ),
      pam_upstream_nt = factor(.data$pam_upstream_nt, levels = c("C", "A", "G", "T")),
      dr_spacer_pair_any = factor(.data$dr_spacer_pair_any, levels = c(FALSE, TRUE))
    )

  stats::glm(
    cbind(edited_colonies, assayed_colonies - edited_colonies) ~
      pam + pam_dist_bin + deepcpf1_score + pam_upstream_nt +
      guide_gc_20 + donor_wt_gc + dr_spacer_pair_any,
    family = stats::binomial(),
    data = model_data
  )
}

tidy_genome_wide_extended_weighted_model <- function(model) {
  broom::tidy(model, conf.int = TRUE) %>%
    dplyr::mutate(
      odds_ratio = exp(.data$estimate),
      odds_ratio_low = exp(.data$conf.low),
      odds_ratio_high = exp(.data$conf.high),
      feature_group = dplyr::case_when(
        stringr::str_starts(.data$term, "pam_upstream_nt") ~ "PAM upstream base",
        stringr::str_starts(.data$term, "pam") ~ "PAM",
        stringr::str_starts(.data$term, "pam_dist_bin") ~ "Edit position",
        .data$term == "deepcpf1_score" ~ "Guide score",
        .data$term %in% c("guide_gc_20", "donor_wt_gc") ~ "GC content",
        .data$term == "dr_spacer_pair_anyTRUE" ~ "DR folding",
        TRUE ~ "Other"
      ),
      term_label = dplyr::case_when(
        .data$term == "pamTTTC" ~ "TTTC vs TTTA",
        .data$term == "pamTTTG" ~ "TTTG vs TTTA",
        .data$term == "pam_dist_bin7-11" ~ "7-11 nt vs 1-6 nt",
        .data$term == "pam_dist_bin12-17" ~ "12-17 nt vs 1-6 nt",
        .data$term == "pam_dist_bin18-23" ~ "18-23 nt vs 1-6 nt",
        .data$term == "pam_dist_binPAM" ~ "PAM-overlap vs 1-6 nt",
        .data$term == "pam_upstream_ntA" ~ "PAM 5' A vs C",
        .data$term == "pam_upstream_ntG" ~ "PAM 5' G vs C",
        .data$term == "pam_upstream_ntT" ~ "PAM 5' T vs C",
        .data$term == "deepcpf1_score" ~ "DeepCpf1 score",
        .data$term == "guide_gc_20" ~ "Guide GC (first 20 nt)",
        .data$term == "donor_wt_gc" ~ "WT donor GC",
        .data$term == "dr_spacer_pair_anyTRUE" ~ "Any DR-spacer pairing",
        TRUE ~ .data$term
      )
    ) %>%
    dplyr::filter(.data$term != "(Intercept)")
}

plot_genome_wide_stage_attrition <- function(stage_summary) {
  plot_data <- stage_summary %>%
    dplyr::filter(.data$round != "overall") %>%
    dplyr::transmute(
      round = .data$round,
      picked_designs = .data$picked_designs,
      assayed_designs = .data$assayed_designs,
      edited_designs = .data$edited_designs,
      picked_colonies = .data$picked_colonies,
      assayed_colonies = .data$assayed_colonies,
      edited_colonies = .data$edited_colonies
    ) %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of("round"),
      names_to = "metric",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      measure = dplyr::if_else(stringr::str_detect(.data$metric, "design"), "Designs", "Colonies"),
      stage = dplyr::case_when(
        stringr::str_starts(.data$metric, "picked") ~ "Picked",
        stringr::str_starts(.data$metric, "assayed") ~ "Assayed",
        TRUE ~ "Edited"
      ),
      stage = factor(.data$stage, levels = c("Picked", "Assayed", "Edited")),
      round = dplyr::recode(.data$round, round_1 = "Round 1", round_2 = "Round 2")
    )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$stage, y = .data$count, fill = .data$round)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.72) +
    ggplot2::facet_wrap(~measure, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("Round 1" = "#6D597A", "Round 2" = "#B56576")) +
    ggplot2::labs(
      x = NULL,
      y = "Count",
      fill = NULL,
      title = "Genome-wide assay attrition"
    ) +
    publication_theme(base_size = 11, legend_position = "top", aspect_ratio = NULL)
}

plot_genome_wide_rate_bars <- function(summary_table, x_col, x_label, title) {
  ggplot2::ggplot(
    summary_table,
    ggplot2::aes(x = .data[[x_col]], y = .data$design_success_rate, fill = .data[[x_col]])
  ) +
    ggplot2::geom_col(width = 0.72, colour = "black", linewidth = 0.2, show.legend = FALSE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$design_success_lower, ymax = .data$design_success_upper),
      width = 0.15
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("n=%d", .data$total_guides)),
      vjust = -0.35,
      size = 3
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.05),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = x_label,
      y = "Design success rate",
      title = title
    ) +
    publication_theme(base_size = 11, legend_position = "none", aspect_ratio = 1)
}

plot_genome_wide_pam_distance_heatmap <- function(summary_table) {
  plot_data <- summary_table %>%
    dplyr::mutate(
      pam = factor(.data$pam, levels = c("TTTC", "TTTA", "TTTG")),
      pam_dist_bin = factor(.data$pam_dist_bin, levels = c("1-6", "7-11", "12-17", "18-23", "PAM")),
      label = sprintf("%s\nn=%d", scales::percent(.data$design_success_rate, accuracy = 1), .data$total_guides)
    )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$pam_dist_bin, y = .data$pam, fill = .data$design_success_rate)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3) +
    ggplot2::scale_fill_gradientn(
      colours = c("#F1D3B3", "#D88C74", "#8C3B3B"),
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Edit position relative to PAM",
      y = "PAM",
      fill = "Design success",
      title = "PAM and edit position dominate design success"
    ) +
    publication_theme(base_size = 11, legend_position = "right", aspect_ratio = 0.8)
}

plot_genome_wide_weighted_model <- function(model_table) {
  plot_data <- model_table %>%
    dplyr::filter(.data$feature_group != "Round") %>%
    dplyr::mutate(
      term_label = factor(.data$term_label, levels = rev(.data$term_label))
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$odds_ratio,
      y = .data$term_label,
      xmin = .data$odds_ratio_low,
      xmax = .data$odds_ratio_high,
      colour = .data$feature_group
    )
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbarh(height = 0.18, linewidth = 0.5) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(
      values = c(
        "PAM" = "#B56576",
        "PAM upstream base" = "#E56B6F",
        "Edit position" = "#6D597A",
        "Guide score" = "#355070",
        "GC content" = "#2A9D8F",
        "DR folding" = "#7A8C5F",
        "Homopolymers" = "#6C9A8B",
        "Variant size" = "#BC6C25",
        "Other" = "#666666"
      )
    ) +
    ggplot2::labs(
      x = "Odds ratio (log scale)",
      y = NULL,
      colour = NULL,
      title = "Weighted multivariable model"
    ) +
    publication_theme(base_size = 11, legend_position = "top", aspect_ratio = NULL)
}

plot_genome_wide_coverage_qc <- function(colony_assays) {
  correlation <- suppressWarnings(
    stats::cor(colony_assays$mean_cov, colony_assays$target_basecov, use = "pairwise.complete.obs")
  )

  ggplot2::ggplot(
    colony_assays %>% dplyr::filter(!is.na(.data$mean_cov), !is.na(.data$target_basecov)),
    ggplot2::aes(x = .data$mean_cov, y = .data$target_basecov, colour = .data$edited_colony)
  ) +
    ggplot2::geom_point(alpha = 0.7, size = 1.8) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted", colour = "grey50") +
    ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "#355070", linewidth = 0.8) +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = "#B56576", "FALSE" = "#6C757D"),
      labels = c("FALSE" = "Not edited", "TRUE" = "Edited")
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.1,
      vjust = 1.2,
      label = sprintf("Pearson r = %.2f", correlation),
      size = 3.4
    ) +
    ggplot2::labs(
      x = "Mean genome coverage",
      y = "Target base coverage",
      colour = NULL,
      title = "Coverage QC"
    ) +
    publication_theme(base_size = 11, legend_position = "top", aspect_ratio = 1)
}

plot_multi_colony_consistency <- function(consistency_summary) {
  plot_data <- consistency_summary %>%
    dplyr::mutate(
      assayed_colonies = factor(.data$assayed_colonies),
      label = sprintf("n=%d", .data$total_designs)
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$assayed_colonies, y = .data$fraction, fill = .data$consistency)
  ) +
    ggplot2::geom_col(width = 0.72, colour = "black", linewidth = 0.2) +
    ggplot2::geom_text(
      data = plot_data %>% dplyr::distinct(.data$assayed_colonies, .data$total_designs, .data$label),
      ggplot2::aes(x = .data$assayed_colonies, y = 1.03, label = .data$label),
      inherit.aes = FALSE,
      size = 3
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "All edited" = "#355070",
        "Mixed" = "#E56B6F",
        "None edited" = "#D9D9D9"
      )
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.08),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = "Assayed colonies per design",
      y = "Fraction of designs",
      fill = NULL,
      title = "Repeated colonies do not always agree"
    ) +
    publication_theme(base_size = 11, legend_position = "top", aspect_ratio = 1)
}

plot_background_variant_burden <- function(colony_background_summary) {
  plot_data <- colony_background_summary %>%
    dplyr::mutate(
      target_edit_status = dplyr::if_else(.data$edited_colony, "Edited target", "Unedited target")
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$target_edit_status, y = .data$high_conf_non_target_count, fill = .data$target_edit_status)
  ) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.14, height = 0.02, alpha = 0.45, size = 1.1) +
    ggplot2::scale_fill_manual(values = c("Edited target" = "#355070", "Unedited target" = "#B56576")) +
    ggplot2::scale_y_continuous(breaks = 0:4, limits = c(0, 4)) +
    ggplot2::labs(
      x = NULL,
      y = "High-confidence non-target variants per colony",
      fill = NULL,
      title = "Background variant burden is sparse"
    ) +
    publication_theme(base_size = 11, legend_position = "none", aspect_ratio = 1)
}
