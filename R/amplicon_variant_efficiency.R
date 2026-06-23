suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) {
    y
  } else {
    x
  }
}

normalize_column_key <- function(x) {
  x <- iconv(as.character(x), to = "ASCII//TRANSLIT")
  x <- tolower(x)
  gsub("[^a-z0-9]+", "", x)
}

copy_column_alias <- function(data, canonical, aliases) {
  if (canonical %in% names(data)) {
    return(data)
  }

  source_index <- match(
    normalize_column_key(aliases),
    normalize_column_key(names(data)),
    nomatch = 0
  )
  source_index <- source_index[source_index > 0]

  if (length(source_index) > 0) {
    data[[canonical]] <- data[[source_index[[1]]]]
  }

  data
}

standardize_amplicon_columns <- function(data) {
  aliases <- list(
    MUTATION = c("MUTATION", "editing_class", "mutation_class"),
    frc_alt = c("frc_alt", "efficiency_pct", "allele_frequency_pct", "af"),
    frc_ref = c("frc_ref", "reference_pct", "ref_pct", "ref_af"),
    source_batch = c("source_batch"),
    source_file = c("source_file"),
    source_row = c("source_row"),
    cas_variant = c("cas_variant"),
    DR_cognate = c("DR_cognate", "direct_repeat", "dr_cognate"),
    Experiment = c("Experiment", "assay", "experiment"),
    promoter = c("promoter", "crRNA_promoter", "guide_promoter"),
    recruit = c("recruit", "donor_recruitment"),
    Guide = c("Guide", "crRNA_id", "guide_id"),
    Sample = c("Sample", "sample"),
    Timepoint = c("Timepoint", "timepoint"),
    Replicate = c("Replicate", "replicate"),
    total_gen_in_liquid = c("total_gen_in_liquid", "generations"),
    CHROM = c("CHROM", "chromosome", "chrom", "chr"),
    POS = c("POS", "position", "pos"),
    REF = c("REF", "ref"),
    ALT = c("ALT", "alt"),
    TYPE = c("TYPE", "variant_type", "type"),
    AO = c("AO", "ao", "alt_depth", "alternate_depth"),
    RO = c("RO", "ro", "ref_depth", "reference_depth"),
    DP = c("DP", "dp", "depth"),
    QUAL = c("QUAL", "qual"),
    pos_mismatch = c("pos_mismatch", "mismatch_position", "mismatch_positions"),
    mismatches = c("mismatches", "mismatch_alt", "mismatch_bases")
  )

  for (canonical in names(aliases)) {
    data <- copy_column_alias(data, canonical, aliases[[canonical]])
  }

  data
}

read_variant_table <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("tsv", "txt")) {
    return(readr::read_tsv(path, show_col_types = FALSE, name_repair = "unique_quiet"))
  }

  if (!ext %in% c("csv")) {
    stop(
      "Unsupported table format: .",
      ext,
      ". Use CSV, TSV, or TXT input.",
      call. = FALSE
    )
  }

  readr::read_csv(path, show_col_types = FALSE, name_repair = "unique_quiet")
}

write_csv_mkdir <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(data, path)
  invisible(path)
}

check_required_files <- function(paths) {
  missing_paths <- names(paths)[!file.exists(unlist(paths))]

  if (length(missing_paths) > 0) {
    stop(
      "Missing required input files: ",
      paste(missing_paths, collapse = ", "),
      call. = FALSE
    )
  }
}

ensure_column <- function(data, column, value = NA) {
  if (!column %in% names(data)) {
    data[[column]] <- value
  }

  data
}

ensure_columns <- function(data, columns) {
  for (column in columns) {
    data <- ensure_column(data, column)
  }

  data
}

standardize_guide_features <- function(guide_features) {
  aliases <- list(
    guide_id = c("guide_id", "crRNA_id", "Guide"),
    target = c("target"),
    variant_position_annotation = c("variant_position_annotation"),
    guide_sequence = c("guide_sequence", "crRNA_sequence", "sequence"),
    guide_length = c("guide_length", "crRNA_length", "sequence_length"),
    pam = c("pam"),
    pre_pam_nt = c("pre_pam_nt"),
    distance_from_pam = c("distance_from_pam", "distance_from_PAM"),
    target_strand = c("target_strand", "strand"),
    guide_gc_pct = c("guide_gc_pct", "crRNA_gc_pct"),
    deepcpf1_score = c("deepcpf1_score"),
    alternate_allele = c("alternate_allele", "ALT"),
    mutation_type = c("mutation_type"),
    donor_sequence = c("donor_sequence"),
    crRNA_promoter = c("crRNA_promoter", "promoter"),
    intended_chromosome = c("intended_chromosome", "chr"),
    intended_position = c("intended_position", "position", "mutation_position"),
    intended_match = c("intended_match", "match")
  )

  for (canonical in names(aliases)) {
    guide_features <- copy_column_alias(
      guide_features,
      canonical,
      aliases[[canonical]]
    )
  }

  required <- c("guide_id")
  missing_required <- setdiff(required, names(guide_features))
  if (length(missing_required) > 0) {
    stop(
      "Guide feature table is missing required columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  guide_features <- ensure_columns(guide_features, setdiff(names(aliases), "guide_id"))

  guide_features %>%
    dplyr::select(dplyr::all_of(names(aliases))) %>%
    dplyr::distinct()
}

amplicon_efficiency_from_long_table <- function(
  long_variant_table,
  guide_features_path = NULL
) {
  data <- if (is.character(long_variant_table)) {
    read_variant_table(long_variant_table)
  } else {
    long_variant_table
  }
  data <- standardize_amplicon_columns(data)

  required <- c("MUTATION", "frc_alt", "frc_ref")
  missing_required <- setdiff(required, names(data))
  if (length(missing_required) > 0) {
    stop(
      "Long variant table is missing required columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  data <- ensure_columns(
    data,
    c(
      "source_batch",
      "source_file",
      "source_row",
      "cas_variant",
      "DR_cognate",
      "Experiment",
      "promoter",
      "recruit",
      "Guide",
      "Sample",
      "Timepoint",
      "Replicate",
      "total_gen_in_liquid",
      "CHROM",
      "POS",
      "REF",
      "ALT",
      "TYPE",
      "AO",
      "RO",
      "DP",
      "QUAL"
    )
  )

  out <- data %>%
    dplyr::mutate(
      mutation_class = as.character(MUTATION),
      allele_frequency_pct = as.numeric(frc_alt),
      ref_pct = as.numeric(frc_ref),
      non_hdr_pct = dplyr::if_else(
        mutation_class == "HDR",
        pmax(0, 100 - allele_frequency_pct - ref_pct),
        0
      )
    ) %>%
    dplyr::transmute(
      dataset = "Amplicon_panel",
      source_batch,
      source_file,
      source_row,
      cas_variant,
      dr_cognate = DR_cognate,
      experiment = Experiment,
      promoter,
      donor_recruitment = recruit,
      guide_id = Guide,
      sample = Sample,
      timepoint = Timepoint,
      replicate = Replicate,
      total_gen_in_liquid,
      mutation_class,
      allele_frequency_pct,
      hdr_efficiency_pct = dplyr::if_else(mutation_class == "HDR", allele_frequency_pct, NA_real_),
      ref_pct,
      non_hdr_pct,
      chromosome = CHROM,
      position = POS,
      ref = REF,
      alt = ALT,
      variant_type = TYPE,
      ao = AO,
      ro = RO,
      dp = DP,
      qual = QUAL
    )

  if (!is.null(guide_features_path) && nzchar(guide_features_path)) {
    guide_features <- read_variant_table(guide_features_path)

    if (all(c("panel", "sequence_type") %in% names(guide_features))) {
      guide_features <- guide_features %>%
        dplyr::filter(.data$panel == "Amplicon", .data$sequence_type == "guide")
    }

    guide_features <- standardize_guide_features(guide_features)
    guide_features <- guide_features %>%
      dplyr::distinct(.data$guide_id, .keep_all = TRUE)

    out <- out %>%
      dplyr::left_join(guide_features, by = "guide_id")
  }

  out
}

normalize_chromosome <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sub("^CHR", "", x)
}

contains_integer <- function(text, expected) {
  mapply(
    function(value, target) {
      if (is.na(value) || is.na(target)) {
        return(FALSE)
      }

      tokens <- regmatches(
        as.character(value),
        gregexpr("[0-9]+", as.character(value))
      )[[1]]
      suppressWarnings(as.integer(tokens)) %in% as.integer(target) |> any()
    },
    text,
    expected,
    USE.NAMES = FALSE
  )
}

normalize_frequency_pct <- function(x) {
  values <- suppressWarnings(as.numeric(x))
  observed <- values[is.finite(values)]

  if (length(observed) > 0 && max(observed) <= 1) {
    values <- values * 100
  }

  values
}

first_non_missing <- function(x, default = NA) {
  observed <- x[!is.na(x)]
  if (length(observed) == 0) default else observed[[1]]
}

derive_mismatch_fields <- function(ref, alt, pos) {
  mismatch_one <- function(ref_value, alt_value, pos_value) {
    if (
      is.na(ref_value) ||
        is.na(alt_value) ||
        is.na(pos_value) ||
        ref_value == "" ||
        alt_value == ""
    ) {
      return(list(mismatches = NA_character_, positions = NA_character_, pos_mismatch = NA_integer_))
    }

    ref_bases <- strsplit(toupper(ref_value), "", fixed = TRUE)[[1]]
    alt_bases <- strsplit(toupper(alt_value), "", fixed = TRUE)[[1]]
    compared_length <- min(length(ref_bases), length(alt_bases))
    mismatch_offsets <- which(
      ref_bases[seq_len(compared_length)] != alt_bases[seq_len(compared_length)]
    ) - 1L

    if (length(mismatch_offsets) == 0) {
      mismatch_offsets <- 0L
      mismatch_bases <- toupper(alt_value)
    } else {
      mismatch_bases <- paste0(alt_bases[mismatch_offsets + 1L], collapse = "")
    }

    list(
      mismatches = mismatch_bases,
      positions = paste0("[", paste(mismatch_offsets, collapse = ","), "]"),
      pos_mismatch = as.integer(pos_value) + mismatch_offsets[[1]]
    )
  }

  derived <- mapply(
    mismatch_one,
    as.character(ref),
    as.character(alt),
    suppressWarnings(as.integer(pos)),
    SIMPLIFY = FALSE
  )

  tibble::tibble(
    derived_mismatches = vapply(derived, `[[`, character(1), "mismatches"),
    derived_positions = vapply(derived, `[[`, character(1), "positions"),
    derived_pos_mismatch = vapply(derived, `[[`, integer(1), "pos_mismatch")
  )
}

standardize_endogenous_designs <- function(design_table) {
  designs <- standardize_guide_features(design_table)
  required <- c("guide_id")
  missing_required <- required[
    vapply(required, function(column) all(is.na(designs[[column]])), logical(1))
  ]

  if (length(missing_required) > 0) {
    stop(
      "Endogenous design table is missing required values for: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  match_parts <- stringr::str_match(
    as.character(designs$intended_match),
    "^([^_]+)_([^_]+)_([^_]+)_([0-9]+)_(.+)$"
  )

  standardized <- designs %>%
    dplyr::mutate(
      crRNA_promoter = dplyr::coalesce(
        as.character(.data$crRNA_promoter),
        match_parts[, 2]
      ),
      expected_chromosome = dplyr::coalesce(
        as.character(.data$intended_chromosome),
        match_parts[, 4]
      ),
      expected_position = dplyr::coalesce(
        suppressWarnings(as.integer(.data$intended_position)),
        suppressWarnings(as.integer(match_parts[, 5]))
      ),
      expected_alternate_allele = toupper(dplyr::coalesce(
        as.character(.data$alternate_allele),
        match_parts[, 6]
      ))
    ) %>%
    dplyr::rename(crRNA_id = "guide_id")

  if (
    all(is.na(standardized$expected_position)) ||
      all(is.na(standardized$expected_alternate_allele))
  ) {
    stop(
      paste(
        "The design table needs an exact intended locus.",
        "Provide either intended_match/match or intended_chromosome (or chr),",
        "intended_position (or position), and alternate_allele (or ALT).",
        "variant_position_annotation is retained as publication metadata but",
        "is not used because it contains the PAM coordinate in the manuscript table."
      ),
      call. = FALSE
    )
  }

  standardized
}

merge_endogenous_design_features <- function(design_assignments, guide_features) {
  assignments <- standardize_guide_features(design_assignments)
  features <- standardize_guide_features(guide_features) %>%
    dplyr::distinct(.data$guide_id, .keep_all = TRUE)

  feature_columns <- c(
    "target",
    "variant_position_annotation",
    "guide_sequence",
    "guide_length",
    "pam",
    "pre_pam_nt",
    "distance_from_pam",
    "target_strand",
    "guide_gc_pct",
    "deepcpf1_score",
    "alternate_allele",
    "mutation_type",
    "donor_sequence"
  )

  assignment_loci <- assignments %>%
    dplyr::select(
      guide_id,
      crRNA_promoter,
      intended_chromosome,
      intended_position,
      intended_match,
      assignment_alternate_allele = alternate_allele
    )

  features %>%
    dplyr::select(guide_id, dplyr::all_of(feature_columns)) %>%
    dplyr::right_join(assignment_loci, by = "guide_id") %>%
    dplyr::mutate(
      alternate_allele = dplyr::coalesce(
        .data$assignment_alternate_allele,
        .data$alternate_allele
      )
    ) %>%
    dplyr::select(-"assignment_alternate_allele")
}

annotate_endogenous_amplicon_variants <- function(
  long_variant_table,
  design_table,
  intended_ao_min = 2,
  unintended_af_min = 50,
  unintended_dp_min = 4,
  min_dp = 0
) {
  variants <- if (is.character(long_variant_table)) {
    read_variant_table(long_variant_table)
  } else {
    long_variant_table
  }
  variants <- standardize_amplicon_columns(variants)

  designs <- if (is.character(design_table)) {
    read_variant_table(design_table)
  } else {
    design_table
  }
  designs <- standardize_endogenous_designs(designs)

  required <- c("Sample", "Guide", "CHROM", "POS", "REF", "ALT", "AO", "RO", "DP")
  missing_required <- setdiff(required, names(variants))
  if (length(missing_required) > 0) {
    stop(
      "Long endogenous variant table is missing required columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  variants <- ensure_columns(
    variants,
    c(
      "TYPE",
      "QUAL",
      "frc_alt",
      "frc_ref",
      "pos_mismatch",
      "positions",
      "mismatches",
      "promoter",
      "recruit",
      "Experiment",
      "Timepoint",
      "Replicate",
      "cas_variant",
      "DR_cognate",
      "total_gen_in_liquid"
    )
  )

  mismatch_fields <- derive_mismatch_fields(
    variants$REF,
    variants$ALT,
    variants$POS
  )
  variants <- dplyr::bind_cols(variants, mismatch_fields)

  designs_for_join <- designs %>%
    dplyr::rename(
      design_target = "target",
      design_variant_position = "variant_position_annotation",
      design_crRNA_sequence = "guide_sequence",
      design_crRNA_length = "guide_length",
      design_pam = "pam",
      design_pre_pam_nt = "pre_pam_nt",
      design_distance_from_pam = "distance_from_pam",
      design_strand = "target_strand",
      design_crRNA_gc_pct = "guide_gc_pct",
      design_deepcpf1_score = "deepcpf1_score",
      design_alternate_allele = "alternate_allele",
      design_mutation_type = "mutation_type",
      design_donor_sequence = "donor_sequence",
      design_crRNA_promoter = "crRNA_promoter"
    )

  variants <- variants %>%
    dplyr::mutate(
      Guide = as.character(.data$Guide),
      Sample = as.character(.data$Sample),
      promoter = dplyr::coalesce(
        as.character(.data$promoter),
        dplyr::case_when(
          as.character(.data$Experiment) %in% c("E1", "E3") ~ "SNR52",
          as.character(.data$Experiment) %in% c("E2", "E4") ~ "RPR1",
          TRUE ~ NA_character_
        )
      ),
      recruit = dplyr::coalesce(
        as.character(.data$recruit),
        dplyr::case_when(
          as.character(.data$Experiment) %in% c("E3", "E4") ~ "LexA-FHA",
          as.character(.data$Experiment) %in% c("E1", "E2") ~ "No",
          TRUE ~ NA_character_
        )
      ),
      AO = suppressWarnings(as.numeric(.data$AO)),
      RO = suppressWarnings(as.numeric(.data$RO)),
      DP = suppressWarnings(as.numeric(.data$DP)),
      mismatches = dplyr::coalesce(
        as.character(.data$mismatches),
        .data$derived_mismatches
      ),
      positions = dplyr::coalesce(
        as.character(.data$positions),
        .data$derived_positions
      ),
      pos_mismatch = dplyr::coalesce(
        suppressWarnings(as.integer(.data$pos_mismatch)),
        .data$derived_pos_mismatch
      ),
      allele_frequency_pct = dplyr::if_else(
        .data$DP > 0,
        100 * .data$AO / .data$DP,
        NA_real_
      ),
      reference_pct = dplyr::if_else(
        .data$DP > 0,
        100 * .data$RO / .data$DP,
        NA_real_
      ),
      frc_alt = .data$allele_frequency_pct,
      frc_ref = .data$reference_pct
    ) %>%
    dplyr::filter(!is.na(.data$DP), .data$DP > min_dp) %>%
    dplyr::select(-dplyr::starts_with("derived_"))

  use_promoter <- any(!is.na(designs_for_join$design_crRNA_promoter))
  annotated <- if (use_promoter) {
    variants %>%
      dplyr::left_join(
        designs_for_join,
        by = c(
          "Guide" = "crRNA_id",
          "promoter" = "design_crRNA_promoter"
        )
      )
  } else {
    variants %>%
      dplyr::left_join(designs_for_join, by = c("Guide" = "crRNA_id"))
  }

  annotated <- annotated %>%
    dplyr::mutate(
      design_matched = !is.na(.data$expected_position),
      chromosome_matches_design = normalize_chromosome(.data$CHROM) ==
        normalize_chromosome(.data$expected_chromosome),
      mismatch_position_matches_design = contains_integer(
        .data$pos_mismatch,
        .data$expected_position
      ),
      mismatch_alt_matches_design = stringr::str_detect(
        toupper(as.character(.data$mismatches)),
        stringr::fixed(.data$expected_alternate_allele)
      ),
      candidate_match = paste(
        .data$promoter,
        .data$Guide,
        .data$CHROM,
        .data$pos_mismatch,
        .data$mismatches,
        sep = "_"
      ),
      simple_variant_matches_design = suppressWarnings(as.integer(.data$POS)) ==
        .data$expected_position &
        toupper(as.character(.data$ALT)) == .data$expected_alternate_allele,
      is_reference_row = toupper(as.character(.data$TYPE)) == "REF" |
        as.character(.data$REF) == as.character(.data$ALT),
      is_substitution_call = tolower(as.character(.data$TYPE)) %in%
        c("snp", "mnp", "complex") &
        nchar(as.character(.data$REF)) == nchar(as.character(.data$ALT)),
      is_intended_variant = .data$design_matched &
        .data$chromosome_matches_design &
        !.data$is_reference_row &
        .data$is_substitution_call &
        (
          (!is.na(.data$intended_match) &
            .data$candidate_match == .data$intended_match) |
            (is.na(.data$intended_match) &
              (
                (.data$mismatch_position_matches_design &
                  .data$mismatch_alt_matches_design) |
                  .data$simple_variant_matches_design
              ))
        ),
      intended_edit_detected = .data$is_intended_variant &
        !is.na(.data$AO) &
        .data$AO >= intended_ao_min,
      unintended_variant_af50_dp4 = !.data$is_intended_variant &
        !.data$is_reference_row &
        !is.na(.data$allele_frequency_pct) &
        !is.na(.data$DP) &
        .data$allele_frequency_pct >= unintended_af_min &
        .data$DP >= unintended_dp_min,
      variant_annotation = dplyr::case_when(
        .data$is_intended_variant & .data$intended_edit_detected ~ "HDR",
        .data$is_intended_variant ~ "HDR_below_read_threshold",
        .data$is_reference_row ~ "REF",
        .data$unintended_variant_af50_dp4 ~ "unintended_AF50_DP4",
        TRUE ~ "other_variant"
      )
    )

  unmatched_guides <- annotated %>%
    dplyr::filter(!is.na(.data$Guide), !.data$design_matched) %>%
    dplyr::distinct(.data$Guide) %>%
    dplyr::pull(.data$Guide)

  if (length(unmatched_guides) > 0) {
    warning(
      "No endogenous design annotation was found for guide IDs: ",
      paste(utils::head(sort(unmatched_guides), 10), collapse = ", "),
      call. = FALSE
    )
  }

  annotated
}

build_endogenous_editing_window <- function(annotated_variants) {
  data <- standardize_amplicon_columns(annotated_variants)
  required <- c(
    "Sample",
    "variant_annotation",
    "is_intended_variant",
    "intended_edit_detected",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "TYPE",
    "AO",
    "RO",
    "DP"
  )
  missing_required <- setdiff(required, names(data))
  if (length(missing_required) > 0) {
    stop(
      "Annotated raw variant table is missing required columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  sample_metadata <- data %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::any_of(c(
          "cas_variant",
          "promoter",
          "DR_cognate",
          "recruit",
          "Guide",
          "Replicate",
          "Timepoint",
          "total_gen_in_liquid",
          "Experiment"
        )),
        dplyr::first
      ),
      expected_chromosome = first_non_missing(.data$expected_chromosome, NA_character_),
      expected_position = first_non_missing(.data$expected_position, NA_integer_),
      expected_alternate_allele = first_non_missing(
        .data$expected_alternate_allele,
        NA_character_
      ),
      .groups = "drop"
    )

  anchors <- data %>%
    dplyr::filter(.data$intended_edit_detected) %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::slice_max(.data$AO, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      Sample,
      anchor_chromosome = .data$CHROM,
      anchor_position = suppressWarnings(as.integer(.data$POS)),
      anchor_ref = .data$REF,
      reference_count = as.numeric(.data$RO)
    )

  observed_window <- data %>%
    dplyr::inner_join(anchors, by = "Sample") %>%
    dplyr::filter(
      normalize_chromosome(.data$CHROM) == normalize_chromosome(.data$anchor_chromosome),
      suppressWarnings(as.integer(.data$POS)) == .data$anchor_position,
      as.character(.data$REF) == as.character(.data$anchor_ref)
    ) %>%
    dplyr::mutate(
      editing_class = dplyr::if_else(
        .data$intended_edit_detected,
        "HDR",
        as.character(.data$TYPE)
      ),
      count = as.numeric(.data$AO),
      synthetic_row = FALSE
    )

  reference_rows <- observed_window %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ALT = .data$REF,
      TYPE = "REF",
      editing_class = "REF",
      count = .data$reference_count,
      AO = .data$reference_count,
      synthetic_row = TRUE,
      intended_edit_detected = FALSE,
      unintended_variant_af50_dp4 = FALSE,
      variant_annotation = "REF"
    )

  detected_window <- dplyr::bind_rows(observed_window, reference_rows) %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::mutate(
      editing_window_total = sum(.data$count, na.rm = TRUE),
      editing_window_pct = dplyr::if_else(
        .data$editing_window_total > 0,
        100 * .data$count / .data$editing_window_total,
        NA_real_
      )
    ) %>%
    dplyr::ungroup()

  missing_samples <- sample_metadata %>%
    dplyr::anti_join(anchors, by = "Sample")

  make_missing_rows <- function(class_name) {
    missing_samples %>%
      dplyr::transmute(
        Sample,
        cas_variant,
        promoter,
        DR_cognate,
        recruit,
        Guide,
        Replicate,
        Timepoint,
        total_gen_in_liquid,
        Experiment,
        CHROM = .data$expected_chromosome,
        POS = .data$expected_position,
        REF = NA_character_,
        ALT = if (class_name == "HDR") {
          .data$expected_alternate_allele
        } else {
          NA_character_
        },
        TYPE = class_name,
        AO = if (class_name == "HDR") 0 else 1,
        RO = 1,
        DP = 1,
        editing_class = class_name,
        count = if (class_name == "HDR") 0 else 1,
        editing_window_total = 1,
        editing_window_pct = if (class_name == "HDR") 0 else 100,
        intended_edit_detected = FALSE,
        unintended_variant_af50_dp4 = FALSE,
        variant_annotation = if (class_name == "HDR") "HDR_not_detected" else "REF",
        synthetic_row = TRUE
      )
  }

  dplyr::bind_rows(
    detected_window,
    make_missing_rows("HDR"),
    make_missing_rows("REF")
  ) %>%
    dplyr::arrange(.data$Sample, factor(.data$editing_class, levels = c("HDR", "REF")))
}

endogenous_amplicon_efficiency_table <- function(editing_window) {
  data <- standardize_amplicon_columns(editing_window)
  required <- c("editing_class", "editing_window_pct")
  missing_required <- setdiff(required, names(data))
  if (length(missing_required) > 0) {
    stop(
      "Editing-window table is missing required columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  sample_summary <- data %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::summarise(
      hdr_pct = sum(.data$editing_window_pct[.data$editing_class == "HDR"], na.rm = TRUE),
      ref_pct = sum(.data$editing_window_pct[.data$editing_class == "REF"], na.rm = TRUE),
      non_hdr_pct_sample = sum(
        .data$editing_window_pct[!.data$editing_class %in% c("HDR", "REF")],
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  data %>%
    dplyr::left_join(sample_summary, by = "Sample") %>%
    dplyr::mutate(
      efficiency_pct = .data$editing_window_pct,
      reference_pct = .data$ref_pct,
      non_hdr_pct = .data$non_hdr_pct_sample
    ) %>%
    dplyr::transmute(
      cas_variant,
      crRNA_promoter = promoter,
      direct_repeat = DR_cognate,
      donor_recruitment = recruit,
      crRNA_id = Guide,
      replicate = Replicate,
      timepoint = Timepoint,
      generations = total_gen_in_liquid,
      editing_class,
      efficiency_pct,
      reference_pct,
      non_hdr_pct,
      sample = Sample,
      chromosome = CHROM,
      position = POS,
      ref = REF,
      alt = ALT,
      ao = AO,
      ro = RO,
      dp = DP,
      intended_edit_detected,
      unintended_variant_af50_dp4
    )
}

simplify_variant_table <- function(df) {
  simplify_variant <- function(pos, ref, alt, type) {
    if (is.na(ref) || is.na(alt) || ref == "" || alt == "") {
      return(data.frame(new_POS = NA_integer_, new_REF = NA_character_, new_ALT = NA_character_, mismatch_details = NA_character_))
    }

    ref_seq <- strsplit(ref, "", fixed = TRUE)[[1]]
    alt_seq <- strsplit(alt, "", fixed = TRUE)[[1]]

    while (length(ref_seq) > 0 && length(alt_seq) > 0 && ref_seq[1] == alt_seq[1]) {
      ref_seq <- ref_seq[-1]
      alt_seq <- alt_seq[-1]
      pos <- pos + 1
    }

    while (length(ref_seq) > 0 && length(alt_seq) > 0 && tail(ref_seq, 1) == tail(alt_seq, 1)) {
      ref_seq <- head(ref_seq, -1)
      alt_seq <- head(alt_seq, -1)
    }

    simplified_ref <- paste(ref_seq, collapse = "")
    simplified_alt <- paste(alt_seq, collapse = "")
    mismatch_details <- NA_character_

    if (type == "del") {
      simplified_alt <- ""
    } else if (type == "ins") {
      simplified_ref <- substr(ref, pos - 1, pos - 1)
      simplified_alt <- paste0(simplified_ref, simplified_alt)
      pos <- pos - 1
    } else if (type %in% c("snp", "mnp", "complex")) {
      min_len <- min(length(ref_seq), length(alt_seq))
      mismatches <- which(ref_seq[seq_len(min_len)] != alt_seq[seq_len(min_len)])

      if (length(mismatches) > 0) {
        mismatch_details <- paste0(
          "Positions: ",
          paste(pos + mismatches - 1, collapse = ","),
          "; REF: ",
          paste(ref_seq[mismatches], collapse = ""),
          "; ALT: ",
          paste(alt_seq[mismatches], collapse = "")
        )
      }
    }

    data.frame(
      new_POS = as.integer(pos),
      new_REF = simplified_ref,
      new_ALT = simplified_alt,
      mismatch_details = mismatch_details
    )
  }

  simplified <- mapply(
    simplify_variant,
    df$POS,
    df$REF,
    df$ALT,
    df$TYPE,
    SIMPLIFY = FALSE
  )

  dplyr::bind_cols(df, dplyr::bind_rows(simplified))
}

complement_nuc <- function(nuc) {
  dplyr::case_when(
    nuc == "A" ~ "T",
    nuc == "T" ~ "A",
    nuc == "C" ~ "G",
    nuc == "G" ~ "C",
    TRUE ~ NA_character_
  )
}

summarise_donor_position_hdr <- function(
  variants,
  target_ref = "CTCAAGATACCCAGACCATATGAA",
  design_col = "Part6",
  timepoint_col = "Part7",
  replicate_col = "Part8",
  position_offset = 213,
  swap_designs = c(`17` = 18, `18` = 17)
) {
  working <- simplify_variant_table(variants)
  target_ref_chars <- strsplit(target_ref, "", fixed = TRUE)[[1]]

  working[[design_col]] <- as.integer(working[[design_col]])
  if (length(swap_designs) > 0) {
    swap_lookup <- as.integer(swap_designs)
    names(swap_lookup) <- names(swap_designs)
    swap_idx <- as.character(working[[design_col]]) %in% names(swap_lookup)
    working[[design_col]][swap_idx] <- swap_lookup[as.character(working[[design_col]][swap_idx])]
  }

  working %>%
    dplyr::filter(TYPE == "snp") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      expected_pos = position_offset + .data[[design_col]],
      correct_ref = target_ref_chars[.data[[design_col]] + 1],
      correct_alt = complement_nuc(correct_ref),
      is_correct_HDR = new_POS == expected_pos & new_REF == correct_ref & new_ALT == correct_alt
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is_correct_HDR) %>%
    dplyr::transmute(
      dist = .data[[design_col]],
      timepoint = .data[[timepoint_col]],
      replicate = .data[[replicate_col]],
      frc_alt = as.numeric(frc_alt),
      frc_ref = as.numeric(frc_ref),
      HDR_pct = frc_alt
    )
}

summarise_fixed_hdr <- function(
  variants,
  expected_pos = 217,
  correct_ref = "A",
  correct_alt = "T",
  guide_length_col = "guide_length",
  strain_col = "strain",
  replicate_col = "replicate",
  timepoint_col = "timepoint"
) {
  working <- simplify_variant_table(variants)

  working %>%
    dplyr::filter(TYPE == "snp") %>%
    dplyr::filter(new_POS == expected_pos, new_REF == correct_ref, new_ALT == correct_alt) %>%
    dplyr::transmute(
      Sample,
      guide_length = .data[[guide_length_col]],
      strain = .data[[strain_col]],
      replicate = .data[[replicate_col]],
      timepoint = .data[[timepoint_col]],
      frc_alt = as.numeric(frc_alt),
      frc_ref = as.numeric(frc_ref),
      HDR_pct = frc_alt
    )
}

prepare_fn_donor_variants <- function(variants) {
  variants <- variants %>%
    dplyr::mutate(
      Part4 = as.character(Part4),
      Part6 = suppressWarnings(as.integer(Part6))
    )

  cleaned <- variants %>%
    dplyr::filter(!(Part7 == "timepoint3" & Part6 %in% c(8, 9, 10, 11, 20, 21, 22, 23)))

  mapping <- tibble::tibble(
    Part3 = "repeat",
    Part4 = as.character(1:16),
    Part8 = rep(1:2, 8),
    Part6 = rep(c(8, 9, 10, 11, 20, 21, 22, 23), each = 2),
    Part7 = "timepoint3"
  )

  updated_rows <- cleaned %>%
    dplyr::filter(Part3 == "repeat", Part4 %in% mapping$Part4) %>%
    dplyr::select(-Part6, -Part7, -Part8) %>%
    dplyr::left_join(mapping, by = c("Part3", "Part4"))

  cleaned %>%
    dplyr::filter(!(Part3 == "repeat" & Part4 %in% mapping$Part4)) %>%
    dplyr::bind_rows(updated_rows)
}

prepare_enas_donor_variants <- function(variants) {
  variants %>%
    dplyr::filter(Part6 == "donor")
}

prepare_short_guide_variants <- function(variants) {
  variants %>%
    dplyr::filter(Part6 == "length") %>%
    dplyr::mutate(
      Part11 = as.character(Part11),
      Part10 = as.character(Part10),
      replicate = dplyr::if_else(Part11 %in% c("1", "2"), Part11, Part10),
      timepoint = dplyr::if_else(Part11 %in% c("1", "2"), as.character(Part12), Part11),
      guide_length = as.character(Part7),
      strain = as.character(Part8)
    )
}

prepare_donor_and_guide_length_outputs <- function(paths, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  check_required_files(paths)

  fn_donor_raw <- read_variant_table(paths$fn_donor_variants)
  enas_donor_raw <- read_variant_table(paths$enas_donor_variants)
  short_guide_raw <- read_variant_table(paths$short_guide_variants)

  fn_donor <- prepare_fn_donor_variants(fn_donor_raw) %>%
    summarise_donor_position_hdr(
      design_col = "Part6",
      timepoint_col = "Part7",
      replicate_col = "Part8",
      swap_designs = c(`17` = 18, `18` = 17)
    ) %>%
    dplyr::mutate(strain = "ys88", cas_variant = "FnCas12a")

  enas_donor <- prepare_enas_donor_variants(enas_donor_raw) %>%
    summarise_donor_position_hdr(
      design_col = "Part8",
      timepoint_col = "Part11",
      replicate_col = "Part10",
      swap_designs = c(`17` = 18, `18` = 17)
    ) %>%
    dplyr::mutate(strain = "ys85", cas_variant = "enAsCas12a")

  short_guide <- prepare_short_guide_variants(short_guide_raw) %>%
    summarise_fixed_hdr(
      expected_pos = 217,
      correct_ref = "A",
      correct_alt = "T",
      guide_length_col = "guide_length",
      strain_col = "strain",
      replicate_col = "replicate",
      timepoint_col = "timepoint"
    ) %>%
    dplyr::mutate(
      cas_variant = dplyr::case_when(
        strain == "ys85" ~ "enAsCas12a",
        strain == "ys88" ~ "FnCas12a",
        TRUE ~ strain
      )
    )

  write_csv_mkdir(fn_donor, file.path(output_dir, "fn_donor_position_hdr_from_variants.csv"))
  write_csv_mkdir(enas_donor, file.path(output_dir, "enas_donor_position_hdr_from_variants.csv"))
  write_csv_mkdir(dplyr::bind_rows(enas_donor, fn_donor), file.path(output_dir, "donor_position_hdr_from_variants.csv"))
  write_csv_mkdir(short_guide, file.path(output_dir, "short_guide_hdr_from_variants.csv"))

  list(
    fn_donor = fn_donor,
    enas_donor = enas_donor,
    donor_position = dplyr::bind_rows(enas_donor, fn_donor),
    short_guide = short_guide
  )
}
