suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) {
    y
  } else {
    x
  }
}

read_variant_table <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("tsv", "txt")) {
    return(readr::read_tsv(path, show_col_types = FALSE, name_repair = "unique_quiet"))
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

amplicon_efficiency_from_long_table <- function(long_variant_table, guide_features_path = NULL) {
  data <- if (is.character(long_variant_table)) {
    read_variant_table(long_variant_table)
  } else {
    long_variant_table
  }

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
    guide_features <- read_variant_table(guide_features_path) %>%
      dplyr::filter(panel == "Amplicon", sequence_type == "guide") %>%
      dplyr::transmute(
        guide_id,
        target,
        variant_position_annotation,
        guide_sequence = sequence,
        guide_length = sequence_length,
        pam,
        pre_pam_nt,
        distance_from_pam = distance_from_PAM,
        target_strand = strand,
        guide_gc_pct,
        deepcpf1_score
      )

    out <- out %>%
      dplyr::left_join(guide_features, by = "guide_id")
  }

  out
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
