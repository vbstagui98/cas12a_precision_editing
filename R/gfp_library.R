build_gfp_transfer_manifest <- function(paths) {
  transfer_counts_dir <- require_dir(
    file.path(paths$cas12a_project_root, "data", "Transfer_counts"),
    "Transfer_counts directory"
  )

  tibble::tribble(
    ~dataset_key, ~data_path, ~loader_type, ~colony_counts_path, ~promoter, ~donor_recruitment, ~cas_variant,
    "enas_rpr_no", file.path(transfer_counts_dir, "20241018_gfplib_enascas12a_rpr1_nolexA", "clean_stats.csv"), "direct", file.path(transfer_counts_dir, "20241018_gfplib_enascas12a_rpr1_nolexA", "colony_counts.csv"), "RPR1", "No", "enAsCas12a",
    "enas_rpr_yes", file.path(transfer_counts_dir, "20241028_gfplib_enascas12a_rpr1_lexA", "clean_stats.csv"), "direct", file.path(transfer_counts_dir, "20241028_gfplib_enascas12a_rpr1_lexA", "colony_counts.csv"), "RPR1", "Yes", "enAsCas12a",
    "enas_snr_no", file.path(transfer_counts_dir, "gfp_lib_snr53_nodonor_20240403", "clean_stats.csv"), "direct", file.path(transfer_counts_dir, "gfp_lib_snr53_nodonor_20240403", "colony_counts.csv"), "SNR52", "No", "enAsCas12a",
    "enas_snr_yes", file.path(transfer_counts_dir, "20241218_enascas12a_snr52_lexA_gfplib"), "enas_snr_yes", file.path(transfer_counts_dir, "20241218_enascas12a_snr52_lexA_gfplib", "data", "raw", "colony_counts.csv"), "SNR52", "Yes", "enAsCas12a",
    "fncas12a_rpr_no", file.path(transfer_counts_dir, "20241123_gfplib_fncas12a_RPR1_nolexA", "data", "processed", "clean_stats.csv"), "direct", file.path(transfer_counts_dir, "20241123_gfplib_fncas12a_RPR1_nolexA", "data", "raw", "colony_counts.csv"), "RPR1", "No", "FnCas12a",
    "fncas12a_rpr_yes", file.path(transfer_counts_dir, "20241204_fncas12a_lexA_RPR1_gfplib", "data", "processed", "clean_stats.csv"), "direct", file.path(transfer_counts_dir, "20241204_fncas12a_lexA_RPR1_gfplib", "data", "raw", "colony_counts.csv"), "RPR1", "Yes", "FnCas12a",
    "fncas12a_snr_no", file.path(transfer_counts_dir, "20241123_gfplib_fncas12a_SNR_nolexA", "data", "processed", "clean_stats.csv"), "direct", file.path(transfer_counts_dir, "20241123_gfplib_fncas12a_SNR_nolexA", "data", "raw", "colony_counts.csv"), "SNR52", "No", "FnCas12a",
    "fncas12a_snr_yes", file.path(transfer_counts_dir, "20241204_fncas12a_lexA_SNR_gfplib", "data", "processed", "clean_stats.csv"), "direct", file.path(transfer_counts_dir, "20241204_fncas12a_lexA_SNR_gfplib", "data", "raw", "colony_counts.csv"), "SNR52", "Yes", "FnCas12a"
  )
}

load_enas_snr_yes_dataset <- function(experiment_dir) {
  experiment_dir <- require_dir(experiment_dir, "enAs SNR52 donor-recruitment experiment")

  char_cols <- readr::cols(.default = readr::col_character())

  pop_stats <- read_csv_clean(
    file.path(experiment_dir, "data", "processed", "pop_stats.csv"),
    col_types = char_cols
  )
  plate_map <- read_csv_clean(
    file.path(experiment_dir, "data", "raw", "platemap.csv"),
    col_types = char_cols
  )
  metadata <- read_delim_clean(
    file.path(experiment_dir, "data", "raw", "guides_donor_selected_info.csv"),
    delim = ";",
    locale = readr::locale(decimal_mark = ","),
    col_types = char_cols
  )

  assert_columns(pop_stats, c("name"))
  assert_columns(plate_map, c("name", "well_position"))
  assert_columns(metadata, c("well_position"))

  pop_stats %>%
    dplyr::inner_join(plate_map, by = "name") %>%
    dplyr::inner_join(metadata, by = "well_position")
}

load_transfer_colony_counts <- function(path) {
  colony_counts <- read_csv_clean(
    path,
    col_types = readr::cols(.default = readr::col_character())
  )
  assert_columns(colony_counts, c("well_position", "true_colony_count"))

  arrange_columns <- c("well_position", intersect("image", names(colony_counts)))

  colony_counts %>%
    dplyr::mutate(
      well_position = as.character(well_position),
      true_colony_count = as.numeric(true_colony_count)
    ) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(arrange_columns))) %>%
    dplyr::group_by(well_position) %>%
    dplyr::mutate(replicate = dplyr::if_else(dplyr::row_number() %% 2L == 1L, 1L, 2L)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      colony_folder = dplyr::any_of("folder"),
      colony_image = dplyr::any_of("image"),
      folder_true_colony_count = true_colony_count
    ) %>%
    dplyr::select(
      well_position,
      replicate,
      dplyr::any_of(c("colony_folder", "colony_image")),
      folder_true_colony_count
    )
}

standardize_gfp_condition <- function(df, colony_counts_df, promoter, donor_recruitment, cas_variant) {
  if (!("timepoint" %in% names(df)) && "time" %in% names(df)) {
    df <- dplyr::rename(df, timepoint = time)
  }

  if (!("true_colony_count" %in% names(df)) && "colony_count" %in% names(df)) {
    df <- dplyr::rename(df, true_colony_count = colony_count)
  }

  if (!("editing_efficiency" %in% names(df)) && "frequency" %in% names(df)) {
    df <- dplyr::mutate(df, editing_efficiency = 1 - as.numeric(frequency))
  }

  if (!("total_gen_in_liquid" %in% names(df))) {
    df$total_gen_in_liquid <- NA_real_
  }

  if (!("replicate" %in% names(df))) {
    df$replicate <- NA_integer_
  }

  assert_columns(df, c("well_position", "timepoint", "editing_efficiency"))

  standardized_df <- df %>%
    dplyr::mutate(
      well_position = as.character(well_position),
      replicate = suppressWarnings(as.integer(replicate)),
      timepoint = normalize_transfer_timepoint(as.character(timepoint)),
      total_gen_in_liquid = parse_decimal_number(total_gen_in_liquid),
      true_colony_count_original = as.numeric(true_colony_count),
      editing_efficiency = as.numeric(editing_efficiency)
    ) %>%
    dplyr::left_join(colony_counts_df, by = c("well_position", "replicate")) %>%
    dplyr::mutate(
      true_colony_count = folder_true_colony_count,
      editing_efficiency = dplyr::if_else(editing_efficiency <= 1, editing_efficiency * 100, editing_efficiency),
      promoter = promoter,
      donor_recruitment = donor_recruitment,
      cas_variant = cas_variant
    )

  missing_folder_counts <- sum(is.na(standardized_df$folder_true_colony_count))
  if (missing_folder_counts > 0) {
    warning(
      glue::glue(
        "Missing mapped folder colony counts for {missing_folder_counts} rows in ",
        "{promoter} / {donor_recruitment} / {cas_variant}."
      ),
      call. = FALSE
    )
  }

  canonical_columns <- c(
    "name",
    "well_position",
    "replicate",
    "total_gen_in_liquid",
    "timepoint",
    "promoter",
    "donor_recruitment",
    "cas_variant",
    "guide_start",
    "cut_site",
    "strand",
    "position",
    "editing_efficiency",
    "true_colony_count",
    "true_colony_count_original",
    "folder_true_colony_count",
    "colony_folder",
    "colony_image"
  )

  missing_columns <- setdiff(canonical_columns, names(standardized_df))
  for (column_name in missing_columns) {
    standardized_df[[column_name]] <- NA
  }

  standardized_df %>%
    dplyr::select(dplyr::all_of(canonical_columns))
}

load_all_gfp_transfer_counts <- function(paths) {
  manifest <- build_gfp_transfer_manifest(paths)

  invisible(purrr::walk(manifest$data_path, require_file))
  invisible(purrr::walk(manifest$colony_counts_path, require_file))

  loaded_list <- purrr::map(
    seq_len(nrow(manifest)),
    function(i) {
      spec_row <- manifest[i, ]
      raw_df <- switch(
        spec_row$loader_type,
        direct = read_csv_clean(
          spec_row$data_path,
          col_types = readr::cols(.default = readr::col_character())
        ),
        enas_snr_yes = load_enas_snr_yes_dataset(spec_row$data_path),
        stop("Unknown loader type: ", spec_row$loader_type)
      )
      colony_df <- load_transfer_colony_counts(spec_row$colony_counts_path)
      standardize_gfp_condition(
        df = raw_df,
        colony_counts_df = colony_df,
        promoter = spec_row$promoter,
        donor_recruitment = spec_row$donor_recruitment,
        cas_variant = spec_row$cas_variant
      )
    }
  )

  names(loaded_list) <- manifest$dataset_key

  list(
    manifest = manifest,
    datasets = loaded_list,
    all_replicates = dplyr::bind_rows(loaded_list, .id = "dataset_key")
  )
}

prune_gfp_library_data <- function(df) {
  df %>%
    dplyr::select(
      well_position,
      promoter,
      donor_recruitment,
      cas_variant,
      total_gen_in_liquid,
      timepoint,
      replicate,
      guide_start,
      cut_site,
      strand,
      position,
      editing_efficiency,
      true_colony_count
    ) %>%
    dplyr::mutate(
      editing_efficiency = as.numeric(editing_efficiency),
      editing_efficiency = dplyr::if_else(editing_efficiency <= 1, 100 * editing_efficiency, editing_efficiency),
      total_gen_in_liquid = parse_decimal_number(total_gen_in_liquid),
      total_gen_in_liquid = dplyr::coalesce(total_gen_in_liquid, 0)
    ) %>%
    dplyr::filter(true_colony_count > 0)
}

average_replicates <- function(df, group_cols = c("well_position", "timepoint"), value_col = "editing_efficiency") {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(avg_value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
}

gfp_guide_order <- function() {
  c("C05", "A01", "C03", "A04", "C04", "H11", "B05", "A12", "G03", "E02", "E05", "E01")
}

prepare_gfp_condition_averages <- function(all_replicates) {
  all_replicates %>%
    prune_gfp_library_data() %>%
    dplyr::group_by(promoter, donor_recruitment, cas_variant, well_position, timepoint) %>%
    dplyr::summarise(avg_value = mean(editing_efficiency, na.rm = TRUE), .groups = "drop")
}

plot_gfp_heatmaps <- function(data, plot_title, fill_high, control_well = "A09", col_max = 100) {
  clean_data <- data %>%
    dplyr::mutate(
      condition_label = paste(promoter, "+", label_recruitment(donor_recruitment))
    )

  timepoints <- unique(clean_data$timepoint)
  output <- vector("list", length(timepoints))
  names(output) <- timepoints

  for (timepoint_name in timepoints) {
    timepoint_data <- clean_data %>%
      dplyr::filter(timepoint == timepoint_name)

    stats_data <- timepoint_data
    if (!is.null(control_well)) {
      stats_data <- stats_data %>% dplyr::filter(well_position != control_well)
    }

    stats_df <- stats_data %>%
      dplyr::group_by(condition_label) %>%
      dplyr::summarise(true_cond_avg = mean(avg_value, na.rm = TRUE), .groups = "drop")

    timepoint_data <- timepoint_data %>%
      dplyr::left_join(stats_df, by = "condition_label") %>%
      dplyr::mutate(
        y_label_final = paste0(condition_label, "\n(Avg: ", round(true_cond_avg, 1), "%)"),
        donor_recruitment = label_recruitment(donor_recruitment)
      )

    y_levels <- timepoint_data %>%
      dplyr::mutate(
        promoter = factor(promoter, levels = c("RPR1", "SNR52")),
        donor_recruitment = factor(donor_recruitment, levels = c("LexA-FHA", "No"))
      ) %>%
      dplyr::arrange(promoter, donor_recruitment, dplyr::desc(true_cond_avg)) %>%
      dplyr::pull(y_label_final) %>%
      unique()

    guide_levels <- gfp_guide_order()
    if (!is.null(control_well)) {
      guide_levels <- c(control_well, guide_levels)
    }

    timepoint_data$y_label_final <- factor(timepoint_data$y_label_final, levels = y_levels)
    timepoint_data$well_position <- factor(timepoint_data$well_position, levels = guide_levels)

    output[[timepoint_name]] <- ggplot2::ggplot(
      timepoint_data,
      ggplot2::aes(x = well_position, y = y_label_final, fill = avg_value)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.6) +
      ggplot2::scale_fill_gradient(
        low = "#F8F9FA",
        high = fill_high,
        limits = c(0, col_max),
        oob = scales::squish,
        name = "Edited %"
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        title = paste(plot_title, "-", timepoint_name),
        x = "Guide Position",
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
        axis.text.y = ggplot2::element_text(face = "bold", color = "black", size = 8),
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 9, face = "bold"),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.3, "cm")
      )
  }

  output
}

compute_gfp_viability <- function(all_replicates, control_well = "A09") {
  all_replicates %>%
    dplyr::group_by(dataset_key, timepoint) %>%
    dplyr::mutate(
      control_colony_count = true_colony_count[well_position == control_well][1],
      viability = true_colony_count / control_colony_count
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-control_colony_count)
}

summarize_first_timepoint <- function(all_replicates, timepoint_name = "scraped", control_well = "A09") {
  all_replicates %>%
    dplyr::filter(timepoint == timepoint_name) %>%
    dplyr::group_by(dataset_key, well_position, promoter, donor_recruitment, cas_variant) %>%
    dplyr::summarise(
      efficiency_avg = mean(editing_efficiency, na.rm = TRUE),
      colony_count_avg = mean(true_colony_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(control = dplyr::if_else(well_position == control_well, "Control", "Target")) %>%
    dplyr::group_by(dataset_key) %>%
    dplyr::mutate(viability = colony_count_avg / colony_count_avg[control == "Control"][1]) %>%
    dplyr::ungroup() %>%
    dplyr::select(-control)
}

plot_gfp_kinetics <- function(all_replicates, cas_variant_name) {
  filtered <- all_replicates %>%
    dplyr::filter(cas_variant == cas_variant_name) %>%
    dplyr::mutate(
      timepoint = factor(timepoint, levels = c("scraped", "6_gen", "11_gen", "20_gen")),
      donor_recruitment = label_recruitment(donor_recruitment)
    )

  averages <- filtered %>%
    dplyr::group_by(promoter, donor_recruitment, timepoint, well_position) %>%
    dplyr::summarise(avg_frc = mean(editing_efficiency, na.rm = TRUE), .groups = "drop")

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = filtered,
      ggplot2::aes(x = timepoint, y = editing_efficiency, color = well_position),
      alpha = 0.5,
      size = 1.5,
      position = ggplot2::position_jitter(width = 0.05)
    ) +
    ggplot2::geom_line(
      data = averages,
      ggplot2::aes(x = timepoint, y = avg_frc, group = well_position, color = well_position),
      linewidth = 0.8
    ) +
    ggplot2::facet_wrap(
      ~ donor_recruitment * promoter,
      labeller = ggplot2::labeller(donor_recruitment = c("LexA-FHA" = "LexA-FHA", "No" = ""))
    ) +
    ggplot2::scale_x_discrete(
      name = "Generations",
      labels = c("scraped" = "Plate", "6_gen" = "6", "11_gen" = "11", "20_gen" = "20")
    ) +
    ggplot2::labs(x = "Timepoint", y = "Efficiency (%)") +
    publication_theme(legend_position = "right", aspect_ratio = 1)
}
