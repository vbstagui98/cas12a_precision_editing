canonical_variant_columns <- function() {
  c(
    "name",
    "parent",
    "parent_frequency",
    "plate",
    "sample",
    "group",
    "var",
    "rep",
    "control",
    "prom_cas",
    "timepoint",
    "generations",
    "efficiency"
  )
}

select_variant_columns <- function(data) {
  required <- c("var", "rep", "control", "prom_cas", "timepoint", "generations", "efficiency")
  assert_columns(data, required)

  missing_columns <- setdiff(canonical_variant_columns(), names(data))
  for (column_name in missing_columns) {
    data[[column_name]] <- NA
  }

  data %>%
    dplyr::mutate(
      var = as.character(var),
      rep = as.character(rep),
      control = as.character(control),
      prom_cas = as.character(prom_cas),
      timepoint = as.character(timepoint),
      generations = as.numeric(generations),
      efficiency = as.numeric(efficiency)
    ) %>%
    dplyr::select(dplyr::all_of(canonical_variant_columns()))
}

build_gfp_on_variants <- function(paths) {
  main_on <- read_csv_clean(file.path(paths$attune_root, "20240912", "TEF_PGK1_vars_GFP_ON_20241104_results.csv")) %>%
    dplyr::mutate(efficiency = gfp * 100) %>%
    dplyr::filter(
      !var %in% "enAsCas12a",
      control == "No",
      timepoint != "T02_fridge",
      !(var == "MbCas12a" & prom_cas == "PGK1")
    ) %>%
    select_variant_columns()

  enas_repeat <- read_csv_clean(file.path(paths$attune_root, "20240918", "editing_20240918_results.csv")) %>%
    dplyr::filter(strain %in% c("ys304", "ys318"), guide_prom == "SNR52", control == "No") %>%
    dplyr::mutate(
      efficiency = gfp * 100,
      var = "enAsCas12a",
      prom_cas = dplyr::if_else(strain == "ys304", "TEF", "PGK1"),
      timepoint = timepoint,
      generations = gen
    ) %>%
    select_variant_columns()

  mbcas_pgk1 <- read_csv_clean(file.path(paths$attune_root, "20241007", "results_PGK1_MbCas12a_20241104_results.csv")) %>%
    dplyr::filter(assay == "GFP_ON", control == "No") %>%
    dplyr::mutate(
      efficiency = gfp * 100,
      prom_cas = "PGK1"
    ) %>%
    select_variant_columns()

  less_nls <- read_csv_clean(file.path(paths$attune_root, "20241004", "editing_20241009_results.csv")) %>%
    dplyr::filter(
      var == "enAsCas12a_less_NLS",
      cas_prom %in% c("PGK1", "TEF"),
      !group %in% c("T03_constitutive_fridge", "T04_constitutive_fridge", "T05_fridge")
    ) %>%
    dplyr::mutate(
      efficiency = gfp * 100,
      prom_cas = cas_prom,
      timepoint = substr(group, 1, 3),
      generations = gen
    ) %>%
    dplyr::filter(!timepoint %in% c("T05", "T06")) %>%
    select_variant_columns()

  dplyr::bind_rows(main_on, enas_repeat, mbcas_pgk1, less_nls) %>%
    dplyr::mutate(
      rep = as.character(rep),
      var = as.character(var)
    )
}

compute_cumulative_generations <- function(data) {
  data %>%
    dplyr::mutate(
      generations = parse_decimal_number(generations),
      generations = dplyr::coalesce(generations, 0),
      tp_num = dplyr::if_else(
        stringr::str_detect(timepoint, "\\d+"),
        as.numeric(stringr::str_extract(timepoint, "\\d+")),
        0
      )
    ) %>%
    dplyr::arrange(rep, prom_cas, var, control, tp_num) %>%
    dplyr::group_by(rep, prom_cas, var, control) %>%
    dplyr::mutate(
      generations = dplyr::if_else(tp_num == 0, 0, generations),
      generations = cumsum(generations)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-tp_num)
}

load_direct_repeat_timecourse <- function(paths) {
  read_csv_clean(
    file.path(paths$cas12a_project_root, "data", "Transfer", "dilutions_20250120_cognate_DR", "clean_stats.csv"),
    col_types = readr::cols(.default = readr::col_character())
  ) %>%
    dplyr::transmute(
      identifier = identifier,
      efficiency = parse_decimal_number(editing_efficiency) * 100,
      rep = as.character(replicate),
      generations = parse_decimal_number(total_gen_in_liquid),
      timepoint = normalize_transfer_timepoint(timepoint),
      prom_cas = cas_promoter,
      var = cas_version,
      control = tolower(control),
      guide_promoter = guide_promoter,
      dr_cognate = dr_cognate
    )
}

load_direct_repeat_colony_counts <- function(paths) {
  colony_counts <- read_csv_clean(file.path(paths$cognate_dr_root, "colony_counts.csv")) %>%
    tidyr::separate(identifier, into = c("strain", "guide_plasmid"), sep = "_", remove = FALSE) %>%
    dplyr::mutate(
      rep = as.character(replicate),
      true_colony_count = as.numeric(true_colony_count)
    )

  duplicated_keys <- colony_counts %>%
    dplyr::count(identifier, rep) %>%
    dplyr::filter(n > 1)

  if (nrow(duplicated_keys) > 0) {
    stop("Direct-repeat colony counts contain duplicated identifier/rep keys.", call. = FALSE)
  }

  colony_counts
}

build_direct_repeat_final_tables <- function(paths) {
  direct_repeat <- load_direct_repeat_timecourse(paths)
  colony_counts <- load_direct_repeat_colony_counts(paths)

  with_colony_counts <- direct_repeat %>%
    dplyr::left_join(
      colony_counts %>%
        dplyr::select(identifier, rep, strain, guide_plasmid, true_colony_count),
      by = c("identifier", "rep")
    )

  missing_counts <- sum(is.na(with_colony_counts$true_colony_count))
  if (missing_counts > 0) {
    warning(
      sprintf("Missing direct-repeat colony counts for %d time-course rows.", missing_counts),
      call. = FALSE
    )
  }

  list(
    raw = direct_repeat,
    with_colony_counts = with_colony_counts,
    syn1 = with_colony_counts %>% dplyr::filter(guide_promoter == "syn1")
  )
}

guide_promoter_strain_metadata <- function() {
  tibble::tribble(
    ~strain, ~var,          ~assay,
    "ys286", "enAsCas12a",  "OFF",
    "ys289", "FnCas12a",    "OFF",
    "ys318", "enAsCas12a",  "ON",
    "ys321", "FnCas12a",    "ON"
  )
}

guide_promoter_plate_guide_map <- function() {
  tibble::tribble(
    ~guide_id, ~guide_promoter, ~control, ~guide_role,
    "A68",     "RPR1",          "No",     "targeting",
    "A69",     "RPR1",          "No",     "targeting",
    "A70",     "SCR1",          "No",     "targeting",
    "A14",     "SCR1",          "No",     "targeting",
    "A71",     "SYN1",          "No",     "targeting",
    "A26",     "SYN1",          "No",     "targeting",
    "p136",    "SNR52",         "No",     "targeting",
    "A15",     "SNR52",         "No",     "targeting",
    "p120",    "SNR52",         "Yes",    "control"
  )
}

build_guide_promoter_timecourse <- function(paths) {
  strain_meta <- guide_promoter_strain_metadata()

  read_csv_clean(file.path(paths$attune_root, "20240918", "editing_20240918_results.csv")) %>%
    dplyr::filter(
      strain %in% strain_meta$strain,
      timepoint != "T03_fridge"
    ) %>%
    dplyr::left_join(strain_meta, by = "strain") %>%
    dplyr::mutate(
      rep = as.character(rep),
      efficiency = dplyr::if_else(assay == "ON", gfp * 100, (1 - gfp) * 100),
      prom_cas = "PGK1",
      guide_promoter = guide_prom,
      generations = gen
    ) %>%
    dplyr::select(
      name,
      sample,
      strain,
      rep,
      guide_promoter,
      control,
      timepoint,
      generations,
      efficiency,
      var,
      assay,
      prom_cas
    )
}

build_pgk1_guide_promoter_plate_counts <- function(paths) {
  strain_meta <- guide_promoter_strain_metadata()
  guide_map <- guide_promoter_plate_guide_map()

  annotation <- read_excel_clean(
    file.path(paths$cas12a_project_root, "plate_counts", "annotate_plate_images.xlsx")
  )

  exclusion_column <- names(annotation)[names(annotation) %in% c("6", "col_6")]
  if (length(exclusion_column) == 1) {
    annotation <- annotation %>%
      dplyr::filter(is.na(.data[[exclusion_column]]) | .data[[exclusion_column]] != "exclude")
  }

  colony_counts <- read_csv_clean(
    file.path(paths$cas12a_project_root, "plate_counts", "plate_counts_GFP", "predict5", "colony_counts.csv")
  )

  plate_counts <- annotation %>%
    dplyr::mutate(file = paste0(file, ".jpg")) %>%
    dplyr::left_join(colony_counts, by = c("file" = "image_name")) %>%
    dplyr::rename(strain = ys, guide_id = guide) %>%
    dplyr::filter(
      strain %in% strain_meta$strain,
      folder == "20240916",
      !is.na(colonies_cutoff)
    ) %>%
    dplyr::mutate(
      rep = as.character(rep),
      folder = as.character(folder)
    ) %>%
    dplyr::left_join(strain_meta, by = "strain") %>%
    dplyr::left_join(guide_map, by = "guide_id")

  if (any(is.na(plate_counts$guide_promoter))) {
    missing_guides <- sort(unique(plate_counts$guide_id[is.na(plate_counts$guide_promoter)]))
    stop(
      sprintf("Missing PGK1 guide-promoter mappings for guide IDs: %s", paste(missing_guides, collapse = ", ")),
      call. = FALSE
    )
  }

  duplicated_keys <- plate_counts %>%
    dplyr::count(strain, rep, guide_promoter, control) %>%
    dplyr::filter(n > 1)

  if (nrow(duplicated_keys) > 0) {
    stop("PGK1 plate counts contain duplicated strain/rep/guide_promoter/control keys.", call. = FALSE)
  }

  control_counts <- plate_counts %>%
    dplyr::filter(control == "Yes") %>%
    dplyr::transmute(
      strain,
      rep,
      control_colony_count = colonies_cutoff,
      control_count_file = file
    )

  plate_counts %>%
    dplyr::left_join(control_counts, by = c("strain", "rep")) %>%
    dplyr::mutate(relative_colony_count = colonies_cutoff / control_colony_count)
}

build_guide_promoter_final_tables <- function(paths) {
  guide_promoter_data <- build_guide_promoter_timecourse(paths)
  plate_counts <- build_pgk1_guide_promoter_plate_counts(paths)

  with_colony_counts <- guide_promoter_data %>%
    dplyr::left_join(
      plate_counts %>%
        dplyr::select(
          strain,
          rep,
          guide_promoter,
          control,
          guide_id,
          guide_role,
          folder,
          file,
          colonies_cutoff,
          control_colony_count,
          control_count_file,
          relative_colony_count
        ) %>%
        dplyr::rename(
          colony_count_folder = folder,
          colony_count_file = file
        ),
      by = c("strain", "rep", "guide_promoter", "control")
    )

  missing_counts <- sum(is.na(with_colony_counts$colonies_cutoff))
  if (missing_counts > 0) {
    warning(
      sprintf("Missing PGK1 guide-promoter colony counts for %d time-course rows.", missing_counts),
      call. = FALSE
    )
  }

  list(
    base = guide_promoter_data,
    plate_counts = plate_counts,
    with_colony_counts = with_colony_counts
  )
}

build_gfp_off_variants <- function(paths) {
  base_off <- read_csv_clean(file.path(paths$downloads_cas_variants_root, "clean_stats.csv")) %>%
    dplyr::rename(
      efficiency = editing_efficiency,
      rep = replicate,
      generations = n_gen,
      prom_cas = cas_promoter,
      var = cas_version
    ) %>%
    dplyr::filter(prom_cas != "GAL") %>%
    dplyr::mutate(
      efficiency = efficiency * 100,
      control = tolower(control)
    ) %>%
    compute_cumulative_generations() %>%
    dplyr::filter(!(var == "MbCas12a" & prom_cas == "PGK1")) %>%
    select_variant_columns()

  mbcas_off <- read_csv_clean(file.path(paths$attune_root, "20241007", "results_PGK1_MbCas12a_20241104_results.csv")) %>%
    dplyr::filter(assay == "GFP_OFF", control == "No") %>%
    dplyr::mutate(
      efficiency = (1 - gfp) * 100,
      prom_cas = "PGK1",
      control = "no"
    ) %>%
    select_variant_columns()

  dr_data <- load_direct_repeat_timecourse(paths)

  redone_tef <- dr_data %>%
    dplyr::filter(
      guide_promoter == "SNR52",
      control == "no",
      dr_cognate == "non_cognate" | var %in% c("AsCas12a", "enAsCas12a")
    ) %>%
    dplyr::select(-guide_promoter, -dr_cognate) %>%
    select_variant_columns()

  off_updated <- dplyr::bind_rows(
    base_off %>% dplyr::filter(prom_cas != "TEF"),
    mbcas_off,
    redone_tef
  ) %>%
    dplyr::distinct()

  list(
    base = dplyr::bind_rows(base_off, mbcas_off) %>% dplyr::distinct(),
    direct_repeat = dr_data,
    updated = off_updated
  )
}

plot_direct_repeat_by_guide_promoter <- function(data, cas_promoter, exclude = character(), nrows = 1) {
  filtered <- data %>%
    dplyr::filter(
      prom_cas == cas_promoter,
      control == "no",
      !var %in% exclude
    )

  summary_df <- filtered %>%
    dplyr::group_by(var, guide_promoter, dr_cognate, timepoint) %>%
    dplyr::summarise(
      generations = mean(generations, na.rm = TRUE),
      efficiency = mean(efficiency, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot2::ggplot(
    filtered,
    ggplot2::aes(
      x = generations,
      y = efficiency,
      color = var,
      shape = dr_cognate,
      alpha = dr_cognate
    )
  ) +
    ggplot2::geom_line(
      data = summary_df,
      ggplot2::aes(group = interaction(var, dr_cognate)),
      linewidth = 0.7
    ) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(~guide_promoter, nrow = nrows) +
    ggplot2::scale_alpha_manual(values = c("cognate" = 1, "non_cognate" = 0.5)) +
    ggplot2::scale_color_manual(values = variant_palette()) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::labs(x = "Generations", y = "Efficiency (%)") +
    publication_theme(legend_position = "top", aspect_ratio = 1)
}

plot_variant_timecourse <- function(data, cas_promoter, exclude = character(), facet = FALSE, direct_repeat = FALSE, nrows = 1) {
  filtered <- data %>%
    dplyr::filter(prom_cas == cas_promoter, !var %in% exclude)

  if (direct_repeat) {
    assert_columns(filtered, c("dr_cognate"))

    summary_df <- filtered %>%
      dplyr::group_by(var, dr_cognate, timepoint) %>%
      dplyr::summarise(
        generations = mean(generations, na.rm = TRUE),
        efficiency = mean(efficiency, na.rm = TRUE),
        .groups = "drop"
      )

    plot <- ggplot2::ggplot(
      filtered,
      ggplot2::aes(
        x = generations,
        y = efficiency,
        color = var,
        shape = dr_cognate,
        alpha = dr_cognate
      )
    ) +
      ggplot2::geom_line(
        data = summary_df,
        ggplot2::aes(group = interaction(var, dr_cognate), linetype = dr_cognate),
        linewidth = 0.7
      ) +
      ggplot2::geom_point(size = 1.3) +
      ggplot2::scale_alpha_manual(values = c("cognate" = 1, "non_cognate" = 0.5))
  } else {
    summary_df <- filtered %>%
      dplyr::group_by(var, timepoint) %>%
      dplyr::summarise(
        generations = mean(generations, na.rm = TRUE),
        efficiency = mean(efficiency, na.rm = TRUE),
        .groups = "drop"
      )

    plot <- ggplot2::ggplot(filtered, ggplot2::aes(generations, efficiency, color = var)) +
      ggplot2::geom_line(
        data = summary_df,
        ggplot2::aes(group = var),
        linewidth = 0.7
      ) +
      ggplot2::geom_point(size = 2)
  }

  plot <- plot +
    ggplot2::scale_color_manual(values = variant_palette()) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::labs(x = "Generations", y = "Efficiency (%)") +
    publication_theme(legend_position = "top", aspect_ratio = 1)

  if (facet) {
    plot <- plot + ggplot2::facet_wrap(~var, nrow = nrows) +
      ggplot2::theme(strip.text = ggplot2::element_blank())
  }

  plot
}

plot_codon_optimized_comparison <- function(data, cas_promoter) {
  plot_data <- data %>%
    dplyr::filter(prom_cas == cas_promoter)

  summary_df <- plot_data %>%
    dplyr::group_by(parent_var, var, timepoint) %>%
    dplyr::summarise(
      generations = mean(generations, na.rm = TRUE),
      efficiency = mean(efficiency, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot2::ggplot(plot_data, ggplot2::aes(generations, efficiency, color = var)) +
    ggplot2::geom_line(data = summary_df, ggplot2::aes(group = var), linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~parent_var, nrow = 1) +
    ggplot2::scale_color_manual(values = variant_palette()) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::labs(x = "Generations", y = "Efficiency (%)") +
    publication_theme(legend_position = "right", aspect_ratio = 1)
}

plot_guide_promoter_timecourse <- function(data, assay_name) {
  filtered <- data %>%
    dplyr::filter(assay == assay_name, control == "No")

  summary_df <- filtered %>%
    dplyr::group_by(var, guide_promoter, timepoint) %>%
    dplyr::summarise(
      generations = mean(generations, na.rm = TRUE),
      efficiency = mean(efficiency, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot2::ggplot(filtered, ggplot2::aes(generations, efficiency, color = var)) +
    ggplot2::geom_line(data = summary_df, ggplot2::aes(group = var), linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~guide_promoter, nrow = 1) +
    ggplot2::scale_color_manual(values = variant_palette()) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::labs(x = "Generations", y = "Efficiency (%)") +
    publication_theme(legend_position = "right", aspect_ratio = 1)
}
