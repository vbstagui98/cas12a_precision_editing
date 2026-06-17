figure_table_columns <- function() {
  c(
    "figure",
    "panel",
    "assay",
    "cas_variant",
    "cas_promoter",
    "guide_promoter",
    "direct_repeat",
    "donor_recruitment",
    "guide_id",
    "well_position",
    "strain",
    "replicate",
    "timepoint",
    "generations",
    "position",
    "guide_length",
    "editing_class",
    "efficiency_pct",
    "reference_pct",
    "non_hdr_pct",
    "colony_count",
    "control_colony_count",
    "viability",
    "experiment",
    "control",
    "source_table"
  )
}

empty_figure_table <- function(n) {
  output <- tibble::tibble(.rows = n)
  for (column_name in figure_table_columns()) {
    output[[column_name]] <- NA_character_
  }
  output
}

as_character_column <- function(x) {
  if (inherits(x, "factor")) {
    return(as.character(x))
  }
  if (inherits(x, "Date") || inherits(x, "POSIXt")) {
    return(as.character(x))
  }
  as.character(x)
}

make_figure_table <- function(data, figure, panel, source_table, mapping = list()) {
  output <- empty_figure_table(nrow(data))
  output$figure <- figure
  output$panel <- panel
  output$source_table <- source_table

  for (target_column in names(mapping)) {
    source_column <- mapping[[target_column]]
    if (
      length(source_column) != 1 ||
        !source_column %in% names(data) ||
        !target_column %in% figure_table_columns()
    ) {
      next
    }
    output[[target_column]] <- as_character_column(data[[source_column]])
  }

  finalize_figure_table(output)
}

numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}

missing_text <- function(x) {
  is.na(x) | !nzchar(as.character(x))
}

finalize_figure_table <- function(data) {
  output <- data

  colony_count <- numeric_or_na(output$colony_count)
  control_colony_count <- numeric_or_na(output$control_colony_count)
  viability <- numeric_or_na(output$viability)

  inferred_control <- !is.na(colony_count) & !is.na(viability) & viability != 0 & is.na(control_colony_count)
  control_colony_count[inferred_control] <- colony_count[inferred_control] / viability[inferred_control]

  inferred_viability <- !is.na(colony_count) & !is.na(control_colony_count) & control_colony_count != 0 & is.na(viability)
  viability[inferred_viability] <- colony_count[inferred_viability] / control_colony_count[inferred_viability]

  output$colony_count <- ifelse(is.na(colony_count), NA_character_, as.character(colony_count))
  output$control_colony_count <- ifelse(is.na(control_colony_count), NA_character_, as.character(control_colony_count))
  output$viability <- ifelse(is.na(viability), NA_character_, as.character(viability))

  output %>%
    dplyr::select(dplyr::all_of(figure_table_columns()))
}

read_result_table <- function(paths, ...) {
  read_csv_clean(file.path(paths$results_dir, ...))
}

first_existing_result_table <- function(paths, candidates) {
  for (candidate in candidates) {
    path <- file.path(paths$results_dir, candidate)
    if (file.exists(path)) {
      return(read_csv_clean(path))
    }
  }

  stop(
    "None of the expected result tables exists: ",
    paste(file.path(paths$results_dir, candidates), collapse = ", "),
    call. = FALSE
  )
}

export_publication_figure_tables <- function(paths, output_dir = NULL) {
  output_dir <- output_dir %||% file.path(paths$results_dir, "publication_figure_tables")
  output_dir <- ensure_dir(output_dir)

  tables <- list(
    figure_1 = build_figure_1_table(paths),
    figure_2 = build_figure_2_table(paths),
    figure_3 = build_figure_3_table(paths),
    figure_4 = build_figure_4_table(paths)
  )

  purrr::iwalk(
    tables,
    function(table_data, table_name) {
      write_csv_checked(table_data, file.path(output_dir, paste0(table_name, "_source_data.csv")))
    }
  )

  summary <- purrr::imap_dfr(
    tables,
    function(table_data, table_name) {
      table_data %>%
        dplyr::count(figure, panel, source_table, name = "rows") %>%
        dplyr::mutate(file = paste0(table_name, "_source_data.csv"), .before = 1)
    }
  )
  write_csv_checked(summary, file.path(output_dir, "source_data_summary.csv"))

  invisible(list(output_dir = output_dir, tables = tables, summary = summary))
}

build_figure_1_table <- function(paths) {
  paper_table <- try_build_figure_1_table_from_paper(paths)
  if (!is.null(paper_table)) {
    return(paper_table)
  }

  gfp_on <- read_result_table(paths, "04_cas_variants_overview", "tables", "gfp_on_variants.csv")
  gfp_off <- read_result_table(paths, "04_cas_variants_overview", "tables", "gfp_off_variants.csv")
  direct_repeat <- read_result_table(paths, "04_cas_variants_overview", "tables", "gfp_off_direct_repeat.csv")
  guide_promoter <- read_result_table(paths, "04_cas_variants_overview", "tables", "guide_promoter_data_with_colony_counts.csv")

  gfp_on_plotted <- gfp_on %>%
    dplyr::filter(
      (prom_cas == "TEF" & !var %in% c("enAsCas12a_less_NLS", "FnCas12a(CO)", "LbCas12a(CO)")) |
        (prom_cas == "PGK1" & !var %in% c("FnCas12a(CO)", "LbCas12a(CO)"))
    )

  gfp_off_plotted <- gfp_off %>%
    dplyr::filter(
      (prom_cas == "TEF" & !var %in% c("FnCas12a(CO)", "LbCas12a(CO)")) |
        (prom_cas == "PGK1" & !var %in% c("enAsCas12a_less_NLS", "FnCas12a(CO)", "LbCas12a(CO)"))
    )

  direct_repeat_plotted <- direct_repeat %>%
    dplyr::filter(
      prom_cas == "TEF",
      !var %in% c("enAsCas12a_less_NLS", "FnCas12a(CO)", "LbCas12a(CO)")
    )

  dplyr::bind_rows(
    make_figure_table(
      gfp_on_plotted,
      figure = "Figure 1",
      panel = "GFP ON Cas-variant time course",
      source_table = "results/04_cas_variants_overview/tables/gfp_on_variants.csv",
      mapping = list(
        assay = "control",
        cas_variant = "var",
        cas_promoter = "prom_cas",
        replicate = "rep",
        timepoint = "timepoint",
        generations = "generations",
        efficiency_pct = "efficiency",
        sample = "sample",
        control = "control"
      )
    ) %>% dplyr::mutate(assay = "GFP_ON"),
    make_figure_table(
      gfp_off_plotted,
      figure = "Figure 1",
      panel = "GFP OFF Cas-variant time course",
      source_table = "results/04_cas_variants_overview/tables/gfp_off_variants.csv",
      mapping = list(
        cas_variant = "var",
        cas_promoter = "prom_cas",
        replicate = "rep",
        timepoint = "timepoint",
        generations = "generations",
        efficiency_pct = "efficiency",
        sample = "sample",
        control = "control"
      )
    ) %>% dplyr::mutate(assay = "GFP_OFF"),
    make_figure_table(
      direct_repeat_plotted,
      figure = "Figure 1",
      panel = "Direct-repeat comparison",
      source_table = "results/04_cas_variants_overview/tables/gfp_off_direct_repeat.csv",
      mapping = list(
        assay = "control",
        cas_variant = "var",
        cas_promoter = "prom_cas",
        guide_promoter = "guide_promoter",
        direct_repeat = "dr_cognate",
        strain = "strain",
        guide_id = "guide_plasmid",
        replicate = "rep",
        timepoint = "timepoint",
        generations = "generations",
        efficiency_pct = "efficiency",
        colony_count = "true_colony_count",
        control = "control"
      )
    ) %>% dplyr::mutate(assay = "GFP_OFF"),
    make_figure_table(
      guide_promoter,
      figure = "Figure 1",
      panel = "Guide-promoter comparison",
      source_table = "results/04_cas_variants_overview/tables/guide_promoter_data_with_colony_counts.csv",
      mapping = list(
        assay = "assay",
        cas_variant = "var",
        cas_promoter = "prom_cas",
        guide_id = "guide_id",
        guide_promoter = "guide_promoter",
        strain = "strain",
        replicate = "rep",
        timepoint = "timepoint",
        generations = "generations",
        efficiency_pct = "efficiency",
        colony_count = "colonies_cutoff",
        control_colony_count = "control_colony_count",
        relative_colony_count = "relative_colony_count",
        sample = "sample",
        control = "control"
      )
    )
  )
}

try_build_figure_1_table_from_paper <- function(paths) {
  root <- paths$paper_notebooks_root %||% ""
  required_files <- c(
    "processed_data/GFP_OFF_vars_colonies_20260311.csv",
    "processed_data/GFP_ON_vars_colonies_20260311.csv",
    "processed_data/PGK1_proms_eff_colonies_20260311.csv",
    "processed_data/DR_data_processed_20260312.csv",
    "processed_data/DR_data_gfp_on_cognate_processed_20260506.csv",
    "processed_data/DR_data_gfp_on_viability_bridge_normalized_20260512.csv",
    "processed_data/annotated_stats_syn1_dr_pgk1_20260505.csv"
  )

  if (!nzchar(root) || !all(file.exists(file.path(root, required_files)))) {
    return(NULL)
  }

  read_paper_csv <- function(relative_path) {
    read_csv_clean(file.path(root, relative_path))
  }

  normalize_timepoint <- function(data) {
    data %>%
      dplyr::mutate(
        timepoint = dplyr::case_when(
          timepoint == "T01" ~ "scraped",
          timepoint == "T02" ~ "6_gen",
          timepoint == "T03" ~ "11_gen",
          timepoint == "T04" ~ "20_gen",
          TRUE ~ as.character(timepoint)
        )
      )
  }

  gfp_off <- read_paper_csv("processed_data/GFP_OFF_vars_colonies_20260311.csv") %>%
    normalize_timepoint() %>%
    dplyr::mutate(dr_cognate = "As_DR")

  gfp_on <- read_paper_csv("processed_data/GFP_ON_vars_colonies_20260311.csv") %>%
    normalize_timepoint() %>%
    dplyr::mutate(dr_cognate = "As_DR")

  guide_promoter_as_dr <- read_paper_csv("processed_data/PGK1_proms_eff_colonies_20260311.csv") %>%
    normalize_timepoint() %>%
    dplyr::filter(control == "No") %>%
    dplyr::mutate(
      guide_promoter = dplyr::recode(guide_promoter, SYN1 = "syn1"),
      dr_cognate = "As_DR"
    )

  guide_promoter_syn1_cognate <- read_paper_csv("processed_data/annotated_stats_syn1_dr_pgk1_20260505.csv") %>%
    dplyr::mutate(
      timepoint = dplyr::case_when(
        timepoint == "scraped" ~ "scraped",
        timepoint == "timepoint_1" ~ "6_gen",
        timepoint == "timepoint_2" ~ "11_gen",
        timepoint == "timepoint_3" ~ "20_gen",
        TRUE ~ as.character(timepoint)
      ),
      rep = as.integer(replicate),
      generations = as.numeric(total_gen_in_liquid),
      prom_cas = as.character(cas_promoter),
      var = as.character(cas_version),
      control = "No",
      guide_promoter = as.character(guide_promoter),
      assay = dplyr::case_when(
        assay == "GFP_OFF" ~ "OFF",
        assay == "GFP_ON" ~ "ON",
        TRUE ~ as.character(assay)
      ),
      efficiency = 100 * as.numeric(editing_efficiency),
      colonies_cutoff = as.numeric(colony_count),
      dr_cognate = dplyr::if_else(dr_cognate == "non_cognate", "As_DR", "cognate")
    ) %>%
    dplyr::group_by(prom_cas, guide_promoter, assay) %>%
    dplyr::mutate(
      control_colony_count = max(colonies_cutoff, na.rm = TRUE),
      viability = colonies_cutoff / control_colony_count
    ) %>%
    dplyr::ungroup()

  direct_repeat_off <- read_paper_csv("processed_data/DR_data_processed_20260312.csv") %>%
    dplyr::mutate(
      dr_cognate = dplyr::if_else(
        dr_cognate == "non_cognate" | var %in% c("AsCas12a", "enAsCas12a"),
        "As_DR",
        "cognate"
      )
    )

  direct_repeat_on_cognate <- read_paper_csv("processed_data/DR_data_gfp_on_cognate_processed_20260506.csv") %>%
    dplyr::mutate(
      dr_cognate = dplyr::if_else(
        dr_cognate == "non_cognate" | var %in% c("AsCas12a", "enAsCas12a"),
        "As_DR",
        "cognate"
      )
    )

  bridge_normalized <- read_paper_csv("processed_data/DR_data_gfp_on_viability_bridge_normalized_20260512.csv")

  gfp_on_plotted <- gfp_on %>%
    dplyr::filter(
      (prom_cas == "TEF" & !var %in% c("enAsCas12a_less_NLS", "FnCas12a(CO)", "LbCas12a(CO)")) |
        (prom_cas == "PGK1" & !var %in% c("FnCas12a(CO)", "LbCas12a(CO)"))
    )

  gfp_off_plotted <- gfp_off %>%
    dplyr::filter(
      (prom_cas == "TEF" & !var %in% c("FnCas12a(CO)", "LbCas12a(CO)")) |
        (prom_cas == "PGK1" & !var %in% c("enAsCas12a_less_NLS", "FnCas12a(CO)", "LbCas12a(CO)"))
    )

  guide_promoter <- dplyr::bind_rows(
    guide_promoter_as_dr %>%
      dplyr::select(rep, guide_promoter, control, timepoint, generations, efficiency, var, assay, prom_cas, viability, dr_cognate, colonies_cutoff, control_colony_count),
    guide_promoter_syn1_cognate %>%
      dplyr::select(rep, guide_promoter, control, timepoint, generations, efficiency, var, assay, prom_cas, viability, dr_cognate, colonies_cutoff, control_colony_count)
  )

  direct_repeat_off_plotted <- direct_repeat_off %>%
    dplyr::filter(
      prom_cas == "TEF",
      control == "no",
      !var %in% c("enAsCas12a_less_NLS", "FnCas12a(CO)", "LbCas12a(CO)")
    )

  direct_repeat_on_plotted <- direct_repeat_on_cognate %>%
    dplyr::filter(
      prom_cas == "TEF",
      control == "no",
      guide_promoter == "SNR52",
      !var %in% c("AsCas12a", "enAsCas12a")
    )

  dplyr::bind_rows(
    make_figure_table(
      gfp_on_plotted,
      figure = "Figure 1",
      panel = "GFP ON Cas-variant time course",
      source_table = "processed_data/GFP_ON_vars_colonies_20260311.csv",
      mapping = list(
        cas_variant = "var",
        cas_promoter = "prom_cas",
        direct_repeat = "dr_cognate",
        replicate = "rep",
        timepoint = "timepoint",
        generations = "generations",
        efficiency_pct = "efficiency",
        colony_count = "colonies_cutoff",
        control_colony_count = "max_colony_count",
        viability = "viability",
        sample = "sample",
        control = "control"
      )
    ) %>% dplyr::mutate(assay = "GFP_ON"),
    make_figure_table(
      gfp_off_plotted,
      figure = "Figure 1",
      panel = "GFP OFF Cas-variant time course",
      source_table = "processed_data/GFP_OFF_vars_colonies_20260311.csv",
      mapping = list(
        cas_variant = "var",
        cas_promoter = "prom_cas",
        direct_repeat = "dr_cognate",
        replicate = "rep",
        timepoint = "timepoint",
        generations = "generations",
        efficiency_pct = "efficiency",
        colony_count = "true_colony_count",
        control_colony_count = "max_colony_count",
        viability = "viability",
        control = "control"
      )
    ) %>% dplyr::mutate(assay = "GFP_OFF"),
    make_figure_table(
      direct_repeat_off_plotted,
      figure = "Figure 1",
      panel = "GFP OFF direct-repeat comparison",
      source_table = "processed_data/DR_data_processed_20260312.csv",
      mapping = list(
        cas_variant = "var",
        cas_promoter = "prom_cas",
        guide_promoter = "guide_promoter",
        direct_repeat = "dr_cognate",
        replicate = "rep",
        timepoint = "timepoint",
        generations = "generations",
        efficiency_pct = "efficiency",
        colony_count = "true_colony_count",
        viability = "viability",
        control = "control"
      )
    ) %>% dplyr::mutate(assay = "GFP_OFF"),
    make_figure_table(
      direct_repeat_on_plotted,
      figure = "Figure 1",
      panel = "GFP ON cognate direct-repeat kinetics",
      source_table = "processed_data/DR_data_gfp_on_cognate_processed_20260506.csv",
      mapping = list(
        cas_variant = "var",
        cas_promoter = "prom_cas",
        guide_promoter = "guide_promoter",
        direct_repeat = "dr_cognate",
        replicate = "rep",
        timepoint = "timepoint",
        generations = "generations",
        efficiency_pct = "efficiency",
        colony_count = "true_colony_count",
        viability = "viability",
        control = "control"
      )
    ) %>% dplyr::mutate(assay = "GFP_ON"),
    make_figure_table(
      bridge_normalized,
      figure = "Figure 1",
      panel = "GFP ON bridge-normalized direct-repeat viability",
      source_table = "processed_data/DR_data_gfp_on_viability_bridge_normalized_20260512.csv",
      mapping = list(
        cas_variant = "var",
        cas_promoter = "prom_cas",
        guide_promoter = "guide_promoter",
        direct_repeat = "dr_cognate",
        replicate = "rep",
        timepoint = "timepoint",
        generations = "generations",
        efficiency_pct = "efficiency",
        colony_count = "adjusted_colony_count",
        viability = "viability",
        variant = "source_dataset",
        control = "control"
      )
    ) %>% dplyr::mutate(assay = "GFP_ON"),
    make_figure_table(
      guide_promoter,
      figure = "Figure 1",
      panel = "Guide-promoter comparison",
      source_table = "processed_data/PGK1_proms_eff_colonies_20260311.csv; processed_data/annotated_stats_syn1_dr_pgk1_20260505.csv",
      mapping = list(
        assay = "assay",
        cas_variant = "var",
        cas_promoter = "prom_cas",
        guide_promoter = "guide_promoter",
        direct_repeat = "dr_cognate",
        replicate = "rep",
        timepoint = "timepoint",
        generations = "generations",
        efficiency_pct = "efficiency",
        colony_count = "colonies_cutoff",
        control_colony_count = "control_colony_count",
        viability = "viability",
        control = "control"
      )
    )
  )
}

build_figure_2_table <- function(paths) {
  all_replicates <- first_existing_result_table(
    paths,
    c(
      file.path("01_gfp_library_transfer_counts", "tables", "gfp_lib_standardized_all_replicates.csv"),
      file.path("01_gfp_library_transfer_counts", "tables", "standardized_all_replicates.csv")
    )
  )
  first_timepoint_summary <- read_result_table(paths, "01_gfp_library_transfer_counts", "tables", "first_timepoint_summary.csv")

  all_replicates <- all_replicates %>%
    dplyr::group_by(dataset_key, timepoint) %>%
    dplyr::mutate(
      control_colony_count = true_colony_count[well_position == "A09"][1],
      viability = true_colony_count / control_colony_count
    ) %>%
    dplyr::ungroup()

  first_timepoint_summary <- first_timepoint_summary %>%
    dplyr::group_by(dataset_key) %>%
    dplyr::mutate(
      control_colony_count = colony_count_avg[well_position == "A09"][1],
      viability = colony_count_avg / control_colony_count
    ) %>%
    dplyr::ungroup()

  dplyr::bind_rows(
    make_figure_table(
      all_replicates,
      figure = "Figure 2",
      panel = "GFP panel replicate kinetics",
      source_table = "results/01_gfp_library_transfer_counts/tables/gfp_lib_standardized_all_replicates.csv",
      mapping = list(
        assay = "dataset_key",
        cas_variant = "cas_variant",
        guide_promoter = "promoter",
        donor_recruitment = "donor_recruitment",
        replicate = "replicate",
        timepoint = "timepoint",
        generations = "total_gen_in_liquid",
        efficiency_pct = "editing_efficiency",
        colony_count = "true_colony_count",
        control_colony_count = "control_colony_count",
        viability = "viability",
        position = "position",
        well_position = "well_position",
        sample = "name"
      )
    ),
    make_figure_table(
      first_timepoint_summary,
      figure = "Figure 2",
      panel = "GFP panel first-timepoint viability",
      source_table = "results/01_gfp_library_transfer_counts/tables/first_timepoint_summary.csv",
      mapping = list(
        assay = "dataset_key",
        cas_variant = "cas_variant",
        guide_promoter = "promoter",
        donor_recruitment = "donor_recruitment",
        efficiency_pct = "efficiency_avg",
        colony_count = "colony_count_avg",
        control_colony_count = "control_colony_count",
        viability = "viability",
        well_position = "well_position"
      )
    )
  ) %>%
    dplyr::mutate(assay = dplyr::if_else(is.na(assay), "GFP_panel", assay))
}

build_figure_3_table <- function(paths) {
  paper_table <- try_build_figure_3_table_from_paper(paths)
  if (!is.null(paper_table)) {
    return(paper_table)
  }

  donor_positions <- read_result_table(paths, "03_donor_location", "tables", "donor_positions.csv")
  guide_length <- read_result_table(paths, "03_donor_location", "tables", "guide_length.csv")
  short_guide_viability <- read_result_table(paths, "03_donor_location", "tables", "short_guide_viability_plot_data.csv")
  donor_pos_fn <- read_result_table(paths, "03_donor_location", "tables", "donor_pos_fn.csv")
  donor_pos_enas <- read_result_table(paths, "03_donor_location", "tables", "donor_pos_enas.csv")

  donor_position_plotted <- donor_positions %>%
    dplyr::filter(
      (cas_variant == "enAsCas12a" & timepoint == "T1") |
        (cas_variant == "FnCas12a" & timepoint == "T2"),
      as.numeric(dist) != 6
    )

  guide_length_plotted <- guide_length %>%
    dplyr::filter(
      (cas_variant == "enAsCas12a" & timepoint == "T1") |
        (cas_variant == "FnCas12a" & timepoint == "T2")
    )

  donor_colonies_plotted <- dplyr::bind_rows(
    donor_pos_fn %>% dplyr::filter(timepoint == "T2"),
    donor_pos_enas %>% dplyr::filter(timepoint == "T1")
  ) %>%
    dplyr::group_by(cas_variant, timepoint) %>%
    dplyr::mutate(
      control_colony_count = max(true_colony_count, na.rm = TRUE),
      viability = true_colony_count / control_colony_count
    ) %>%
    dplyr::ungroup()

  short_guide_viability <- short_guide_viability %>%
    dplyr::mutate(efficiency_pct = 100 * as.numeric(editing_efficiency))

  donor_position_plotted <- donor_position_plotted %>%
    dplyr::mutate(non_hdr_pct = pmax(0, 100 - as.numeric(frc_alt) - as.numeric(frc_ref)))

  guide_length_plotted <- guide_length_plotted %>%
    dplyr::mutate(non_hdr_pct = pmax(0, 100 - as.numeric(frc_alt) - as.numeric(frc_ref)))

  dplyr::bind_rows(
    make_figure_table(
      donor_position_plotted,
      figure = "Figure 3",
      panel = "Donor position HDR",
      source_table = "results/03_donor_location/tables/donor_positions.csv",
      mapping = list(
        cas_variant = "cas_variant",
        strain = "strain",
        replicate = "replicate",
        timepoint = "timepoint",
        efficiency_pct = "hdr_pct",
        reference_pct = "frc_ref",
        non_hdr_pct = "non_hdr_pct",
        position = "dist"
      )
    ) %>% dplyr::mutate(assay = "donor_position", editing_class = "HDR"),
    make_figure_table(
      guide_length_plotted,
      figure = "Figure 3",
      panel = "Guide length HDR",
      source_table = "results/03_donor_location/tables/guide_length.csv",
      mapping = list(
        cas_variant = "cas_variant",
        strain = "strain",
        replicate = "replicate",
        timepoint = "timepoint",
        efficiency_pct = "hdr_pct",
        reference_pct = "frc_ref",
        non_hdr_pct = "non_hdr_pct",
        guide_length = "guide_length"
      )
    ) %>% dplyr::mutate(assay = "shorter_guide", editing_class = "HDR"),
    make_figure_table(
      short_guide_viability,
      figure = "Figure 3",
      panel = "Guide length viability",
      source_table = "results/03_donor_location/tables/short_guide_viability_plot_data.csv",
      mapping = list(
        cas_variant = "cas_variant",
        strain = "strain",
        replicate = "replicate",
        timepoint = "timepoint",
        generations = "total_gen_in_liquid",
        efficiency_pct = "efficiency_pct",
        viability = "viability",
        guide_length = "guide_length"
      )
    ) %>% dplyr::mutate(assay = "shorter_guide", editing_class = "HDR"),
    make_figure_table(
      donor_colonies_plotted,
      figure = "Figure 3",
      panel = "Donor position colony counts",
      source_table = "results/03_donor_location/tables/donor_pos_fn.csv; results/03_donor_location/tables/donor_pos_enas.csv",
      mapping = list(
        cas_variant = "cas_variant",
        strain = "strain",
        replicate = "replicate",
        timepoint = "timepoint",
        efficiency_pct = "hdr_pct",
        reference_pct = "frc_ref",
        colony_count = "true_colony_count",
        control_colony_count = "control_colony_count",
        viability = "viability",
        position = "dist"
      )
    ) %>% dplyr::mutate(assay = "donor_position", editing_class = "HDR")
  )
}

try_build_figure_3_table_from_paper <- function(paths) {
  root <- paths$paper_notebooks_root %||% ""
  required_files <- c(
    "processed_data/donor_position_processed_with_viability_20260421.csv",
    "processed_data/guide_length_processed_20260313.csv"
  )

  if (!nzchar(root) || !all(file.exists(file.path(root, required_files)))) {
    return(NULL)
  }

  read_paper_csv <- function(relative_path) {
    read_csv_clean(file.path(root, relative_path))
  }

  donor_position <- read_paper_csv("processed_data/donor_position_processed_with_viability_20260421.csv") %>%
    dplyr::filter(is_selected_timepoint) %>%
    dplyr::mutate(
      control_colony_count = true_colony_count / as.numeric(viability)
    )

  guide_length <- read_paper_csv("processed_data/guide_length_processed_20260313.csv") %>%
    dplyr::filter(is_selected_timepoint)

  dplyr::bind_rows(
    make_figure_table(
      donor_position,
      figure = "Figure 3",
      panel = "Donor position",
      source_table = "processed_data/donor_position_processed_with_viability_20260421.csv",
      mapping = list(
        cas_variant = "cas_variant",
        strain = "strain",
        replicate = "replicate",
        timepoint = "timepoint",
        efficiency_pct = "hdr_pct",
        reference_pct = "frc_ref",
        non_hdr_pct = "non_hdr_pct",
        colony_count = "true_colony_count",
        control_colony_count = "control_colony_count",
        viability = "viability",
        position = "dist"
      )
    ) %>% dplyr::mutate(assay = "donor_position", editing_class = "HDR"),
    make_figure_table(
      guide_length,
      figure = "Figure 3",
      panel = "Guide length",
      source_table = "processed_data/guide_length_processed_20260313.csv",
      mapping = list(
        cas_variant = "cas_variant",
        strain = "strain",
        replicate = "replicate",
        timepoint = "timepoint",
        efficiency_pct = "hdr_pct",
        reference_pct = "frc_ref",
        non_hdr_pct = "non_hdr_pct",
        viability = "viability",
        generations = "total_gen_in_liquid",
        guide_length = "guide_length",
        sample = "sample"
      )
    ) %>% dplyr::mutate(assay = "shorter_guide", editing_class = "HDR")
  )
}

build_figure_4_table <- function(paths) {
  paper_table <- try_build_figure_4_table_from_paper(paths)
  if (!is.null(paper_table)) {
    return(paper_table)
  }

  merged_data <- read_result_table(paths, "02_enas_amplicon_colony_counts", "tables", "merged_data.csv")

  hdr_vs_counts <- merged_data %>%
    dplyr::filter(timepoint == "T1", mutation == "HDR")

  composition_data <- merged_data %>%
    dplyr::filter(mutation == "HDR", timepoint == "T2") %>%
    dplyr::mutate(non_hdr = 100 - frc_alt - frc_ref) %>%
    dplyr::group_by(promoter, recruit, experiment, timepoint, guide) %>%
    dplyr::summarise(
      hdr = mean(frc_alt, na.rm = TRUE),
      non_hdr = mean(non_hdr, na.rm = TRUE),
      unedited = mean(frc_ref, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(c(hdr, non_hdr, unedited), names_to = "class", values_to = "pct")

  non_hdr_type_g10 <- merged_data %>%
    dplyr::filter(timepoint == "T2", guide == "G10", experiment == "E1") %>%
    dplyr::group_by(mutation) %>%
    dplyr::summarise(frc = mean(frc_alt, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(mutation %in% c("del", "ins"))

  hdr_kinetics <- merged_data %>%
    dplyr::filter(mutation == "HDR")

  dplyr::bind_rows(
    make_figure_table(
      hdr_vs_counts,
      figure = "Figure 4",
      panel = "Amplicon HDR versus colony count",
      source_table = "results/02_enas_amplicon_colony_counts/tables/merged_data.csv",
      mapping = list(
        assay = "aggregate",
        guide_id = "guide",
        guide_promoter = "promoter",
        donor_recruitment = "recruit",
        replicate = "replicate",
        timepoint = "timepoint",
        editing_class = "mutation",
        efficiency_pct = "frc_alt",
        reference_pct = "frc_ref",
        colony_count = "colonies_cutoff",
        condition_label = "promoter",
        experiment = "experiment",
        sample = "sample",
        variant = "variant"
      )
    ),
    make_figure_table(
      composition_data,
      figure = "Figure 4",
      panel = "Amplicon HDR/unedited/non-HDR composition",
      source_table = "results/02_enas_amplicon_colony_counts/tables/merged_data.csv",
      mapping = list(
        guide_id = "guide",
        guide_promoter = "promoter",
        donor_recruitment = "recruit",
        timepoint = "timepoint",
        editing_class = "class",
        efficiency_pct = "pct",
        condition_label = "promoter",
        experiment = "experiment"
      )
    ) %>% dplyr::mutate(assay = "amplicon_panel"),
    make_figure_table(
      non_hdr_type_g10,
      figure = "Figure 4",
      panel = "G10 non-HDR type",
      source_table = "results/02_enas_amplicon_colony_counts/tables/merged_data.csv",
      mapping = list(
        editing_class = "mutation",
        efficiency_pct = "frc"
      )
    ) %>% dplyr::mutate(assay = "amplicon_panel", guide_id = "G10", timepoint = "T2", experiment = "E1"),
    make_figure_table(
      hdr_kinetics,
      figure = "Figure 4",
      panel = "Amplicon HDR kinetics",
      source_table = "results/02_enas_amplicon_colony_counts/tables/merged_data.csv",
      mapping = list(
        assay = "aggregate",
        guide_id = "guide",
        guide_promoter = "promoter",
        donor_recruitment = "recruit",
        replicate = "replicate",
        timepoint = "timepoint",
        editing_class = "mutation",
        efficiency_pct = "frc_alt",
        reference_pct = "frc_ref",
        condition_label = "promoter",
        experiment = "experiment",
        sample = "sample",
        variant = "variant"
      )
    )
  )
}

try_build_figure_4_table_from_paper <- function(paths) {
  root <- paths$paper_notebooks_root %||% ""
  relative_path <- "processed_data/genomic_ampli_enas_fn_harmonized_20260318.csv"

  if (!nzchar(root) || !file.exists(file.path(root, relative_path))) {
    return(NULL)
  }

  amplicons <- read_csv_clean(file.path(root, relative_path)) %>%
    dplyr::mutate(
      colony_count = as.numeric(colony_count),
      total_gen_in_liquid = as.numeric(total_gen_in_liquid),
      non_hdr_pct = dplyr::if_else(
        mutation == "HDR",
        pmax(0, 100 - as.numeric(frc_alt) - as.numeric(frc_ref)),
        NA_real_
      )
    ) %>%
    dplyr::group_by(cas_variant, promoter, recruit) %>%
    dplyr::mutate(
      control_colony_count = if (all(is.na(colony_count))) NA_real_ else max(colony_count, na.rm = TRUE),
      viability = dplyr::if_else(
        is.na(colony_count) | is.na(control_colony_count) | control_colony_count == 0,
        NA_real_,
        colony_count / control_colony_count
      )
    ) %>%
    dplyr::ungroup()

  make_figure_table(
    amplicons,
    figure = "Figure 4",
    panel = "Amplicon long variant and colony data",
    source_table = relative_path,
    mapping = list(
      assay = "aggregate",
      cas_variant = "cas_variant",
      guide_promoter = "promoter",
      direct_repeat = "dr_cognate",
      donor_recruitment = "recruit",
      guide_id = "guide",
      strain = "ys",
      replicate = "replicate",
      timepoint = "timepoint",
      generations = "total_gen_in_liquid",
      editing_class = "mutation",
      efficiency_pct = "frc_alt",
      reference_pct = "frc_ref",
      non_hdr_pct = "non_hdr_pct",
      colony_count = "colony_count",
      control_colony_count = "control_colony_count",
      viability = "viability",
      experiment = "experiment"
    )
  )
}
