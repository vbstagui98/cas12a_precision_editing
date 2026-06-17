suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

check_flow_packages <- function() {
  packages <- c("flowCore", "flowWorkspace", "openCyto")
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop(
      paste(
        "Missing required Bioconductor packages for FCS processing:",
        paste(missing_packages, collapse = ", "),
        "\nInstall them with BiocManager before running the raw FCS step."
      ),
      call. = FALSE
    )
  }
}

read_csv_keep_names <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)
  semicolon_count <- stringr::str_count(first_line, ";")
  comma_count <- stringr::str_count(first_line, ",")

  if (semicolon_count > comma_count) {
    data <- readr::read_delim(path, delim = ";", show_col_types = FALSE, name_repair = "minimal")
  } else {
    data <- readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  }

  repair_blank_index_columns(data)
}

repair_blank_index_columns <- function(data) {
  names(data)[names(data) == ""] <- paste0("source_index_", seq_len(sum(names(data) == "")))
  data
}

normalise_well_name <- function(x) {
  x <- as.character(x)
  sub("([A-H])([1-9])$", "\\10\\2", x)
}

read_platemap_for_fcs <- function(path) {
  platemap <- read_csv_keep_names(path)

  if (!"name" %in% names(platemap)) {
    stop("The plate map must contain a `name` column with FCS file names.", call. = FALSE)
  }

  if ("WellPosition" %in% names(platemap)) {
    platemap$WellPosition <- normalise_well_name(platemap$WellPosition)
  }

  platemap$name <- as.character(platemap$name)
  platemap
}

read_optional_table <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(NULL)
  }
  read_csv_keep_names(path)
}

make_common_logicle_transform <- function(fs, channels) {
  data_1frame <- fs[[1]]
  flowCore::exprs(data_1frame) <- flowCore::fsApply(fs, function(x) {
    flowCore::exprs(x)
  })

  flowCore::estimateLogicle(data_1frame, channels = channels)
}

gate_gfp_population <- function(fstrans, gfp_channel, gate_range) {
  gate_gfp <- flowCore::fsApply(fstrans, function(x) {
    openCyto:::.mindensity2(
      x,
      channels = gfp_channel,
      filterId = "GFP",
      gate_range = gate_range
    )
  })

  gs <- flowWorkspace::GatingSet(fstrans)
  flowWorkspace::gs_pop_add(gs, gate_gfp, parent = "root")
  flowWorkspace::recompute(gs)

  as_tibble(flowWorkspace::gs_pop_get_count_fast(gs, statistic = "freq"))
}

read_gate_fcs <- function(
  fcs_dir,
  platemap_path,
  min_events = 500,
  transform_channels = c("FSC.A", "SSC.A", "FSC.H", "SSC.H", "BL1.A"),
  gfp_channel = "BL1.A",
  gate_range = c(1, 2.5)
) {
  check_flow_packages()

  platemap <- read_platemap_for_fcs(platemap_path)

  fs <- flowCore::read.flowSet(
    path = fcs_dir,
    pattern = ".fcs",
    alter.names = TRUE,
    emptyValue = FALSE
  )

  metadata <- dplyr::left_join(flowCore::pData(fs), platemap, by = "name")
  flowCore::pData(fs) <- metadata

  event_counts <- as.data.frame(flowCore::fsApply(fs, nrow))
  names(event_counts) <- "n_events"
  event_counts$name <- rownames(event_counts)

  keep_samples <- event_counts %>%
    dplyr::filter(n_events > min_events) %>%
    dplyr::pull(name)

  if (length(keep_samples) == 0) {
    stop("No FCS files passed the event-count filter.", call. = FALSE)
  }

  newfs <- fs[keep_samples]

  trans <- make_common_logicle_transform(newfs, transform_channels)
  fstrans <- flowCore::transform(newfs, trans)
  pop_stats <- gate_gfp_population(fstrans, gfp_channel, gate_range)

  list(
    pop_stats = pop_stats,
    event_counts = as_tibble(event_counts),
    platemap = platemap
  )
}

standardise_pop_stats <- function(pop_stats) {
  pop_stats <- repair_blank_index_columns(pop_stats)
  pop_stats <- as_tibble(pop_stats)

  required <- c("name", "Population", "Parent", "Frequency", "ParentFrequency")
  missing_required <- setdiff(required, names(pop_stats))
  if (length(missing_required) > 0) {
    stop(
      "Population stats are missing columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  pop_stats %>%
    dplyr::transmute(
      name = as.character(name),
      population = as.character(Population),
      parent = as.character(Parent),
      gfp_positive_fraction = as.numeric(Frequency),
      parent_frequency = as.numeric(ParentFrequency)
    )
}

join_optional_by_first_key <- function(data, extra, keys) {
  if (is.null(extra)) {
    return(data)
  }

  key <- keys[keys %in% names(data) & keys %in% names(extra)][1]
  if (is.na(key)) {
    warning("Could not join optional table; no shared key found.")
    return(data)
  }

  dplyr::left_join(data, extra, by = key)
}

editing_efficiency_from_gfp <- function(gfp_fraction, assay) {
  assay <- toupper(as.character(assay))

  dplyr::case_when(
    assay %in% c("GFP_ON", "ON", "RESTORE", "RESTORE_GFP") ~ gfp_fraction,
    assay %in% c("GFP_OFF", "OFF", "LOSS", "GFP_LOSS") ~ 1 - gfp_fraction,
    TRUE ~ 1 - gfp_fraction
  )
}

normalise_control_count <- function(colony_table, control_well = "A09") {
  if (!"true_colony_count" %in% names(colony_table)) {
    return(colony_table)
  }

  group_cols <- intersect(c("timepoint", "replicate", "dataset_key", "cas_variant"), names(colony_table))
  well_col <- intersect(c("WellPosition", "well_position", "well_attune"), names(colony_table))[1]

  if (is.na(well_col)) {
    return(
      colony_table %>%
        dplyr::mutate(
          control_colony_count = NA_real_,
          viability = NA_real_
        )
    )
  }

  controls <- colony_table %>%
    dplyr::filter(.data[[well_col]] == control_well) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(control_colony_count = mean(true_colony_count, na.rm = TRUE), .groups = "drop")

  colony_table %>%
    dplyr::left_join(controls, by = group_cols) %>%
    dplyr::mutate(
      viability = dplyr::if_else(
        !is.na(control_colony_count) & control_colony_count > 0,
        true_colony_count / control_colony_count,
        NA_real_
      )
    )
}

annotate_fcs_efficiency <- function(
  pop_stats,
  platemap_path,
  metadata_path = NULL,
  guide_features_path = NULL,
  colony_counts_path = NULL,
  assay_col = "assay",
  control_well = "A09"
) {
  platemap <- read_platemap_for_fcs(platemap_path)
  metadata <- read_optional_table(metadata_path)
  guide_features <- read_optional_table(guide_features_path)
  colony_counts <- read_optional_table(colony_counts_path)

  out <- standardise_pop_stats(pop_stats) %>%
    dplyr::left_join(platemap, by = "name")

  out <- join_optional_by_first_key(out, metadata, c("identifier", "SampleID", "WellPosition", "well_attune", "name"))
  out <- join_optional_by_first_key(out, guide_features, c("WellPosition", "SampleID", "guide", "guide_id"))
  out <- join_optional_by_first_key(out, colony_counts, c("name", "WellPosition", "well_attune", "identifier"))

  if (!assay_col %in% names(out)) {
    out[[assay_col]] <- "GFP_OFF"
  }

  if ("colony_count" %in% names(out) && !"true_colony_count" %in% names(out)) {
    out$true_colony_count <- suppressWarnings(as.numeric(out$colony_count))
  }

  out <- out %>%
    dplyr::mutate(
      editing_efficiency = editing_efficiency_from_gfp(gfp_positive_fraction, .data[[assay_col]])
    )

  normalise_control_count(out, control_well = control_well)
}

run_fcs_efficiency <- function(
  fcs_dir,
  platemap_path,
  output_dir,
  metadata_path = NULL,
  guide_features_path = NULL,
  colony_counts_path = NULL,
  min_events = 500,
  assay_col = "assay",
  control_well = "A09",
  transform_channels = c("FSC.A", "SSC.A", "FSC.H", "SSC.H", "BL1.A"),
  gfp_channel = "BL1.A",
  gate_range = c(1, 2.5)
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  gated <- read_gate_fcs(
    fcs_dir = fcs_dir,
    platemap_path = platemap_path,
    min_events = min_events,
    transform_channels = transform_channels,
    gfp_channel = gfp_channel,
    gate_range = gate_range
  )

  annotated <- annotate_fcs_efficiency(
    pop_stats = gated$pop_stats,
    platemap_path = platemap_path,
    metadata_path = metadata_path,
    guide_features_path = guide_features_path,
    colony_counts_path = colony_counts_path,
    assay_col = assay_col,
    control_well = control_well
  )

  readr::write_csv(gated$pop_stats, file.path(output_dir, "pop_stats.csv"))
  readr::write_csv(gated$event_counts, file.path(output_dir, "event_counts.csv"))
  readr::write_csv(annotated, file.path(output_dir, "fcs_efficiency.csv"))

  invisible(annotated)
}
