#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

source("config/paths.R")
source("R/common.R")
source("R/amplicon_variant_efficiency.R")

load_paper_packages(extra = character())

paths <- get_analysis_paths()

publication_dir <- file.path(getwd(), "publication_source_data")
raw_variant_dir <- file.path(publication_dir, "raw_variant_tables")
ensure_dir(raw_variant_dir)

copy_if_present <- function(from, to) {
  require_file(from, basename(from))
  ensure_dir(dirname(to))
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) {
    stop(glue::glue("Could not copy {from} to {to}"), call. = FALSE)
  }
  normalizePath(to, mustWork = TRUE)
}

first_nonempty_file <- function(candidates, label) {
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    stop(
      glue::glue(
        "Could not find {label}. Checked:\n{paste(candidates, collapse = '\n')}"
      ),
      call. = FALSE
    )
  }

  sizes <- file.info(existing)$size
  nonempty <- existing[!is.na(sizes) & sizes > 0]
  if (length(nonempty) == 0) {
    stop(
      glue::glue(
        "All candidate files for {label} are empty. Checked:\n{paste(existing, collapse = '\n')}"
      ),
      call. = FALSE
    )
  }

  normalizePath(nonempty[[1]], mustWork = TRUE)
}

clean_for_export <- function(data) {
  names(data) <- clean_column_names(names(data))
  numeric_names <- grepl("^[0-9]+$", names(data))
  names(data)[numeric_names] <- paste0("unnamed_", names(data)[numeric_names])

  data %>%
    tibble::as_tibble() %>%
    dplyr::relocate(
      dplyr::any_of(c(
        "raw_variant_table",
        "source_dataset",
        "source_file",
        "cas_variant",
        "assay",
        "sample"
      ))
    )
}

add_simplified_variant_columns <- function(data) {
  required <- c("POS", "REF", "ALT", "TYPE")
  if (!all(required %in% names(data))) {
    return(data)
  }

  simplify_variant_table(data)
}

normalize_sample_annotation_columns <- function(data) {
  data %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::matches("^Part[0-9]+$|^\\.\\.\\.[0-9]+$"),
        as.character
      )
    )
}

amplicon_dataset_name <- function(path) {
  basename(dirname(dirname(dirname(path))))
}

export_variant_table <- function(data, path) {
  write_csv_checked(clean_for_export(data), path)
}

figure_workbook <- file.path(
  paths$results_dir,
  "publication_figure_tables",
  "publication_figure_source_tables_clean_20260617.xlsx"
)
copied_workbook <- copy_if_present(
  figure_workbook,
  file.path(publication_dir, "publication_figure_source_tables_clean_20260617.xlsx")
)

presentation_amplicon_root <- Sys.getenv(
  "CAS12A_PRESENTATIONS_AMPLICON_ROOT",
  unset = file.path(
    dirname(dirname(paths$paper_notebooks_root)),
    "Mac (2)",
    "Documents",
    "PhD",
    "Presentations",
    "9-month_figs",
    "amplicon_seq"
  )
)

fn_donor_path <- first_nonempty_file(
  c(
    file.path(paths$presentations_amplicon_root, "20250312_ys88_donor_mut", "data", "raw", "master_df_filtered.csv"),
    file.path(presentation_amplicon_root, "20250312_ys88_donor_mut", "data", "raw", "master_df_filtered.csv")
  ),
  "FnCas12a donor-position raw variant table"
)

enas_donor_path <- first_nonempty_file(
  c(
    file.path(paths$presentations_amplicon_root, "20250611_amplicons_ys85_misc", "data", "raw", "master_df_filtered.csv"),
    file.path(presentation_amplicon_root, "20250611_amplicons_ys85_misc", "data", "raw", "master_df_filtered.csv")
  ),
  "enAsCas12a donor-position raw variant table"
)

short_guide_path <- first_nonempty_file(
  c(
    file.path(paths$presentations_amplicon_root, "20250611_amplicons_shorterguide", "data", "raw", "master_df_filtered.csv"),
    file.path(presentation_amplicon_root, "20250611_amplicons_shorterguide", "data", "raw", "master_df_filtered.csv")
  ),
  "short-guide raw variant table"
)

genomic_panel_path <- first_nonempty_file(
  c(
    file.path(paths$paper_notebooks_root, "processed_data", "genomic_ampli_enas_fn_harmonized_20260318.csv"),
    file.path(paths$results_dir, "publication_figure_tables", "figure_4_source_data.csv")
  ),
  "genomic amplicon panel long variant table"
)

fn_donor_variants <- read_variant_table(fn_donor_path) %>%
  prepare_fn_donor_variants() %>%
  add_simplified_variant_columns() %>%
  normalize_sample_annotation_columns() %>%
  dplyr::mutate(
    raw_variant_table = "donor_position",
    source_dataset = amplicon_dataset_name(fn_donor_path),
    cas_variant = "FnCas12a",
    assay = "different_donor",
    source_file = basename(fn_donor_path),
    .before = 1
  )

enas_donor_variants <- read_variant_table(enas_donor_path) %>%
  prepare_enas_donor_variants() %>%
  add_simplified_variant_columns() %>%
  normalize_sample_annotation_columns() %>%
  dplyr::mutate(
    raw_variant_table = "donor_position",
    source_dataset = amplicon_dataset_name(enas_donor_path),
    cas_variant = "enAsCas12a",
    assay = "different_donor",
    source_file = basename(enas_donor_path),
    .before = 1
  )

short_guide_variants <- read_variant_table(short_guide_path) %>%
  prepare_short_guide_variants() %>%
  add_simplified_variant_columns() %>%
  normalize_sample_annotation_columns() %>%
  dplyr::mutate(
    raw_variant_table = "short_guide",
    source_dataset = amplicon_dataset_name(short_guide_path),
    assay = "short_guide",
    cas_variant = dplyr::case_when(
      strain == "ys85" ~ "enAsCas12a",
      strain == "ys88" ~ "FnCas12a",
      TRUE ~ as.character(strain)
    ),
    source_file = basename(short_guide_path),
    .before = 1
  )

genomic_panel_variants <- read_variant_table(genomic_panel_path) %>%
  normalize_sample_annotation_columns() %>%
  dplyr::mutate(
    raw_variant_table = "genomic_amplicon_panel",
    source_dataset = "cas12a_paper_notebooks_20260312/processed_data",
    assay = "genomic_amplicon_panel",
    source_file = basename(genomic_panel_path),
    .before = 1
  )

export_paths <- c(
  donor_position = file.path(raw_variant_dir, "donor_position_all_detected_variants.csv"),
  short_guide = file.path(raw_variant_dir, "short_guide_all_detected_variants.csv"),
  genomic_panel = file.path(raw_variant_dir, "genomic_panel_all_detected_variants.csv")
)

donor_position_variants <- dplyr::bind_rows(fn_donor_variants, enas_donor_variants)

export_variant_table(
  donor_position_variants,
  export_paths[["donor_position"]]
)
export_variant_table(short_guide_variants, export_paths[["short_guide"]])
export_variant_table(genomic_panel_variants, export_paths[["genomic_panel"]])

manifest <- tibble::tribble(
  ~file, ~rows, ~columns, ~source,
  "publication_figure_source_tables_clean_20260617.xlsx", NA_integer_, NA_integer_, "results/publication_figure_tables/publication_figure_source_tables_clean_20260617.xlsx",
  "raw_variant_tables/donor_position_all_detected_variants.csv", nrow(donor_position_variants), ncol(clean_for_export(donor_position_variants)), paste(c(
    "20250312_ys88_donor_mut/data/raw/master_df_filtered.csv",
    "20250611_amplicons_ys85_misc/data/raw/master_df_filtered.csv"
  ), collapse = "; "),
  "raw_variant_tables/short_guide_all_detected_variants.csv", nrow(short_guide_variants), ncol(clean_for_export(short_guide_variants)), "20250611_amplicons_shorterguide/data/raw/master_df_filtered.csv",
  "raw_variant_tables/genomic_panel_all_detected_variants.csv", nrow(genomic_panel_variants), ncol(clean_for_export(genomic_panel_variants)), "cas12a_paper_notebooks_20260312/processed_data/genomic_ampli_enas_fn_harmonized_20260318.csv"
)

readr::write_csv(manifest, file.path(publication_dir, "manifest.csv"))

readr::write_lines(
  c(
    "# Publication source data",
    "",
    "This folder contains the figure-source workbook and clean long variant tables committed for the Cas12a precision-editing publication repository.",
    "",
    "## Files",
    "",
    "- `publication_figure_source_tables_clean_20260617.xlsx`: one worksheet per final figure source-data table.",
    "- `raw_variant_tables/donor_position_all_detected_variants.csv`: all detected variants from the donor-position amplicon sequencing runs after applying the same donor sample annotations used for the efficiency summaries.",
    "- `raw_variant_tables/short_guide_all_detected_variants.csv`: all detected variants from the short-guide amplicon sequencing run after applying the same guide-length annotations used for the efficiency summaries.",
    "- `raw_variant_tables/genomic_panel_all_detected_variants.csv`: all detected variants from `genomic_ampli_enas_fn_harmonized_20260318.csv`, the long table used for the genomic amplicon panel figure.",
    "- `manifest.csv`: row counts and source paths used for the export.",
    "",
    "Regenerate these files from the project root with:",
    "",
    "```bash",
    "CAS12A_PAPER_NOTEBOOKS_ROOT=\"/path/to/cas12a_paper_notebooks_20260312\" \\",
    "CAS12A_PRESENTATIONS_AMPLICON_ROOT=\"/path/to/9-month_figs/amplicon_seq\" \\",
    "Rscript scripts/export_publication_source_data.R",
    "```"
  ),
  file.path(publication_dir, "README.md")
)

cat("Wrote publication source data to:", normalizePath(publication_dir), "\n")
print(manifest)
