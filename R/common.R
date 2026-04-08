`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) {
    y
  } else {
    x
  }
}

load_paper_packages <- function(extra = character()) {
  packages <- unique(c(
    "dplyr",
    "ggplot2",
    "ggpmisc",
    "patchwork",
    "purrr",
    "readr",
    "readxl",
    "stringr",
    "tibble",
    "tidyr",
    "glue",
    extra
  ))

  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop(
      glue::glue(
        "Missing required R packages: {toString(missing_packages)}. ",
        "Run `Rscript scripts/install_packages.R` first."
      ),
      call. = FALSE
    )
  }

  invisible(lapply(
    packages,
    function(package_name) {
      suppressPackageStartupMessages(
        library(package_name, character.only = TRUE)
      )
    }
  ))
}

clean_column_names <- function(x) {
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
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

require_file <- function(path, label = NULL) {
  if (!file.exists(path)) {
    stop(glue::glue("Missing {(label %||% 'file')}: {path}"), call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

require_dir <- function(path, label = NULL) {
  if (!dir.exists(path)) {
    stop(glue::glue("Missing {(label %||% 'directory')}: {path}"), call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}

assert_columns <- function(data, required_columns, data_name = deparse(substitute(data))) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      glue::glue(
        "`{data_name}` is missing required columns: {toString(missing_columns)}"
      ),
      call. = FALSE
    )
  }
}

read_csv_clean <- function(path, ...) {
  data <- readr::read_csv(require_file(path), show_col_types = FALSE, ...)
  names(data) <- clean_column_names(names(data))
  tibble::as_tibble(data)
}

read_delim_clean <- function(path, delim, ...) {
  data <- readr::read_delim(
    require_file(path),
    delim = delim,
    show_col_types = FALSE,
    ...,
    name_repair = "minimal"
  )
  names(data) <- clean_column_names(names(data))
  tibble::as_tibble(data)
}

read_excel_clean <- function(path, sheet = 1, ...) {
  data <- readxl::read_excel(require_file(path), sheet = sheet, ...)
  names(data) <- clean_column_names(names(data))
  tibble::as_tibble(data)
}

parse_decimal_number <- function(x) {
  values <- trimws(as.character(x))
  values[values %in% c("", "NA", "NaN", "NULL")] <- NA_character_

  vapply(
    values,
    function(value) {
      if (is.na(value)) {
        return(NA_real_)
      }

      normalized <- gsub("\\s+", "", value)
      has_comma <- grepl(",", normalized, fixed = TRUE)
      has_dot <- grepl(".", normalized, fixed = TRUE)

      if (has_comma && has_dot) {
        last_comma <- max(gregexpr(",", normalized, fixed = TRUE)[[1]])
        last_dot <- max(gregexpr(".", normalized, fixed = TRUE)[[1]])

        if (last_comma > last_dot) {
          normalized <- gsub(".", "", normalized, fixed = TRUE)
          normalized <- sub(",", ".", normalized, fixed = TRUE)
        } else {
          normalized <- gsub(",", "", normalized, fixed = TRUE)
        }
      } else if (has_comma) {
        normalized <- sub(",", ".", normalized, fixed = TRUE)
      }

      suppressWarnings(as.numeric(normalized))
    },
    numeric(1)
  )
}

write_csv_checked <- function(data, path) {
  ensure_dir(dirname(path))
  readr::write_csv(data, path)
  invisible(path)
}

safe_ggsave <- function(plot, filename, width = 7, height = 5, dpi = 300, ...) {
  ensure_dir(dirname(filename))
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    ...
  )
  invisible(filename)
}

analysis_output_dir <- function(paths, notebook_slug) {
  ensure_dir(file.path(paths$results_dir, notebook_slug))
}

publication_theme <- function(base_size = 11, legend_position = "right", aspect_ratio = 1) {
  plot_theme <- ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_position,
      axis.ticks.length = grid::unit(2, "pt"),
      plot.margin = ggplot2::margin(5, 10, 5, 5),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
      strip.background = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.4, "lines")
    )

  if (!is.null(aspect_ratio)) {
    plot_theme <- plot_theme + ggplot2::theme(aspect.ratio = aspect_ratio)
  }

  plot_theme
}

variant_palette <- function() {
  c(
    "AsCas12a" = "#CC79A7",
    "enAsCas12a" = "#E69F00",
    "enAsCas12a_less_NLS" = "#F4A261",
    "FnCas12a" = "#009E73",
    "FnCas12a(CO)" = "#0072B2",
    "MbCas12a" = "#D55E00",
    "LbCas12a" = "#56B4E9",
    "LbCas12a(CO)" = "#2A9D8F"
  )
}

label_recruitment <- function(x) {
  dplyr::case_when(
    x %in% c("Yes", "LexA-FHA") ~ "LexA-FHA",
    TRUE ~ as.character(x)
  )
}

normalize_transfer_timepoint <- function(x) {
  dplyr::case_when(
    x %in% c("Timepoint_1", "timepoint_1", "T1") ~ "6_gen",
    x %in% c("Timepoint_2", "timepoint_2", "T2") ~ "11_gen",
    x %in% c("Timepoint_3", "timepoint_3", "T3") ~ "20_gen",
    x == "scraped" ~ "scraped",
    TRUE ~ as.character(x)
  )
}
