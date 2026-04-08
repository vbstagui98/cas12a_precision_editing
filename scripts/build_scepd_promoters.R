#!/usr/bin/env Rscript

script_file <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_file[grepl("^--file=", script_file)][1])
script_dir <- if (!is.na(script_file) && nzchar(script_file)) dirname(normalizePath(script_file, mustWork = FALSE)) else getwd()

source(file.path(script_dir, "..", "R", "common.R"))
source(file.path(script_dir, "..", "R", "yeast_variant_annotation.R"))

load_paper_packages(c("dplyr", "readr", "tibble"))

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/build_scepd_promoters.R",
      "--scepd raw_scEPD_file.tsv",
      "--output scepd_promoters.bed",
      "[--upstream 500]",
      "[--downstream 100]",
      sep = " "
    ),
    "\n\n",
    "Accepted scEPD input styles:\n",
    "- Headered table with recognizable chr/tss/strand/id columns.\n",
    "- Headerless SGA-like file: chr feature tss strand score id.\n",
    "- BED6 file with 1-bp promoter/TSS intervals.\n\n",
    "Output:\n",
    "- BED6 file where promoter intervals are derived from scEPDnew TSS entries.\n",
    "- On + strand: [TSS-upstream, TSS+downstream]\n",
    "- On - strand: [TSS-downstream, TSS+upstream]\n",
    sep = ""
  )
}

parse_args <- function(args) {
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = 0)
  }

  parsed <- list(
    upstream = 500,
    downstream = 100
  )

  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --name value pairs.", call. = FALSE)
  }

  for (i in seq(1, length(args), by = 2)) {
    key <- args[[i]]
    value <- args[[i + 1]]

    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }

    key <- gsub("-", "_", sub("^--", "", key))
    parsed[[key]] <- value
  }

  required <- c("scepd", "output")
  missing <- required[!required %in% names(parsed)]

  if (length(missing) > 0) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  parsed$upstream <- as.integer(parsed$upstream)
  parsed$downstream <- as.integer(parsed$downstream)

  if (any(is.na(c(parsed$upstream, parsed$downstream))) || parsed$upstream < 0 || parsed$downstream < 0) {
    stop("--upstream and --downstream must be non-negative integers.", call. = FALSE)
  }

  parsed
}

scepd_has_header <- function(path) {
  fields <- normalize_annotation_colname(first_line_fields(path))
  any(fields %in% c("chr", "chrom", "chromosome")) &&
    any(fields %in% c("tss", "position", "pos", "start", "chromstart"))
}

read_scepd_entries <- function(path) {
  require_file(path, "scEPD input")

  has_header <- scepd_has_header(path)
  df <- read_delim_base(path, header = has_header)

  if (has_header) {
    chr_idx <- pick_column(df, c("chr", "chrom", "chromosome"))
    tss_idx <- pick_column(df, c("tss", "position", "pos", "site"))
    strand_idx <- pick_column(df, c("strand"))
    id_idx <- pick_column(df, c("id", "name", "promoter", "promoterid", "gene", "geneid"))

    if (any(is.na(c(chr_idx, tss_idx, strand_idx)))) {
      bed_start_idx <- pick_column(df, c("start", "chromstart"))
      bed_end_idx <- pick_column(df, c("end", "chromend"))

      if (all(!is.na(c(chr_idx, bed_start_idx, bed_end_idx, strand_idx)))) {
        bed <- tibble::as_tibble(df) %>%
          dplyr::mutate(
            chr_raw = as.character(.data[[names(df)[[chr_idx]]]]),
            chr_norm = normalize_yeast_chr(chr_raw),
            bed_start = as.integer(.data[[names(df)[[bed_start_idx]]]]),
            bed_end = as.integer(.data[[names(df)[[bed_end_idx]]]]),
            strand = as.character(.data[[names(df)[[strand_idx]]]]),
            promoter_id = if (!is.na(id_idx)) as.character(.data[[names(df)[[id_idx]]]]) else NA_character_
          )

        widths <- bed$bed_end - bed$bed_start
        if (any(is.na(widths) | widths != 1L)) {
          stop("BED-like scEPD input must contain 1-bp intervals centered on the TSS.", call. = FALSE)
        }

        bed$tss <- ifelse(bed$strand == "-", bed$bed_end, bed$bed_start + 1L)
        return(bed[, c("chr_raw", "chr_norm", "tss", "strand", "promoter_id")])
      }

      stop("Could not identify chr/tss/strand columns in the scEPD input.", call. = FALSE)
    }

    entries <- tibble::as_tibble(df) %>%
      dplyr::mutate(
        chr_raw = as.character(.data[[names(df)[[chr_idx]]]]),
        chr_norm = normalize_yeast_chr(chr_raw),
        tss = as.integer(.data[[names(df)[[tss_idx]]]]),
        strand = as.character(.data[[names(df)[[strand_idx]]]]),
        promoter_id = if (!is.na(id_idx)) as.character(.data[[names(df)[[id_idx]]]]) else NA_character_
      ) %>%
      dplyr::select(chr_raw, chr_norm, tss, strand, promoter_id)

    return(entries)
  }

  if (ncol(df) >= 6 && !grepl("^[0-9]+$", as.character(df[[2]][1])) && grepl("^[0-9]+$", as.character(df[[3]][1])) && as.character(df[[4]][1]) %in% c("+", "-")) {
    entries <- tibble::tibble(
      chr_raw = as.character(df[[1]]),
      chr_norm = normalize_yeast_chr(df[[1]]),
      tss = as.integer(df[[3]]),
      strand = as.character(df[[4]]),
      promoter_id = as.character(df[[6]])
    )

    return(entries)
  }

  if (ncol(df) >= 6 && grepl("^[0-9]+$", as.character(df[[2]][1])) && grepl("^[0-9]+$", as.character(df[[3]][1])) && as.character(df[[6]][1]) %in% c("+", "-")) {
    bed_start <- as.integer(df[[2]])
    bed_end <- as.integer(df[[3]])
    widths <- bed_end - bed_start

    if (any(is.na(widths) | widths != 1L)) {
      stop("Headerless BED-like scEPD input must contain 1-bp intervals centered on the TSS.", call. = FALSE)
    }

    strand <- as.character(df[[6]])
    entries <- tibble::tibble(
      chr_raw = as.character(df[[1]]),
      chr_norm = normalize_yeast_chr(df[[1]]),
      tss = ifelse(strand == "-", bed_end, bed_start + 1L),
      strand = strand,
      promoter_id = as.character(df[[4]])
    )

    return(entries)
  }

  stop(
    "Unsupported scEPD input format. Use a headered chr/tss/strand table, SGA-like file, or BED6 1-bp TSS file.",
    call. = FALSE
  )
}

build_promoter_bed <- function(entries, upstream, downstream) {
  entries %>%
    dplyr::mutate(
      promoter_id = dplyr::if_else(
        is.na(promoter_id) | !nzchar(promoter_id),
        paste(chr_norm, tss, strand, sep = ":"),
        promoter_id
      ),
      strand = ifelse(strand %in% c("+", "-"), strand, "*"),
      start_1based = dplyr::if_else(
        strand == "-",
        pmax(1L, tss - downstream),
        pmax(1L, tss - upstream)
      ),
      end_1based = dplyr::if_else(
        strand == "-",
        tss + upstream,
        tss + downstream
      ),
      chrom = chr_norm,
      chrom_start = start_1based - 1L,
      chrom_end = end_1based,
      score = 0L
    ) %>%
    dplyr::select(chrom, chrom_start, chrom_end, promoter_id, score, strand)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  message("Reading scEPD input: ", args$scepd)
  entries <- read_scepd_entries(args$scepd)

  if (any(is.na(entries$chr_norm) | !nzchar(entries$chr_norm))) {
    stop("Some scEPD chromosomes could not be normalized to yeast chromosome names.", call. = FALSE)
  }

  if (any(is.na(entries$tss))) {
    stop("Some TSS positions could not be parsed as integers.", call. = FALSE)
  }

  if (any(!entries$strand %in% c("+", "-"))) {
    stop("scEPD entries must have strand values '+' or '-'.", call. = FALSE)
  }

  promoters <- build_promoter_bed(entries, args$upstream, args$downstream)

  ensure_dir(dirname(args$output))
  readr::write_tsv(promoters, args$output, col_names = FALSE)

  message("Wrote promoter BED: ", normalizePath(args$output, mustWork = FALSE))
  message("Promoter intervals: ", nrow(promoters))
  message("Window: upstream=", args$upstream, " downstream=", args$downstream)
}

main()
