get_script_path <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_arg <- file_arg[grepl("^--file=", file_arg)]

  if (length(file_arg) == 0) {
    return(normalizePath(getwd(), mustWork = FALSE))
  }

  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
}

get_project_root <- function() {
  script_path <- get_script_path()

  if (dir.exists(script_path)) {
    return(script_path)
  }

  normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
}

source_project_file <- function(relative_path) {
  source(file.path(get_project_root(), relative_path), local = FALSE)
}

normalize_annotation_colname <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(x))
}

guess_delim <- function(path) {
  lines <- readLines(path, n = 5, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  if (length(lines) == 0) {
    stop("Input file is empty: ", path, call. = FALSE)
  }

  first <- lines[[1]]

  if (grepl("\t", first, fixed = TRUE)) {
    return("\t")
  }

  if (grepl(",", first, fixed = TRUE)) {
    return(",")
  }

  ""
}

read_delim_base <- function(path, header = TRUE, comment_char = "", sep = NULL) {
  actual_sep <- sep %||% guess_delim(path)

  if (actual_sep == "\t") {
    return(
      utils::read.delim(
        file = path,
        header = header,
        stringsAsFactors = FALSE,
        comment.char = comment_char,
        quote = "",
        check.names = FALSE
      )
    )
  }

  if (actual_sep == ",") {
    return(
      utils::read.csv(
        file = path,
        header = header,
        stringsAsFactors = FALSE,
        comment.char = comment_char,
        quote = "",
        check.names = FALSE
      )
    )
  }

  utils::read.table(
    file = path,
    sep = actual_sep,
    header = header,
    stringsAsFactors = FALSE,
    comment.char = comment_char,
    quote = "",
    check.names = FALSE
  )
}

first_line_fields <- function(path) {
  sep <- guess_delim(path)
  first_line <- readLines(path, n = 1, warn = FALSE)

  if (!nzchar(first_line)) {
    stop("Input file is empty: ", path, call. = FALSE)
  }

  fields <- strsplit(first_line, if (sep == "") "[[:space:]]+" else sep)[[1]]
  fields[nzchar(fields)]
}

pick_column <- function(df, candidates) {
  normalized <- normalize_annotation_colname(names(df))
  match_idx <- match(candidates, normalized)
  match_idx <- match_idx[!is.na(match_idx)]

  if (length(match_idx) == 0) {
    return(NA_integer_)
  }

  match_idx[[1]]
}

normalize_yeast_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_

  normalized <- toupper(x)
  normalized <- sub("^REF\\|", "", normalized)
  normalized <- sub("\\|$", "", normalized)
  normalized <- sub("^CHROMOSOME", "", normalized)
  normalized <- sub("^CHR", "", normalized)
  normalized <- sub("\\.\\d+$", "", normalized)

  refseq_map <- c(
    "NC_001133" = "I",
    "NC_001134" = "II",
    "NC_001135" = "III",
    "NC_001136" = "IV",
    "NC_001137" = "V",
    "NC_001138" = "VI",
    "NC_001139" = "VII",
    "NC_001140" = "VIII",
    "NC_001141" = "IX",
    "NC_001142" = "X",
    "NC_001143" = "XI",
    "NC_001144" = "XII",
    "NC_001145" = "XIII",
    "NC_001146" = "XIV",
    "NC_001147" = "XV",
    "NC_001148" = "XVI",
    "NC_001224" = "M"
  )

  arabic_map <- c(
    "1" = "I",
    "2" = "II",
    "3" = "III",
    "4" = "IV",
    "5" = "V",
    "6" = "VI",
    "7" = "VII",
    "8" = "VIII",
    "9" = "IX",
    "10" = "X",
    "11" = "XI",
    "12" = "XII",
    "13" = "XIII",
    "14" = "XIV",
    "15" = "XV",
    "16" = "XVI"
  )

  mito_aliases <- c("M", "MT", "MITO", "MITOCHONDRIA", "MITOCHONDRIAL")

  normalized[normalized %in% names(refseq_map)] <- refseq_map[normalized[normalized %in% names(refseq_map)]]
  normalized[normalized %in% names(arabic_map)] <- arabic_map[normalized[normalized %in% names(arabic_map)]]
  normalized[normalized %in% mito_aliases] <- "M"

  roman_chrs <- c(
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII",
    "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"
  )

  normalized[normalized %in% roman_chrs] <- paste0("chr", normalized[normalized %in% roman_chrs])
  normalized
}

extract_fasta_token <- function(header) {
  token <- strsplit(sub("^>", "", header), "[[:space:]]+")[[1]][1]

  refseq_hit <- regmatches(token, regexpr("NC_[0-9]+(\\.[0-9]+)?", token))
  if (length(refseq_hit) == 1 && nzchar(refseq_hit)) {
    return(refseq_hit)
  }

  token
}

read_fasta_sequences <- function(path) {
  require_file(path, "genome FASTA")

  lines <- readLines(path, warn = FALSE)
  header_idx <- which(startsWith(lines, ">"))

  if (length(header_idx) == 0) {
    stop("FASTA contains no header lines: ", path, call. = FALSE)
  }

  sequences <- list()
  raw_names <- character()

  for (i in seq_along(header_idx)) {
    start_line <- header_idx[[i]]
    end_line <- if (i < length(header_idx)) header_idx[[i + 1]] - 1 else length(lines)
    header <- lines[[start_line]]
    sequence <- paste(lines[(start_line + 1):end_line], collapse = "")
    raw_name <- extract_fasta_token(header)
    chr_norm <- normalize_yeast_chr(raw_name)

    if (is.na(chr_norm) || !nzchar(chr_norm)) {
      next
    }

    if (!is.null(sequences[[chr_norm]])) {
      stop("FASTA contains duplicate chromosome after normalization: ", chr_norm, call. = FALSE)
    }

    sequences[[chr_norm]] <- toupper(sequence)
    raw_names[[chr_norm]] <- raw_name
  }

  list(
    sequences = sequences,
    raw_names = raw_names,
    lengths = vapply(sequences, nchar, integer(1))
  )
}

parse_gff_attributes <- function(x) {
  fields <- strsplit(as.character(x), ";", fixed = TRUE)[[1]]
  fields <- trimws(fields)
  fields <- fields[nzchar(fields)]

  out <- list()

  for (field in fields) {
    if (!grepl("=", field, fixed = TRUE)) {
      next
    }

    key <- sub("=.*$", "", field)
    value <- sub("^[^=]*=", "", field)
    out[[key]] <- utils::URLdecode(value)
  }

  out
}

pick_gff_attribute <- function(attributes, keys) {
  for (key in keys) {
    value <- attributes[[key]]
    if (!is.null(value) && nzchar(value)) {
      return(value)
    }
  }

  NA_character_
}

first_non_missing <- function(...) {
  values <- c(...)
  values <- values[!is.na(values) & nzchar(values)]

  if (length(values) == 0) {
    return(NA_character_)
  }

  values[[1]]
}

parse_info_value <- function(info_string, key) {
  pieces <- strsplit(as.character(info_string), ";", fixed = TRUE)[[1]]
  target <- paste0(key, "=")
  hit <- pieces[startsWith(pieces, target)]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  sub(target, "", hit[[1]], fixed = TRUE)
}

parse_snpeff_ann <- function(ann_string) {
  if (is.na(ann_string) || !nzchar(ann_string)) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  entries <- strsplit(ann_string, ",", fixed = TRUE)[[1]]

  parsed <- lapply(
    entries,
    function(entry) {
      fields <- strsplit(entry, "|", fixed = TRUE)[[1]]
      length(fields) <- 16

      data.frame(
        allele = fields[[1]] %||% NA_character_,
        effect = fields[[2]] %||% NA_character_,
        impact = fields[[3]] %||% NA_character_,
        gene_name = fields[[4]] %||% NA_character_,
        gene_id = fields[[5]] %||% NA_character_,
        feature_type = fields[[6]] %||% NA_character_,
        feature_id = fields[[7]] %||% NA_character_,
        transcript_biotype = fields[[8]] %||% NA_character_,
        rank = fields[[9]] %||% NA_character_,
        hgvs_c = fields[[10]] %||% NA_character_,
        hgvs_p = fields[[11]] %||% NA_character_,
        cdna_pos = fields[[12]] %||% NA_character_,
        cds_pos = fields[[13]] %||% NA_character_,
        aa_pos = fields[[14]] %||% NA_character_,
        distance = fields[[15]] %||% NA_character_,
        warnings = fields[[16]] %||% NA_character_,
        stringsAsFactors = FALSE
      )
    }
  )

  dplyr::bind_rows(parsed)
}

collapse_unique_values <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = ";")
}
