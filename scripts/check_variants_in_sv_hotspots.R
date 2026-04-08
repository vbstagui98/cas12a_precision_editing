#!/usr/bin/env Rscript

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/check_variants_in_sv_hotspots.R",
      "--hotspots SCORE_SV_hotspots.bed",
      "--scores s288c_NGG_PAM_sites_editing_SCORE.tsv",
      "--variants variants.tsv",
      "--output annotated_variants.tsv",
      "[--score-column SCORE.SV.minmax]",
      "[--score-threshold 0.9]",
      sep = " "
    ),
    "\n\n",
    "Variant file format:\n",
    "- Headered file with columns like chr,pos,ref,alt (case-insensitive), or\n",
    "- Headerless file where the first 4 columns are chr pos ref alt.\n\n",
    "Output:\n",
    "- One row per variant-hotspot overlap.\n",
    "- Variants with no overlap still get one row with hotspot fields set to NA.\n",
    sep = ""
  )
}

parse_args <- function(args) {
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = 0)
  }

  defaults <- list(
    score_column = "SCORE.SV.minmax",
    score_threshold = 0.9
  )

  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --name value pairs.", call. = FALSE)
  }

  parsed <- defaults

  for (i in seq(1, length(args), by = 2)) {
    key <- args[[i]]
    value <- args[[i + 1]]

    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }

    key <- gsub("-", "_", sub("^--", "", key))
    parsed[[key]] <- value
  }

  required <- c("hotspots", "scores", "variants", "output")
  missing <- required[!required %in% names(parsed)]

  if (length(missing) > 0) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  parsed$score_threshold <- as.numeric(parsed$score_threshold)

  if (is.na(parsed$score_threshold)) {
    stop("--score-threshold must be numeric.", call. = FALSE)
  }

  parsed
}

guess_sep <- function(path) {
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

normalize_colname <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(x))
}

normalize_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  x <- sub("^CHROMOSOME", "", x)
  x <- sub("^CHR", "", x)

  arabic_to_roman <- c(
    "1" = "I", "2" = "II", "3" = "III", "4" = "IV",
    "5" = "V", "6" = "VI", "7" = "VII", "8" = "VIII",
    "9" = "IX", "10" = "X", "11" = "XI", "12" = "XII",
    "13" = "XIII", "14" = "XIV", "15" = "XV", "16" = "XVI"
  )

  mito_aliases <- c("M", "MT", "MITO", "MITOCHONDRIA", "MITOCHONDRIAL")
  x[x %in% names(arabic_to_roman)] <- arabic_to_roman[x[x %in% names(arabic_to_roman)]]
  x[x %in% mito_aliases] <- "M"
  x
}

read_delim_base <- function(path, header = TRUE) {
  sep <- guess_sep(path)

  if (sep == "\t") {
    return(
      utils::read.delim(
        file = path,
        header = header,
        stringsAsFactors = FALSE,
        comment.char = "",
        quote = "",
        check.names = FALSE
      )
    )
  }

  if (sep == ",") {
    return(
      utils::read.csv(
        file = path,
        header = header,
        stringsAsFactors = FALSE,
        comment.char = "",
        quote = "",
        check.names = FALSE
      )
    )
  }

  utils::read.table(
    file = path,
    sep = sep,
    header = header,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = "",
    check.names = FALSE
  )
}

variant_file_has_header <- function(path) {
  sep <- guess_sep(path)
  first_line <- readLines(path, n = 1, warn = FALSE)

  if (!nzchar(first_line)) {
    stop("Variant file is empty.", call. = FALSE)
  }

  fields <- strsplit(first_line, if (sep == "") "[[:space:]]+" else sep)[[1]]
  fields <- fields[nzchar(fields)]
  normalized <- normalize_colname(fields)

  expected <- c("chr", "chrom", "chromosome")
  any(normalized %in% expected)
}

pick_column <- function(df, candidates) {
  normalized <- normalize_colname(names(df))
  match_idx <- match(candidates, normalized)
  match_idx <- match_idx[!is.na(match_idx)]

  if (length(match_idx) == 0) {
    return(NA_integer_)
  }

  match_idx[[1]]
}

read_variants <- function(path) {
  has_header <- variant_file_has_header(path)
  df <- read_delim_base(path, header = has_header)

  if (!has_header) {
    if (ncol(df) < 4) {
      stop("Headerless variant file must have at least 4 columns: chr pos ref alt.", call. = FALSE)
    }

    names(df)[1:4] <- c("chr", "pos", "ref", "alt")
  }

  chr_idx <- pick_column(df, c("chr", "chrom", "chromosome"))
  pos_idx <- pick_column(df, c("pos", "position", "site"))
  ref_idx <- pick_column(df, c("ref", "reference"))
  alt_idx <- pick_column(df, c("alt", "alternate", "alternative"))

  if (any(is.na(c(chr_idx, pos_idx, ref_idx, alt_idx)))) {
    stop("Could not identify chr/pos/ref/alt columns in the variant file.", call. = FALSE)
  }

  df$variant_index <- seq_len(nrow(df))
  df$variant_chr <- df[[chr_idx]]
  df$variant_pos <- as.integer(df[[pos_idx]])
  df$variant_ref <- df[[ref_idx]]
  df$variant_alt <- df[[alt_idx]]
  df$variant_chr_norm <- normalize_chr(df$variant_chr)
  df$variant_id <- paste(df$variant_chr, df$variant_pos, df$variant_ref, df$variant_alt, sep = ":")

  if (any(is.na(df$variant_pos))) {
    stop("Some variant positions could not be parsed as integers.", call. = FALSE)
  }

  df
}

read_hotspots <- function(path) {
  bed <- read_delim_base(path, header = FALSE)

  if (ncol(bed) < 3) {
    stop("Hotspot BED file must have at least 3 columns: chr start end.", call. = FALSE)
  }

  hotspots <- data.frame(
    hotspot_index = seq_len(nrow(bed)),
    hotspot_chr_raw = bed[[1]],
    hotspot_start_0based = as.integer(bed[[2]]),
    hotspot_end_1based = as.integer(bed[[3]]),
    stringsAsFactors = FALSE
  )

  if (any(is.na(hotspots$hotspot_start_0based) | is.na(hotspots$hotspot_end_1based))) {
    stop("Hotspot BED start/end columns must be integers.", call. = FALSE)
  }

  hotspots$hotspot_chr_norm <- normalize_chr(hotspots$hotspot_chr_raw)
  hotspots$hotspot_start_1based <- hotspots$hotspot_start_0based + 1L
  hotspots$hotspot_width_bp <- hotspots$hotspot_end_1based - hotspots$hotspot_start_1based + 1L
  hotspots
}

read_scores <- function(path, score_column) {
  sep <- guess_sep(path)
  header_line <- readLines(path, n = 1, warn = FALSE)
  header_fields <- strsplit(header_line, if (sep == "") "[[:space:]]+" else sep)[[1]]
  header_fields <- header_fields[nzchar(header_fields)]
  required <- c("chr", "pos", score_column)
  optional <- c("strand", "guide")
  requested <- unique(c(required, optional))
  missing <- required[!required %in% header_fields]

  if (length(missing) > 0) {
    stop(
      "Score file is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  col_classes <- rep("NULL", length(header_fields))
  names(col_classes) <- header_fields
  col_classes["chr"] <- "character"
  col_classes["pos"] <- "numeric"
  col_classes[score_column] <- "numeric"

  if ("strand" %in% header_fields) {
    col_classes["strand"] <- "character"
  }

  if ("guide" %in% header_fields) {
    col_classes["guide"] <- "character"
  }

  if (sep == "\t") {
    scores <- utils::read.delim(
      file = path,
      header = TRUE,
      stringsAsFactors = FALSE,
      comment.char = "",
      quote = "",
      check.names = FALSE,
      colClasses = unname(col_classes)
    )
  } else if (sep == ",") {
    scores <- utils::read.csv(
      file = path,
      header = TRUE,
      stringsAsFactors = FALSE,
      comment.char = "",
      quote = "",
      check.names = FALSE,
      colClasses = unname(col_classes)
    )
  } else {
    scores <- utils::read.table(
      file = path,
      sep = sep,
      header = TRUE,
      stringsAsFactors = FALSE,
      comment.char = "",
      quote = "",
      check.names = FALSE,
      colClasses = unname(col_classes)
    )
  }

  scores$chr_norm <- normalize_chr(scores$chr)
  scores$pos <- as.integer(scores$pos)
  scores$score_sv <- as.numeric(scores[[score_column]])
  scores$score_key <- paste(scores$chr_norm, scores$pos, sep = "::")

  if (any(is.na(scores$pos))) {
    stop("Some score positions could not be parsed as integers.", call. = FALSE)
  }

  if (any(is.na(scores$score_sv))) {
    stop("Some values in the selected score column are not numeric.", call. = FALSE)
  }

  scores
}

summarise_hotspot_scores <- function(hotspots, scores, score_threshold) {
  scores_by_chr <- split(scores, scores$chr_norm)
  summaries <- vector("list", nrow(hotspots))

  for (i in seq_len(nrow(hotspots))) {
    chr_scores <- scores_by_chr[[hotspots$hotspot_chr_norm[[i]]]]

    if (is.null(chr_scores) || nrow(chr_scores) == 0) {
      summaries[[i]] <- data.frame(
        hotspot_score_site_count = 0L,
        hotspot_mean_sv_score = NA_real_,
        hotspot_max_sv_score = NA_real_,
        hotspot_high_score = NA,
        hotspot_best_score_pos = NA_integer_,
        hotspot_best_score_strand = NA_character_,
        hotspot_best_score_guide = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }

    in_interval <- chr_scores$pos >= hotspots$hotspot_start_1based[[i]] &
      chr_scores$pos <= hotspots$hotspot_end_1based[[i]]
    hits <- chr_scores[in_interval, , drop = FALSE]

    if (nrow(hits) == 0) {
      summaries[[i]] <- data.frame(
        hotspot_score_site_count = 0L,
        hotspot_mean_sv_score = NA_real_,
        hotspot_max_sv_score = NA_real_,
        hotspot_high_score = FALSE,
        hotspot_best_score_pos = NA_integer_,
        hotspot_best_score_strand = NA_character_,
        hotspot_best_score_guide = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }

    best_idx <- which.max(hits$score_sv)
    best_row <- hits[best_idx, , drop = FALSE]

    summaries[[i]] <- data.frame(
      hotspot_score_site_count = nrow(hits),
      hotspot_mean_sv_score = mean(hits$score_sv, na.rm = TRUE),
      hotspot_max_sv_score = max(hits$score_sv, na.rm = TRUE),
      hotspot_high_score = max(hits$score_sv, na.rm = TRUE) >= score_threshold,
      hotspot_best_score_pos = best_row$pos[[1]],
      hotspot_best_score_strand = if ("strand" %in% names(best_row)) best_row$strand[[1]] else NA_character_,
      hotspot_best_score_guide = if ("guide" %in% names(best_row)) best_row$guide[[1]] else NA_character_,
      stringsAsFactors = FALSE
    )
  }

  cbind(hotspots, do.call(rbind, summaries))
}

summarise_exact_position_scores <- function(variants, scores, score_threshold) {
  best_order <- order(scores$score_key, -scores$score_sv)
  best_rows <- scores[best_order, , drop = FALSE]
  best_rows <- best_rows[!duplicated(best_rows$score_key), , drop = FALSE]

  counts <- as.data.frame(table(scores$score_key), stringsAsFactors = FALSE)
  names(counts) <- c("score_key", "exact_score_site_count")

  exact_summary <- merge(
    counts,
    best_rows[, c("score_key", "score_sv", intersect(c("strand", "guide"), names(best_rows))), drop = FALSE],
    by = "score_key",
    all.x = TRUE,
    sort = FALSE
  )

  variant_keys <- paste(variants$variant_chr_norm, variants$variant_pos, sep = "::")
  match_idx <- match(variant_keys, exact_summary$score_key)

  out <- data.frame(
    exact_score_site_count = ifelse(is.na(match_idx), 0L, as.integer(exact_summary$exact_score_site_count[match_idx])),
    exact_max_sv_score = ifelse(is.na(match_idx), NA_real_, exact_summary$score_sv[match_idx]),
    exact_high_score = ifelse(
      is.na(match_idx),
      FALSE,
      exact_summary$score_sv[match_idx] >= score_threshold
    ),
    exact_best_score_strand = NA_character_,
    exact_best_score_guide = NA_character_,
    stringsAsFactors = FALSE
  )

  if ("strand" %in% names(exact_summary)) {
    out$exact_best_score_strand <- ifelse(
      is.na(match_idx),
      NA_character_,
      as.character(exact_summary$strand[match_idx])
    )
  }

  if ("guide" %in% names(exact_summary)) {
    out$exact_best_score_guide <- ifelse(
      is.na(match_idx),
      NA_character_,
      as.character(exact_summary$guide[match_idx])
    )
  }

  out
}

annotate_variants <- function(variants, hotspots_scored, exact_scores) {
  hotspots_by_chr <- split(hotspots_scored, hotspots_scored$hotspot_chr_norm)
  annotated_rows <- list()
  row_counter <- 1L

  for (i in seq_len(nrow(variants))) {
    chr_hotspots <- hotspots_by_chr[[variants$variant_chr_norm[[i]]]]
    overlaps <- NULL

    if (!is.null(chr_hotspots) && nrow(chr_hotspots) > 0) {
      overlaps <- chr_hotspots[
        chr_hotspots$hotspot_start_1based <= variants$variant_pos[[i]] &
          chr_hotspots$hotspot_end_1based >= variants$variant_pos[[i]],
        ,
        drop = FALSE
      ]
    }

    variant_base <- cbind(
      variants[i, , drop = FALSE],
      exact_scores[i, , drop = FALSE],
      stringsAsFactors = FALSE
    )

    if (is.null(overlaps) || nrow(overlaps) == 0) {
      overlap_stub <- data.frame(
        in_sv_hotspot = FALSE,
        overlapping_hotspot_count = 0L,
        hotspot_index = NA_integer_,
        hotspot_chr_raw = NA_character_,
        hotspot_chr_norm = NA_character_,
        hotspot_start_0based = NA_integer_,
        hotspot_start_1based = NA_integer_,
        hotspot_end_1based = NA_integer_,
        hotspot_width_bp = NA_integer_,
        hotspot_score_site_count = NA_integer_,
        hotspot_mean_sv_score = NA_real_,
        hotspot_max_sv_score = NA_real_,
        hotspot_high_score = NA,
        hotspot_best_score_pos = NA_integer_,
        hotspot_best_score_strand = NA_character_,
        hotspot_best_score_guide = NA_character_,
        stringsAsFactors = FALSE
      )

      annotated_rows[[row_counter]] <- cbind(variant_base, overlap_stub, stringsAsFactors = FALSE)
      row_counter <- row_counter + 1L
      next
    }

    overlaps$in_sv_hotspot <- TRUE
    overlaps$overlapping_hotspot_count <- nrow(overlaps)

    keep_cols <- c(
      "in_sv_hotspot",
      "overlapping_hotspot_count",
      "hotspot_index",
      "hotspot_chr_raw",
      "hotspot_chr_norm",
      "hotspot_start_0based",
      "hotspot_start_1based",
      "hotspot_end_1based",
      "hotspot_width_bp",
      "hotspot_score_site_count",
      "hotspot_mean_sv_score",
      "hotspot_max_sv_score",
      "hotspot_high_score",
      "hotspot_best_score_pos",
      "hotspot_best_score_strand",
      "hotspot_best_score_guide"
    )

    for (j in seq_len(nrow(overlaps))) {
      annotated_rows[[row_counter]] <- cbind(
        variant_base,
        overlaps[j, keep_cols, drop = FALSE],
        stringsAsFactors = FALSE
      )
      row_counter <- row_counter + 1L
    }
  }

  do.call(rbind, annotated_rows)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  hotspots <- read_hotspots(args$hotspots)
  scores <- read_scores(args$scores, args$score_column)
  variants <- read_variants(args$variants)

  hotspots_scored <- summarise_hotspot_scores(
    hotspots = hotspots,
    scores = scores,
    score_threshold = args$score_threshold
  )

  exact_scores <- summarise_exact_position_scores(
    variants = variants,
    scores = scores,
    score_threshold = args$score_threshold
  )

  annotated <- annotate_variants(
    variants = variants,
    hotspots_scored = hotspots_scored,
    exact_scores = exact_scores
  )

  utils::write.table(
    annotated,
    file = args$output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  unique_variants_in_hotspots <- unique(annotated$variant_id[annotated$in_sv_hotspot %in% TRUE])
  unique_variants_in_high_score_hotspots <- unique(
    annotated$variant_id[
      annotated$in_sv_hotspot %in% TRUE &
        annotated$hotspot_high_score %in% TRUE
    ]
  )

  message("Loaded ", nrow(variants), " variants.")
  message("Loaded ", nrow(hotspots), " hotspots.")
  message("Loaded ", nrow(scores), " scored positions.")
  message("Variants overlapping hotspots: ", length(unique_variants_in_hotspots))
  message("Variants overlapping high-score hotspots: ", length(unique_variants_in_high_score_hotspots))
  message("Wrote annotations to: ", args$output)
}

main()
