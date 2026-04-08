#!/usr/bin/env Rscript

script_file <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_file[grepl("^--file=", script_file)][1])
script_dir <- if (!is.na(script_file) && nzchar(script_file)) dirname(normalizePath(script_file, mustWork = FALSE)) else getwd()

source(file.path(script_dir, "..", "R", "common.R"))
source(file.path(script_dir, "..", "R", "yeast_variant_annotation.R"))

load_paper_packages(c("dplyr", "readr", "stringr", "tibble"))

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/annotate_yeast_variants.R",
      "--variants vars.tsv",
      "--genome-fasta sacCer3.fa",
      "--gff sacCer3.gff3",
      "--promoter-bed scepd_promoters.bed",
      "--output vars.annotated.tsv",
      "--snpeff-db sacCer3",
      "[--snpeff-jar /path/to/snpEff.jar]",
      "[--snpeff-command snpEff]",
      "[--snpeff-config /path/to/snpEff.config]",
      "[--snpeff-data-dir /path/to/data]",
      "[--java-opts -Xmx4g]",
      "[--skip-snpeff]",
      sep = " "
    ),
    "\n\n",
    "Expected variant columns:\n",
    "- Headered file with CHR,pos,var,ref,alt (case-insensitive), or\n",
    "- Headerless file where the first 5 columns are chr pos var ref alt.\n\n",
    "Output columns include:\n",
    "- SpCas9 targetability relative to NGG PAMs.\n",
    "- CDS/promoter/intergenic region class.\n",
    "- Flattened SnpEff ANN fields.\n\n",
    "Notes:\n",
    "- All references must be from the same yeast assembly.\n",
    "- Promoters are taken from a BED file derived from scEPDnew TSS entries.\n",
    "- Region class precedence is coding > promoter > intergenic.\n",
    sep = ""
  )
}

parse_args <- function(args) {
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = 0)
  }

  parsed <- list(
    java_opts = "-Xmx4g",
    snpeff_command = "snpEff",
    skip_snpeff = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]

    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }

    key <- gsub("-", "_", sub("^--", "", key))

    if (key %in% c("skip_snpeff")) {
      parsed[[key]] <- TRUE
      i <- i + 1
      next
    }

    if (i == length(args)) {
      stop("Missing value for argument: --", key, call. = FALSE)
    }

    parsed[[key]] <- args[[i + 1]]
    i <- i + 2
  }

  required <- c("variants", "genome_fasta", "gff", "promoter_bed", "output")
  missing <- required[!required %in% names(parsed)]

  if (length(missing) > 0) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (!isTRUE(parsed$skip_snpeff) && !"snpeff_db" %in% names(parsed)) {
    stop("Missing required argument: --snpeff-db (or use --skip-snpeff).", call. = FALSE)
  }

  parsed
}

variant_file_has_header <- function(path) {
  fields <- normalize_annotation_colname(first_line_fields(path))
  any(fields %in% c("chr", "chrom", "chromosome"))
}

read_variants <- function(path) {
  require_file(path, "variant table")

  has_header <- variant_file_has_header(path)
  df <- read_delim_base(path, header = has_header)

  if (!has_header) {
    if (ncol(df) < 5) {
      stop(
        "Headerless variant file must have at least 5 columns: chr pos var ref alt.",
        call. = FALSE
      )
    }

    names(df)[1:5] <- c("chr", "pos", "var", "ref", "alt")
  }

  chr_idx <- pick_column(df, c("chr", "chrom", "chromosome"))
  pos_idx <- pick_column(df, c("pos", "position", "site"))
  var_idx <- pick_column(df, c("var", "variant", "name", "id"))
  ref_idx <- pick_column(df, c("ref", "reference"))
  alt_idx <- pick_column(df, c("alt", "alternate", "alternative"))

  if (any(is.na(c(chr_idx, pos_idx, ref_idx, alt_idx)))) {
    stop("Could not identify chr/pos/ref/alt columns in the variant file.", call. = FALSE)
  }

  variants <- tibble::as_tibble(df) %>%
    dplyr::mutate(
      input_row = dplyr::row_number(),
      variant_chr = as.character(.data[[names(df)[[chr_idx]]]]),
      variant_chr_norm = normalize_yeast_chr(variant_chr),
      variant_pos = as.integer(.data[[names(df)[[pos_idx]]]]),
      variant_name = if (!is.na(var_idx)) as.character(.data[[names(df)[[var_idx]]]]) else sprintf("variant_%06d", input_row),
      variant_ref = toupper(as.character(.data[[names(df)[[ref_idx]]]])),
      variant_alt = toupper(as.character(.data[[names(df)[[alt_idx]]]])),
      variant_start = variant_pos,
      variant_end = variant_pos + pmax(nchar(variant_ref), 1L) - 1L,
      variant_width_bp = variant_end - variant_start + 1L,
      variant_id = paste(variant_chr, variant_pos, variant_ref, variant_alt, sep = ":")
    )

  if (any(is.na(variants$variant_chr_norm) | !nzchar(variants$variant_chr_norm))) {
    stop("Some chromosomes could not be normalized to yeast chromosome names.", call. = FALSE)
  }

  if (any(is.na(variants$variant_pos))) {
    stop("Some variant positions could not be parsed as integers.", call. = FALSE)
  }

  if (any(is.na(variants$variant_ref) | !nzchar(variants$variant_ref))) {
    stop("Some REF alleles are missing.", call. = FALSE)
  }

  if (any(is.na(variants$variant_alt) | !nzchar(variants$variant_alt))) {
    stop("Some ALT alleles are missing.", call. = FALSE)
  }

  variants
}

read_cds_features <- function(path) {
  require_file(path, "GFF annotation")

  gff <- utils::read.delim(
    file = path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "#",
    quote = "",
    fill = TRUE
  )

  if (ncol(gff) < 9) {
    stop("GFF file must have at least 9 columns.", call. = FALSE)
  }

  names(gff)[1:9] <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

  cds <- tibble::as_tibble(gff) %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::mutate(
      chr_norm = normalize_yeast_chr(seqid),
      feature_start = as.integer(start),
      feature_end = as.integer(end),
      parsed_attributes = lapply(attributes, parse_gff_attributes),
      feature_id = vapply(
        parsed_attributes,
        function(attrs) {
          first_non_missing(
            pick_gff_attribute(attrs, c("gene", "gene_name", "Name", "locus_tag", "Parent", "ID"))
          )
        },
        character(1)
      )
    ) %>%
    dplyr::select(chr_norm, feature_start, feature_end, strand, feature_id)

  cds$feature_id[is.na(cds$feature_id) | !nzchar(cds$feature_id)] <- paste0(
    cds$chr_norm[is.na(cds$feature_id) | !nzchar(cds$feature_id)],
    ":",
    cds$feature_start[is.na(cds$feature_id) | !nzchar(cds$feature_id)],
    "-",
    cds$feature_end[is.na(cds$feature_id) | !nzchar(cds$feature_id)]
  )

  cds
}

promoter_bed_has_header <- function(path) {
  fields <- normalize_annotation_colname(first_line_fields(path))
  any(fields %in% c("chr", "chrom", "chromosome")) &&
    any(fields %in% c("start", "chromstart")) &&
    any(fields %in% c("end", "chromend"))
}

read_promoter_bed <- function(path) {
  require_file(path, "promoter BED")

  has_header <- promoter_bed_has_header(path)
  df <- read_delim_base(path, header = has_header)

  if (!has_header) {
    if (ncol(df) < 3) {
      stop("Headerless BED file must have at least 3 columns.", call. = FALSE)
    }

    default_names <- c("chrom", "start", "end", "name", "score", "strand")
    names(df)[seq_len(min(ncol(df), length(default_names)))] <- default_names[seq_len(min(ncol(df), length(default_names)))]
  }

  chr_idx <- pick_column(df, c("chr", "chrom", "chromosome"))
  start_idx <- pick_column(df, c("start", "chromstart"))
  end_idx <- pick_column(df, c("end", "chromend"))
  name_idx <- pick_column(df, c("name", "id", "promoterid", "gene", "geneid"))
  strand_idx <- pick_column(df, c("strand"))

  if (any(is.na(c(chr_idx, start_idx, end_idx)))) {
    stop("Could not identify chrom/start/end columns in the promoter BED file.", call. = FALSE)
  }

  promoters <- tibble::as_tibble(df) %>%
    dplyr::mutate(
      promoter_chr = as.character(.data[[names(df)[[chr_idx]]]]),
      chr_norm = normalize_yeast_chr(promoter_chr),
      bed_start = as.integer(.data[[names(df)[[start_idx]]]]),
      bed_end = as.integer(.data[[names(df)[[end_idx]]]]),
      feature_start = bed_start + 1L,
      feature_end = bed_end,
      strand = if (!is.na(strand_idx)) as.character(.data[[names(df)[[strand_idx]]]]) else "*",
      feature_id = if (!is.na(name_idx)) as.character(.data[[names(df)[[name_idx]]]]) else NA_character_
    ) %>%
    dplyr::select(chr_norm, feature_start, feature_end, strand, feature_id)

  promoters$feature_id[is.na(promoters$feature_id) | !nzchar(promoters$feature_id)] <- paste0(
    promoters$chr_norm[is.na(promoters$feature_id) | !nzchar(promoters$feature_id)],
    ":",
    promoters$feature_start[is.na(promoters$feature_id) | !nzchar(promoters$feature_id)],
    "-",
    promoters$feature_end[is.na(promoters$feature_id) | !nzchar(promoters$feature_id)]
  )

  if (any(is.na(promoters$feature_start) | is.na(promoters$feature_end))) {
    stop("Promoter BED start/end columns contain non-integer values.", call. = FALSE)
  }

  if (any(is.na(promoters$chr_norm) | !nzchar(promoters$chr_norm))) {
    stop("Some promoter BED chromosomes could not be normalized to yeast chromosome names.", call. = FALSE)
  }

  promoters
}

build_spcas9_pam_index <- function(fasta) {
  lapply(
    fasta$sequences,
    function(sequence) {
      bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
      n <- length(bases)

      if (n < 3) {
        return(list(plus = integer(), minus = integer()))
      }

      plus <- which(bases[2:(n - 1)] == "G" & bases[3:n] == "G")
      minus <- which(bases[1:(n - 2)] == "C" & bases[2:(n - 1)] == "C")

      list(plus = plus, minus = minus)
    }
  )
}

annotate_spcas9 <- function(variants, pam_index) {
  plus_counts <- integer(nrow(variants))
  minus_counts <- integer(nrow(variants))

  for (i in seq_len(nrow(variants))) {
    chr_norm <- variants$variant_chr_norm[[i]]
    chr_pams <- pam_index[[chr_norm]]

    if (is.null(chr_pams)) {
      stop("Variant chromosome missing from genome FASTA: ", chr_norm, call. = FALSE)
    }

    plus_min <- variants$variant_start[[i]] + 1L
    plus_max <- variants$variant_end[[i]] + 20L
    minus_min <- variants$variant_start[[i]] - 22L
    minus_max <- variants$variant_end[[i]] - 3L

    plus_counts[[i]] <- sum(chr_pams$plus >= plus_min & chr_pams$plus <= plus_max)
    minus_counts[[i]] <- sum(chr_pams$minus >= minus_min & chr_pams$minus <= minus_max)
  }

  tibble::tibble(
    spcas9_ngg_pam_count_plus = plus_counts,
    spcas9_ngg_pam_count_minus = minus_counts,
    spcas9_ngg_pam_count = plus_counts + minus_counts,
    spcas9_targetable = (plus_counts + minus_counts) > 0
  )
}

annotate_interval_hits <- function(variants, features, label_prefix) {
  feature_index <- split(seq_len(nrow(features)), features$chr_norm)
  hit_counts <- integer(nrow(variants))
  hit_ids <- rep(NA_character_, nrow(variants))

  for (i in seq_len(nrow(variants))) {
    chr_hits <- feature_index[[variants$variant_chr_norm[[i]]]]

    if (length(chr_hits) == 0) {
      next
    }

    overlaps <- chr_hits[
      features$feature_start[chr_hits] <= variants$variant_end[[i]] &
        features$feature_end[chr_hits] >= variants$variant_start[[i]]
    ]

    hit_counts[[i]] <- length(overlaps)

    if (length(overlaps) > 0) {
      hit_ids[[i]] <- collapse_unique_values(features$feature_id[overlaps])
    }
  }

  tibble::tibble(
    !!paste0(label_prefix, "_count") := hit_counts,
    !!paste0("overlaps_", label_prefix) := hit_counts > 0,
    !!paste0(label_prefix, "_ids") := hit_ids
  )
}

choose_region_class <- function(overlaps_cds, overlaps_promoter) {
  dplyr::case_when(
    overlaps_cds ~ "coding",
    overlaps_promoter ~ "promoter",
    TRUE ~ "intergenic"
  )
}

split_shell_words <- function(x) {
  pieces <- strsplit(trimws(x), "[[:space:]]+")[[1]]
  pieces[nzchar(pieces)]
}

write_variant_vcf <- function(variants, fasta, path) {
  missing_chr <- setdiff(unique(variants$variant_chr_norm), names(fasta$raw_names))

  if (length(missing_chr) > 0) {
    stop(
      "These variant chromosomes are not present in the FASTA: ",
      paste(missing_chr, collapse = ", "),
      call. = FALSE
    )
  }

  vcf_chr <- unname(fasta$raw_names[variants$variant_chr_norm])

  vcf_lines <- c(
    "##fileformat=VCFv4.2",
    "##source=annotate_yeast_variants.R",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    paste(
      vcf_chr,
      variants$variant_pos,
      paste0("row", variants$input_row),
      variants$variant_ref,
      variants$variant_alt,
      ".",
      "PASS",
      ".",
      sep = "\t"
    )
  )

  writeLines(vcf_lines, con = path)
  invisible(path)
}

snpeff_summary_template <- function(n) {
  tibble::tibble(
    snpeff_ann_count = rep(0L, n),
    snpeff_primary_effect = rep(NA_character_, n),
    snpeff_primary_impact = rep(NA_character_, n),
    snpeff_primary_gene = rep(NA_character_, n),
    snpeff_primary_gene_id = rep(NA_character_, n),
    snpeff_primary_feature_id = rep(NA_character_, n),
    snpeff_primary_hgvs_c = rep(NA_character_, n),
    snpeff_primary_hgvs_p = rep(NA_character_, n),
    snpeff_all_effects = rep(NA_character_, n),
    snpeff_all_impacts = rep(NA_character_, n),
    snpeff_all_genes = rep(NA_character_, n),
    snpeff_all_gene_ids = rep(NA_character_, n),
    snpeff_all_feature_ids = rep(NA_character_, n)
  )
}

summarize_snpeff_ann <- function(ann_df) {
  if (nrow(ann_df) == 0) {
    return(snpeff_summary_template(1))
  }

  primary <- ann_df[1, , drop = FALSE]

  tibble::tibble(
    snpeff_ann_count = nrow(ann_df),
    snpeff_primary_effect = first_non_missing(primary$effect),
    snpeff_primary_impact = first_non_missing(primary$impact),
    snpeff_primary_gene = first_non_missing(primary$gene_name),
    snpeff_primary_gene_id = first_non_missing(primary$gene_id),
    snpeff_primary_feature_id = first_non_missing(primary$feature_id),
    snpeff_primary_hgvs_c = first_non_missing(primary$hgvs_c),
    snpeff_primary_hgvs_p = first_non_missing(primary$hgvs_p),
    snpeff_all_effects = collapse_unique_values(ann_df$effect),
    snpeff_all_impacts = collapse_unique_values(ann_df$impact),
    snpeff_all_genes = collapse_unique_values(ann_df$gene_name),
    snpeff_all_gene_ids = collapse_unique_values(ann_df$gene_id),
    snpeff_all_feature_ids = collapse_unique_values(ann_df$feature_id)
  )
}

run_snpeff <- function(variants, fasta, args) {
  summary_df <- snpeff_summary_template(nrow(variants))

  if (isTRUE(args$skip_snpeff)) {
    return(summary_df)
  }

  temp_dir <- tempfile("snpeff_")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  input_vcf <- file.path(temp_dir, "variants.vcf")
  output_vcf <- file.path(temp_dir, "variants.snpeff.vcf")
  stderr_file <- file.path(temp_dir, "variants.snpeff.stderr.txt")

  write_variant_vcf(variants, fasta, input_vcf)

  if ("snpeff_jar" %in% names(args)) {
    require_file(args$snpeff_jar, "snpEff jar")
    command <- "java"
    command_args <- c(split_shell_words(args$java_opts), "-jar", args$snpeff_jar, "ann")
  } else {
    command <- args$snpeff_command

    if (!nzchar(Sys.which(command))) {
      stop(
        "Could not find `", command, "` in PATH. Provide --snpeff-jar or install snpEff.",
        call. = FALSE
      )
    }

    command_args <- c("ann")
  }

  if ("snpeff_config" %in% names(args)) {
    require_file(args$snpeff_config, "snpEff config")
    command_args <- c(command_args, "-c", args$snpeff_config)
  }

  if ("snpeff_data_dir" %in% names(args)) {
    require_dir(args$snpeff_data_dir, "snpEff data directory")
    command_args <- c(command_args, "-dataDir", args$snpeff_data_dir)
  }

  command_args <- c(
    command_args,
    "-noStats",
    "-noLog",
    "-nodownload",
    "-ud",
    "0",
    args$snpeff_db,
    input_vcf
  )

  exit_status <- system2(
    command = command,
    args = command_args,
    stdout = output_vcf,
    stderr = stderr_file
  )

  stderr_output <- if (file.exists(stderr_file)) readLines(stderr_file, warn = FALSE) else character()
  status <- exit_status %||% 0

  if (status != 0) {
    stop(
      "snpEff failed with status ", status, ":\n",
      paste(stderr_output, collapse = "\n"),
      call. = FALSE
    )
  }

  if (!file.exists(output_vcf)) {
    stop(
      "snpEff finished without producing an output VCF.\n",
      paste(stderr_output, collapse = "\n"),
      call. = FALSE
    )
  }

  body_lines <- readLines(output_vcf, warn = FALSE)
  body_lines <- body_lines[!startsWith(body_lines, "#")]

  if (length(body_lines) == 0) {
    return(summary_df)
  }

  snpeff_map <- vector("list", nrow(variants))
  names(snpeff_map) <- paste0("row", variants$input_row)

  for (line in body_lines) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]

    if (length(fields) < 8) {
      next
    }

    record_id <- fields[[3]]
    info_field <- fields[[8]]
    ann_string <- parse_info_value(info_field, "ANN")
    snpeff_map[[record_id]] <- summarize_snpeff_ann(parse_snpeff_ann(ann_string))
  }

  summaries <- lapply(
    paste0("row", variants$input_row),
    function(record_id) {
      snpeff_map[[record_id]] %||% snpeff_summary_template(1)
    }
  )

  dplyr::bind_rows(summaries)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  message("Reading variants: ", args$variants)
  variants <- read_variants(args$variants)

  message("Reading genome FASTA: ", args$genome_fasta)
  fasta <- read_fasta_sequences(args$genome_fasta)

  message("Scanning SpCas9 NGG PAM sites across the genome.")
  pam_index <- build_spcas9_pam_index(fasta)
  spcas9_annotations <- annotate_spcas9(variants, pam_index)

  message("Reading CDS annotations: ", args$gff)
  cds_features <- read_cds_features(args$gff)

  message("Reading scEPD-derived promoter BED: ", args$promoter_bed)
  promoter_features <- read_promoter_bed(args$promoter_bed)

  cds_annotations <- annotate_interval_hits(variants, cds_features, "cds")
  promoter_annotations <- annotate_interval_hits(variants, promoter_features, "promoter")

  snpeff_annotations <- run_snpeff(variants, fasta, args)

  annotated <- dplyr::bind_cols(
    variants,
    spcas9_annotations,
    cds_annotations,
    promoter_annotations,
    snpeff_annotations
  ) %>%
    dplyr::mutate(
      region_class = choose_region_class(overlaps_cds, overlaps_promoter)
    )

  ensure_dir(dirname(args$output))
  readr::write_tsv(annotated, args$output, na = "")

  message("Wrote annotated table: ", normalizePath(args$output, mustWork = FALSE))
  message("Variants: ", nrow(annotated))
  message("SpCas9 targetable: ", sum(annotated$spcas9_targetable))
  message("Coding: ", sum(annotated$region_class == "coding"))
  message("Promoter: ", sum(annotated$region_class == "promoter"))
  message("Intergenic: ", sum(annotated$region_class == "intergenic"))
}

main()
