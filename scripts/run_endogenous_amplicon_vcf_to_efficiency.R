#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  "scripts/run_endogenous_amplicon_vcf_to_efficiency.R"
}

source(file.path(
  dirname(normalizePath(script_file, mustWork = FALSE)),
  "run_amplicon_vcf_to_efficiency.R"
))
