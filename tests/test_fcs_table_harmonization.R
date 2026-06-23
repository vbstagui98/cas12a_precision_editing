#!/usr/bin/env Rscript

source("R/fcs_efficiency.R")

platemap <- tibble::tibble(
  name = c("sample.fcs", "control.fcs"),
  WellPosition = c("A01", "A09"),
  assay = c("GFP_OFF", "GFP_OFF"),
  crRNA_id = c("G03", "control"),
  replicate = c("D1", "D1"),
  timepoint = c("T1", "T1")
)
features <- tibble::tibble(
  crRNA_id = "G03",
  deepcpf1_score = 0.4609433
)
colony_counts <- tibble::tibble(
  name = c("sample.fcs", "control.fcs"),
  colony_count = c(50, 100)
)
pop_stats <- tibble::tibble(
  name = c("sample.fcs", "control.fcs"),
  Population = c("GFP+", "GFP+"),
  Parent = c("root", "root"),
  Frequency = c(0.2, 0.1),
  ParentFrequency = c(1, 1)
)

platemap_path <- tempfile(fileext = ".csv")
features_path <- tempfile(fileext = ".csv")
colony_counts_path <- tempfile(fileext = ".csv")
readr::write_csv(platemap, platemap_path)
readr::write_csv(features, features_path)
readr::write_csv(colony_counts, colony_counts_path)

out <- annotate_fcs_efficiency(
  pop_stats = pop_stats,
  platemap_path = platemap_path,
  guide_features_path = features_path,
  colony_counts_path = colony_counts_path,
  control_well = "A09"
)

unlink(c(platemap_path, features_path, colony_counts_path))

sample_row <- out[out$name == "sample.fcs", ]
stopifnot(
  sample_row$guide_id[[1]] == "G03",
  sample_row$deepcpf1_score[[1]] == 0.4609433,
  sample_row$efficiency_pct[[1]] == 80,
  sample_row$colony_count[[1]] == 50,
  sample_row$control_colony_count[[1]] == 100,
  sample_row$relative_viability[[1]] == 0.5,
  sample_row$viability_normalization_sample[[1]] == "control well A09"
)

cat("FCS table harmonization tests passed\n")
