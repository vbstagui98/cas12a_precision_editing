notebooks <- c(
  "notebooks/01_gfp_library_transfer_counts.qmd",
  "notebooks/02_enas_amplicon_colony_counts.qmd",
  "notebooks/03_donor_location.qmd",
  "notebooks/04_cas_variants_overview.qmd",
  "notebooks/05_cas_variants_viability.qmd",
  "notebooks/06_genome_wide_design_rules.qmd"
)

for (notebook in notebooks) {
  message("Rendering ", notebook)
  status <- system2("quarto", c("render", notebook))
  if (!identical(status, 0L)) {
    stop("Failed while rendering ", notebook, call. = FALSE)
  }
}
