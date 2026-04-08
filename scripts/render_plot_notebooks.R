notebooks <- c(
  "notebooks/plots/01_plot_gfp_library_transfer_counts.qmd",
  "notebooks/plots/02_plot_enas_amplicon_colony_counts.qmd",
  "notebooks/plots/03_plot_donor_location.qmd",
  "notebooks/plots/04_plot_cas_variants_overview.qmd",
  "notebooks/plots/05_plot_cas_variants_viability.qmd"
)

for (notebook in notebooks) {
  message("Rendering ", notebook)
  status <- system2("quarto", c("render", notebook))
  if (!identical(status, 0L)) {
    stop("Failed while rendering ", notebook, call. = FALSE)
  }
}
