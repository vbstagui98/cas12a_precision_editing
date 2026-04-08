default_analysis_paths <- function() {
  list(
    cas12a_project_root = "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Writing/Cas12a_precision_editing",
    figures_root = "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Writing/Cas12a_precision_editing/Figures",
    attune_root = "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/ATTUNE_DATA",
    thesis_root = "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/Master_thesis/NGS_library",
    genome_wide_root = "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide",
    presentations_amplicon_root = "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Presentations/9-month_figs/amplicon_seq",
    amplicon_root = "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/amplicon_seq",
    cognate_dr_root = "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/dilutions_20250120_cognate_DR",
    short_guide_root = "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/20250219_ys85ys88_shorter_guide",
    downloads_root = "/Users/u0174312/Downloads",
    downloads_cas_variants_root = "/Users/u0174312/Downloads/20240724_cas12aversions",
    results_dir = normalizePath(file.path(getwd(), "results"), mustWork = FALSE)
  )
}

get_analysis_paths <- function(overrides = list()) {
  env_lookup <- c(
    cas12a_project_root = "CAS12A_PROJECT_ROOT",
    figures_root = "CAS12A_FIGURES_ROOT",
    attune_root = "CAS12A_ATTUNE_ROOT",
    thesis_root = "CAS12A_THESIS_ROOT",
    genome_wide_root = "CAS12A_GENOME_WIDE_ROOT",
    presentations_amplicon_root = "CAS12A_PRESENTATIONS_AMPLICON_ROOT",
    amplicon_root = "CAS12A_AMPLICON_ROOT",
    cognate_dr_root = "CAS12A_COGNATE_DR_ROOT",
    short_guide_root = "CAS12A_SHORT_GUIDE_ROOT",
    downloads_root = "CAS12A_DOWNLOADS_ROOT",
    downloads_cas_variants_root = "CAS12A_DOWNLOADS_VARIANTS_ROOT",
    results_dir = "CAS12A_RESULTS_DIR"
  )

  paths <- default_analysis_paths()

  for (name in names(env_lookup)) {
    env_value <- Sys.getenv(env_lookup[[name]], unset = "")
    if (nzchar(env_value)) {
      paths[[name]] <- env_value
    }
  }

  utils::modifyList(paths, overrides)
}
