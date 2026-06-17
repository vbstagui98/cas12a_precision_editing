default_analysis_paths <- function() {
  project_root <- normalizePath(getwd(), mustWork = FALSE)

  list(
    cas12a_project_root = file.path(project_root, "data", "cas12a_precision_editing"),
    figures_root = file.path(project_root, "figures"),
    attune_root = file.path(project_root, "data", "attune"),
    thesis_root = file.path(project_root, "data", "thesis"),
    genome_wide_root = file.path(project_root, "data", "genome_wide"),
    presentations_amplicon_root = file.path(project_root, "data", "amplicon_seq"),
    amplicon_root = file.path(project_root, "data", "amplicon_seq"),
    cognate_dr_root = file.path(project_root, "data", "cognate_dr"),
    short_guide_root = file.path(project_root, "data", "short_guide"),
    downloads_root = file.path(project_root, "data", "downloads"),
    downloads_cas_variants_root = file.path(project_root, "data", "downloads", "cas12a_variants"),
    results_dir = normalizePath(file.path(project_root, "results"), mustWork = FALSE)
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
