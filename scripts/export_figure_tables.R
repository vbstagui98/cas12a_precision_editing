script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project_root <- if (!is.na(script_file) && nzchar(script_file)) {
  normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = FALSE)
}

source(file.path(project_root, "config", "paths.R"))
source(file.path(project_root, "R", "common.R"))
source(file.path(project_root, "R", "figure_tables.R"))

load_paper_packages()

paths <- get_analysis_paths(
  list(results_dir = normalizePath(file.path(project_root, "results"), mustWork = FALSE))
)

output <- export_publication_figure_tables(paths)

message("Wrote publication figure source tables to: ", output$output_dir)
print(output$summary)
