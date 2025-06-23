#' Run pipeline with custom tumour definitions
#'
#' @param tumours A tibble with columns "tumour", "treatment", and "exact"
#' @param dir Path to directory where the RDS file should be saved.
#'
#' @export
run_pipeline <- function(tumours, dir = "") {
  stopifnot(all(c("tumour", "treatment", "exact") %in% colnames(tumours)))

  # Create the directory if it doesn't exist
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  # Save tumour definitions
  tumour_path <- file.path(dir, "tumour_defs.rds")
  saveRDS(tumours, tumour_path)

  # Set env var so _targets.R can read it
  Sys.setenv(TUMOUR_PATH = normalizePath(tumour_path))

  # Run the pipeline
  targets::tar_make()
}
