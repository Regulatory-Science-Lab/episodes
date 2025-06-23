#' Kaplan-Meier Fit by Prior Lines
#'
#' Fits a Kaplan-Meier survival curve based on time from treatment start to end,
#' stratified by the number of prior lines of treatment.
#'
#' @param drug_transitions A data frame containing treatment transitions per patient.
#'        Must include `patientid`, `state`, `linestartdate`, and `lineenddate`.
#' @param lines A data frame containing covariates such as `prior_lines` for each patient.
#'        Must include `patientid` and `prior_lines`.
#' @param tumour A string identifier for tumour type (used to name the saved output file).
#' @param save_path Optional. Path to directory where the `.RData` file will be saved.
#'        Default is PREDiCText personal folder.
#'
#' @return The `survfit` object from the Kaplan-Meier fit.
#' @export
km_fit <- function(drug_transitions, lines, tumour = "lung", save_path = "H://PREDiCText//nirupama//weibull_estimates//survfit_km") {
  # Prepare survival data
  os_data <- drug_transitions %>%
    dplyr::arrange(linestartdate) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(event = dplyr::if_else(dplyr::row_number() == dplyr::n() & state == "Death", 1L, 0L)) %>%
    dplyr::summarise(
      linestartdate = min(linestartdate),
      lineenddate   = max(lineenddate),
      event         = max(event),
      .groups       = "drop"
    ) %>%
    dplyr::mutate(time = as.numeric(lineenddate - linestartdate))

  # Add covariates
  lines_covariates <- lines %>%
    dplyr::select(patientid, prior_lines) %>%
    dplyr::distinct()

  os_data <- os_data %>%
    dplyr::left_join(lines_covariates, by = "patientid")

  # Fit survival model
  surv_obj <- survival::Surv(time = os_data$time, event = os_data$event)
  km_fit <- survival::survfit(surv_obj ~ prior_lines, data = os_data)

  # Save output
  filename <- file.path(save_path, paste0("km_fit_", tumour, ".RData"))
  save(km_fit, file = filename)

  return(km_fit)
}
