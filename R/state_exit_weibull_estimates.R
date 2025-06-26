#' Estimate Weibull parameters for state transitions
#'
#' Fits Weibull models for different treatment states (e.g., on treatment, off treatment, progression)
#' using survival time and covariates such as line of therapy.
#'
#' @param drug_transitions A dataframe of treatment episodes with start/end dates and states.
#' @param lines A dataframe of treatment line data including prior lines and next line dates.
#'
#' @return A tibble with Weibull shape and scale estimates for each transition state.
#'
#' @importFrom dplyr filter mutate left_join distinct group_by ungroup bind_rows select arrange desc
#' @importFrom tibble tibble
#' @importFrom survival Surv
#' @importFrom flexsurv flexsurvreg
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom stats as.formula
#' @export
state_exit_weibull_estimates <- function(drug_transitions, lines) {

  get_surv_data <- function(state_filter, join_LoT = FALSE) {
    df <- drug_transitions %>%
      dplyr::filter(state == state_filter) %>%
      dplyr::mutate(time = as.numeric(lineenddate - linestartdate)) %>%
      dplyr::filter(time > 0)

    df <- df %>%
      dplyr::left_join(
        lines %>% dplyr::select(patientid, prior_lines) %>% dplyr::distinct(),
        by = "patientid"
      )

    if (join_LoT) {
      LoT_progression <- lines %>%
        dplyr::group_by(patientid) %>%
        dplyr::mutate(new_line = ifelse(!is.na(next_linestart), "LoT", "clinical")) %>%
        dplyr::distinct(patientid, .keep_all = TRUE) %>%
        dplyr::select(patientid, new_line)

      df <- df %>% dplyr::left_join(LoT_progression, by = "patientid")
    }
    return(df)
  }

  fit_and_extract <- function(df, state_label, covariates = NULL) {
    formula <- if (!is.null(covariates)) {
      stats::as.formula(paste("Surv(time, event) ~", paste(covariates, collapse = "+")))
    } else {
      survival::Surv(time, event) ~ 1
    }

    n_events <- sum(df$event == 1, na.rm = TRUE)
    if (n_events == 0) {
      message(glue::glue("No observed events for transition {state_label}. Returning NA."))
      return(tibble::tibble(
        state_transition = state_label,
        shape = NA_real_, shape_lci = NA_real_, shape_uci = NA_real_,
        scale = NA_real_, scale_lci = NA_real_, scale_uci = NA_real_
      ))
    }

    fit <- tryCatch({
      flexsurv::flexsurvreg(formula, data = df, dist = "weibull")
    }, error = function(e) NULL)

    if (is.null(fit)) {
      covariates <- "prior_lines"

      # Drop covariate if it doesn't have >1 level
      covariates <- covariates[sapply(covariates, function(var) {
        n_vals <- length(unique(df[[var]]))
        n_vals > 1 && !all(is.na(df[[var]]))
      })]

      formula <- stats::as.formula(
        paste("Surv(time, event) ~", if (length(covariates) > 0) paste(covariates, collapse = "+") else "1")
      )

      fit <- tryCatch({
        flexsurv::flexsurvreg(formula, data = df, dist = "weibull")
      }, error = function(e) NULL)
    }

    if (!is.null(fit)) {
      tibble::tibble(
        state_transition = state_label,
        shape = fit$res["shape", "est"],
        shape_lci = fit$res["shape", "L95%"],
        shape_uci = fit$res["shape", "U95%"],
        scale = fit$res["scale", "est"],
        scale_lci = fit$res["scale", "L95%"],
        scale_uci = fit$res["scale", "U95%"]
      )
    } else {
      message(glue::glue("Weibull model failed for {state_label}. Returning NA."))
      tibble::tibble(
        state_transition = state_label,
        shape = NA_real_, shape_lci = NA_real_, shape_uci = NA_real_,
        scale = NA_real_, scale_lci = NA_real_, scale_uci = NA_real_
      )
    }
  }

  results <- dplyr::bind_rows(
    fit_and_extract(get_surv_data("On_Treatment_Target_Line"), "On_Treatment_Target_Line", covariates = c("prior_lines")),
    fit_and_extract(get_surv_data("off_treatment", join_LoT = TRUE), "off_treatment", covariates = c("prior_lines", "new_line")),
    fit_and_extract(get_surv_data("progression", join_LoT = TRUE), "progression", covariates = c("prior_lines", "new_line"))
  )

  results <- results %>%
    dplyr::arrange(dplyr::desc(state_transition == "On_Treatment_Target_Line"), state_transition)

  results <- na.omit(results)

  return(results)
}
