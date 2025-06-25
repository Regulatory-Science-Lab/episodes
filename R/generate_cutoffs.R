#' Generate Cutoff Values for Treatment Gaps, Progression, and Death
#'
#' This function computes data-driven cutoff thresholds for defining:
#' - the gap between treatment lines,
#' - the time to progression after chemotherapy, and
#' - the time to death in patients with only one line and a known progression.
#'
#' These cutoffs are useful for labeling transitions (e.g., censoring vs. death vs. progression)
#' in real-world oncology datasets.
#'
#' @param drug_episodes A data.frame of drug episodes, including columns:
#'        \code{linestartdate}, \code{lineenddate}, \code{progression_date},
#'        \code{death_date}, and \code{chemo_line}.
#'
#' @return A list with three numeric cutoff values:
#' \describe{
#'   \item{\code{next_line_cutoff}}{Median number of days to next line of treatment}
#'   \item{\code{progression_cutoff}}{75th percentile of time to progression}
#'   \item{\code{death_cutoff}}{75th percentile of time to death among single-line patients}
#' }
#'
#' @details
#' The function calculates:
#' \itemize{
#'   \item The **next line cutoff** as the 50th percentile (median) of days between end of one line and start of next.
#'   \item The **progression cutoff** as the 75th percentile of time to progression from end of chemo line.
#'   \item The **death cutoff** as the 75th percentile of time to death, calculated only for patients with a single line and a known progression date.
#' }
#'
#' These thresholds can be used as rule-based criteria in modeling progression-free survival or treatment state transitions.
#'
#' @importFrom dplyr group_by ungroup filter mutate lead n
#' @importFrom stats quantile
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   cutoffs <- generate_cutoffs(drug_episodes)
#'   cutoffs$progression_cutoff
#' }
generate_cutoffs <- function(drug_episodes) {
  drug_episodes <- drug_episodes %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(
      next_linestart = dplyr::lead(linestartdate),
      days_to_next_line = as.numeric(next_linestart - lineenddate)
    ) %>%
    dplyr::ungroup()

  next_line_cutoff <- stats::quantile(drug_episodes$days_to_next_line, 0.50, names = FALSE, na.rm = TRUE)

  drug_episodes_progression <- drug_episodes %>%
    dplyr::group_by(patientid) %>%
    dplyr::filter(linenumber == chemo_line) %>%
    dplyr::mutate(
      days_to_progression = as.numeric(progression_date - lineenddate)
    ) %>%
    dplyr::filter(days_to_progression > 0)

  progression_cutoff <- stats::quantile(drug_episodes_progression$days_to_progression, 0.75, names = FALSE, na.rm = TRUE)

  drug_episodes_death <- drug_episodes %>%
    dplyr::group_by(patientid) %>%
    dplyr::filter(dplyr::n() == 1, linenumber == chemo_line, !is.na(progression_date)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      days_to_death = as.numeric(death_date - lineenddate)
    ) %>%
    dplyr::filter(days_to_death > 0)

  death_cutoff <- stats::quantile(drug_episodes_death$days_to_death, 0.75, names = FALSE, na.rm = TRUE)

  return(list(
    next_line_cutoff = next_line_cutoff,
    progression_cutoff = progression_cutoff,
    death_cutoff = death_cutoff
  ))
}
