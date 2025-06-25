#' Summarize Patient Transitions Across States
#'
#' This function computes summary counts of transitions between treatment states
#' (e.g., on-treatment to progression, progression to death), optionally joined
#' with demographic gender data from Flatiron.
#'
#' @param drug_transitions A data frame containing patient state transitions.
#'   Must include at least `patientid`, `linestartdate`, `lineenddate`, and `state`.
#' @param flatiron Logical. If `TRUE`, the function will join Flatiron demographic
#'   data (from a hardcoded path) to annotate patient gender.
#'
#' @return A data frame with total transitions between defined states and
#'   patient sex distribution. Includes:
#' \describe{
#'   \item{N_on_to_prog}{Transitions from on-treatment to progression}
#'   \item{N_on_to_death}{Transitions from on-treatment to death}
#'   \item{N_on_to_off}{Transitions from on-treatment to off-treatment}
#'   \item{N_off_to_prog}{Transitions from off-treatment to progression}
#'   \item{N_off_to_death}{Transitions from off-treatment to death}
#'   \item{N_prog_to_death}{Transitions from progression to death}
#'   \item{N_total}{Total number of unique patients}
#'   \item{N_prop_F}{Proportion of female patients}
#'   \item{N_prop_M}{Proportion of male patients}
#' }
#'
#' @import dplyr tidyr haven
#' @export

state_numbers_summary <- function(drug_transitions, flatiron = TRUE) {

  if (flatiron) {
  flatiron_data <- haven::read_dta("H://PREDiCText//lingyi//Flatiron_exploratory_round1and4//Demographics//Demographics_allDataSets_round1and4.dta")

  flatiron_data <- flatiron_data %>%
    dplyr::select(patientid, gender) %>%
    dplyr::distinct(patientid, .keep_all = TRUE)

  drug_transitions <- drug_transitions %>%
    dplyr::left_join(flatiron_data, by = "patientid")
  }

  transition_counts <- drug_transitions %>%
    dplyr::arrange(patientid, linestartdate, lineenddate) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(next_state = dplyr::lead(state)) %>%
    dplyr::ungroup() %>%
    dplyr::count(state, next_state, sort = TRUE) %>%
    dplyr::arrange(state, next_state)

  summary_counts <- transition_counts %>%
    dplyr::mutate(
      label = dplyr::case_when(
        state == "On_Treatment_Target_Line" & next_state == "progression" ~ "N_on_to_prog",
        state == "On_Treatment_Target_Line" & next_state == "Death" ~ "N_on_to_death",
        state == "On_Treatment_Target_Line" & next_state == "off_treatment" ~ "N_on_to_off",
        state == "off_treatment" & next_state == "progression" ~ "N_off_to_prog",
        state == "off_treatment" & next_state == "Death" ~ "N_off_to_death",
        state == "progression" &  next_state == "Death" ~ "N_prog_to_death",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(label))

  summary_table <- summary_counts %>%
    dplyr::select(label, n) %>%
    tidyr::pivot_wider(names_from = label, values_from = n)

  # Overall count by proportion of sex
  n_sex_prop <- drug_transitions %>%
    dplyr::distinct(patientid, .keep_all = TRUE) %>%
    dplyr::group_by(gender) %>%
    dplyr::summarise(n = n()/length(unique(drug_transitions$patientid)))

  summary_table$N_total <- length(unique(drug_transitions$patientid))
  summary_table$N_prop_F <- n_sex_prop$n[n_sex_prop$gender == "F"]
  summary_table$N_prop_M <- n_sex_prop$n[n_sex_prop$gender == "M"]

  desired_cols <- c(
    "N_on_to_death", "N_on_to_off", "N_on_to_prog",
    "N_off_to_death", "N_off_to_prog", "N_prog_to_death",
    "N_total", "N_prop_F", "N_prop_M"
  )

  summary_table <- summary_table %>%
    tibble::as_tibble()

  for (col in desired_cols) {
    if (!col %in% names(summary_table)) {
      summary_table[[col]] <- 0
    }
  }

  summary_table <- summary_table %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      !!!purrr::map(setdiff(desired_cols, names(.)), ~ 0)
    ) %>%
    dplyr::select(all_of(desired_cols))

  return(summary_table)
}
