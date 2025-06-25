#' Generate a Monthly State and Death Transition Table
#'
#' Creates a monthly summary table showing number of patients in different states
#' stratified by gender, and includes death transitions per state.
#'
#' @param drug_transitions A data frame with patient state transitions.
#' @param flatiron Logical, whether to join with Flatiron demographics (default = TRUE).
#'
#' @return A data frame summarizing transitions by month relative to treatment.
#' @export
#'
#' @import dplyr tidyr lubridate haven
#'
death_table <- function(drug_transitions, flatiron = TRUE) {

  if (flatiron) {
    flatiron_data <- haven::read_dta("H://PREDiCText//lingyi//Flatiron_exploratory_round1and4//Demographics//Demographics_allDataSets_round1and4.dta") %>%
      dplyr::select(patientid, gender) %>%
      dplyr::distinct(patientid, .keep_all = TRUE)
  }

  transition_counts <- drug_transitions %>%
    dplyr::arrange(patientid, linestartdate, lineenddate) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(next_state = dplyr::lead(state)) %>%
    dplyr::ungroup()

  transition_counts_summary <- transition_counts %>%
    dplyr::mutate(month = lubridate::floor_date(linestartdate, "month")) %>%
    dplyr::select(patientid, state, month) %>%
    dplyr::arrange(patientid, month) %>%
    dplyr::group_by(patientid) %>%
    tidyr::complete(month = seq.Date(min(month), max(month), by = "1 month")) %>%
    tidyr::fill(state, .direction = "down") %>%
    dplyr::mutate(
      death_flag = state == "Death",
      after_death = cumsum(death_flag) > 1
    ) %>%
    dplyr::filter(!after_death) %>%
    dplyr::ungroup() %>%
    dplyr::filter(state != "Death") %>%
    dplyr::group_by(patientid, month) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup()

  if (flatiron) {
    transition_counts_summary <- transition_counts_summary %>%
      dplyr::left_join(flatiron_data, by = "patientid") %>%
      dplyr::mutate(gender = dplyr::na_if(trimws(gender), ""))
  }

  transition_counts_summary <- transition_counts_summary %>%
    dplyr::arrange(patientid, month) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(month_relative = dplyr::row_number()) %>%
    dplyr::ungroup()

  summary_counts <- transition_counts_summary %>%
    dplyr::mutate(
      label = dplyr::case_when(
        state == "On_Treatment_Target_Line" & gender == "F" ~ "N_on_F",
        state == "On_Treatment_Target_Line" & gender == "M" ~ "N_on_M",
        state == "off_treatment" & gender == "F" ~ "N_off_F",
        state == "off_treatment" & gender == "M" ~ "N_off_M",
        state == "progression" & gender == "F" ~ "N_p_F",
        state == "progression" & gender == "M" ~ "N_p_M",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(label))

  death_table <- summary_counts %>%
    dplyr::group_by(month_relative, label) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = label, values_from = n, values_fill = 0) %>%
    dplyr::arrange(month_relative)

  # Deaths
  death_transitions <- transition_counts %>%
    dplyr::group_by(patientid) %>%
    dplyr::arrange(linestartdate) %>%
    dplyr::mutate(on_treatment = dplyr::first(linestartdate)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(next_state == "Death") %>%
    dplyr::mutate(
      death_month = round(lubridate::time_length(
        lubridate::interval(on_treatment, lineenddate), "month"
      )) + 1
    ) %>%
    dplyr::select(patientid, state, death_month)

  if (flatiron) {
    death_transitions <- death_transitions %>%
      dplyr::left_join(flatiron_data, by = "patientid")
  }

  summary_counts_death <- death_transitions %>%
    dplyr::mutate(
      label = dplyr::case_when(
        state == "On_Treatment_Target_Line" ~ "N_on_D",
        state == "off_treatment" ~ "N_off_D",
        state == "progression" ~ "N_p_D",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(label))

  summary_table_death <- summary_counts_death %>%
    dplyr::group_by(death_month, label) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = label, values_from = n, values_fill = 0)

  # Join death counts
  death_table <- death_table %>%
    dplyr::left_join(summary_table_death, by = c("month_relative" = "death_month")) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ tidyr::replace_na(., 0)))

  # Order columns
  desired_cols <- c("month_relative", "N_off_F", "N_on_F", "N_on_M", "N_p_F", "N_off_M", "N_p_M",
                   "N_off_D", "N_on_D", "N_p_D")

  death_table <- death_table %>%
    tibble::as_tibble()

  for (col in desired_cols) {
    if (!col %in% names(death_table)) {
      death_table[[col]] <- 0
    }
  }


  death_table <- death_table %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      !!!purrr::map(setdiff(desired_cols, names(.)), ~ 0)
    ) %>%
    dplyr::select(all_of(desired_cols))


  return(death_table)
}
