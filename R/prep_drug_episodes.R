#' Prepare Drug Episode Data from Flatiron drug episodes data
#'
#' This function loads, filters, and processes drug episode data for cancer patients,
#' focusing on advanced or metastatic treatments. It cleans and
#' merges treatment episodes, filters based on the specified treatment given at line 1, 2 or 3, and integrates
#' progression, mortality and last contact data.
#'
#' @param tumour Character. Tumour type to load. Must match the drug episodes files in the Flatiron exploratory round folder. (e.g., "nsclc").
#' @param treatment Character. Regex string of drug names to filter for (e.g., "cisplatin|carboplatin").
#' @param drug_separator A character string separating multiple drugs in `linename`. Default is ",".
#' @param overlap_threshold Integer: number of shared drugs required to collapse lines. Default = 1.
#' @param exact boolean: to specify whether we want an exact match to the drug name or a combination therapy
#'
#' @return A cleaned and filtered `data.frame` of drug episodes with progression, death, and follow-up dates integrated
#' @importFrom dplyr filter mutate select summarise group_by ungroup arrange left_join distinct pull lead
#' @importFrom tidyr drop_na
#' @importFrom stringr str_detect
#' @importFrom lubridate ymd year
#' @importFrom haven read_dta
#' @importFrom glue glue
#' @importFrom knitr kable
#' @importFrom dplyr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   episodes <- prep_episode_data(tumour = "nsclc", treatment = "cisplatin|carboplatin")
#'   head(episodes)
#' }

prep_episode_data <- function(tumour = "nsclc", treatment = "cisplatin|carboplatin", drug_separator = ",",
                              overlap_threshold = 1, exact = FALSE) {
  if(tumour == "thyroid") {
   tumour_data <- read.csv("H://PREDiCText//nirupama//weibull_estimates//thyroid_episodes.csv")
   drug_episodes <- tumour_data %>%
     tidyr::drop_na(linenumber) %>%
     dplyr::group_by(patientid, linenumber) %>%
     dplyr::summarise(
       linename = dplyr::first(linename),
       linestartdate = min(linestartdate),
       lineenddate = max(episodedate),
       .groups = 'drop'
     ) %>%
     dplyr::distinct() %>%
     dplyr::arrange(patientid, linenumber) %>%
     dplyr::mutate(
       linestartdate = lubridate::ymd(linestartdate),
       lineenddate   = lubridate::ymd(lineenddate)
     )

  } else{
  # Load round 1 and round 4 drug episode data
  tumour_ep_round1 <- haven::read_dta(glue::glue("H:\\PREDiCText\\lingyi\\Flatiron_exploratory_round1and4\\DrugEpisodes\\DrugEpisode_{tumour}_round1.dta"))
  tumour_ep_round4 <- haven::read_dta(glue::glue("H:\\PREDiCText\\lingyi\\Flatiron_exploratory_round1and4\\DrugEpisodes\\DrugEpisode_{tumour}_round4.dta"))

  tumour_data <- dplyr::bind_rows(tumour_ep_round1, tumour_ep_round4)
  n_all <- length(unique(tumour_data$patientid))
  message(glue::glue("The total number of {tumour} patients is {n_all}"))

  # Filter for advanced/metastatic lines and extract relevant info
  drug_episodes <- tumour_data %>%
    dplyr::filter(stringr::str_detect(tolower(linesetting), "advanced|metastatic")) %>%
    dplyr::select(patientid, linenumber, linename, linesetting,
                  linestartdate, lineenddate, episodedate, detaileddrugcategory) %>%
    tidyr::drop_na(linenumber) %>%
    dplyr::group_by(patientid, linenumber) %>%
    dplyr::summarise(
      linename = dplyr::first(linename),
      linestartdate = min(linestartdate),
      lineenddate = max(episodedate),
      .groups = 'drop'
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(patientid, linenumber) %>%
    dplyr::mutate(
      linestartdate = lubridate::ymd(linestartdate),
      lineenddate   = lubridate::ymd(lineenddate)
    )
  }

  n_adv <- length(unique(drug_episodes$patientid))
  message(glue::glue("The total number of {tumour} patients is {n_adv}"))

  # First-line patients with specified treatment
  if (exact) {
    ep_patients_first_line <- drug_episodes %>%
      dplyr::filter(tolower(linename) == tolower(treatment), linenumber == 1) %>%
      dplyr::pull(patientid) %>%
      unique()
  } else {
    ep_patients_first_line <- drug_episodes %>%
      dplyr::filter(stringr::str_detect(tolower(linename), treatment), linenumber == 1) %>%
      dplyr::pull(patientid) %>%
      unique()
  }

  drug_episodes_first_line <- drug_episodes %>%
    dplyr::filter(patientid %in% ep_patients_first_line)
  n_first_line <- length(unique(ep_patients_first_line))
  message(glue::glue("The total number of {tumour} first-line patients is {n_first_line}"))

  # Second-line patients
  if (exact) {
    ep_patients_second_line <- drug_episodes %>%
      dplyr::filter(tolower(linename) == tolower(treatment), linenumber == 2) %>%
      dplyr::pull(patientid) %>%
      unique()
  } else {
    ep_patients_second_line <- drug_episodes %>%
      dplyr::filter(stringr::str_detect(tolower(linename), treatment), linenumber == 2) %>%
      dplyr::pull(patientid) %>%
      unique()
  }

  drug_episodes_second_line <- drug_episodes %>%
    dplyr::filter(patientid %in% ep_patients_second_line & !(patientid %in% ep_patients_first_line))
  n_second_line <- length(unique(drug_episodes_second_line$patientid))
  message(glue::glue("The total number of {tumour} second-line patients is {n_second_line}"))

  # Third-line patients
  if (exact) {
    ep_patients_third_line <- drug_episodes %>%
      dplyr::filter(tolower(linename) == tolower(treatment), linenumber == 3) %>%
      dplyr::pull(patientid) %>%
      unique()
  } else {
    ep_patients_third_line <- drug_episodes %>%
      dplyr::filter(stringr::str_detect(tolower(linename), treatment), linenumber == 3) %>%
      dplyr::pull(patientid) %>%
      unique()
  }

  drug_episodes_third_line <- drug_episodes %>%
    dplyr::filter(patientid %in% ep_patients_third_line &
                    !(patientid %in% ep_patients_first_line) &
                    !(patientid %in% ep_patients_second_line))
  n_third_line <- length(unique(drug_episodes_third_line$patientid))
  message(glue::glue("The total number of {tumour} third-line patients is {n_third_line}"))

  # Summary table
  summary_table <- data.frame(
    Line_of_Therapy = c("First Line", "Second Line", "Third Line"),
    Number_of_Patients = c(n_first_line, n_second_line, n_third_line)
  )
  print(knitr::kable(summary_table, caption = glue::glue("Number of patients who received {treatment} by line of therapy")))

  # Combine all selected episodes
  drug_episodes <- dplyr::bind_rows(drug_episodes_first_line, drug_episodes_second_line, drug_episodes_third_line)

  # Remove overlapping lines
  drug_episodes <- collapse_overlapping_lines(ep_data = drug_episodes)
  n_overlap <- length(unique(drug_episodes$patientid))
  message(glue::glue("The total number of {tumour} patients with overlapping lines that were filtered out is {n_overlap}"))

  # Remove single doses
  drug_episodes <- single_dose_clean(ep_data = drug_episodes, drug_separator = drug_separator, overlap_threshold = overlap_threshold)
  n_single_dose <- length(unique(drug_episodes$patientid))
  message(glue::glue("The total number of {tumour} patients with single doses that were filtered out is {n_single_dose}"))

  # Tag chemotherapy line
  chemo_lines <- drug_episodes %>%
    dplyr::filter(stringr::str_detect(tolower(linename), treatment)) %>%
    dplyr::group_by(patientid) %>%
    dplyr::summarise(chemo_line = min(linenumber), .groups = "drop") %>%
    dplyr::filter(chemo_line %in% c(1, 2, 3))

  # Add subsequent line
  chemo_and_next <- drug_episodes %>%
    dplyr::inner_join(chemo_lines, by = "patientid", relationship = "many-to-many") %>%
    dplyr::group_by(patientid) %>%
    dplyr::filter(linenumber %in% c(chemo_line, chemo_line + 1)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(patientid, linenumber) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(prior_lines = chemo_line - 1)

  # Load death and follow-up data
  last_contact_data <- haven::read_dta("H:\\PREDiCText\\ek\\Flatiron_data_prep\\1_Dates_E_LastContact.dta")
  mort_dat <- haven::read_dta("H:\\PREDiCText\\lingyi\\Flatiron_exploratory_round1and4\\DeathDates\\DeathDates_AllDataSets_round1and4.dta") %>%
    dplyr::arrange(parsemonthofdeath) %>%
    dplyr::group_by(patientid) %>%
    dplyr::slice(1)

  # Merge death info and coalesce date
  last_contact_data <- last_contact_data %>%
    dplyr::left_join(mort_dat, by = "patientid") %>%
    dplyr::select(patientid, Date_LastContact, Date_LastClinicalNote, Date_TxEnd_PostReport, parsemonthofdeath) %>%
    dplyr::mutate(parsemonthofdeath = lubridate::ymd(parsemonthofdeath)) %>%
    dplyr::mutate(Date_LastFollowUp = pmax(Date_LastContact, Date_LastClinicalNote, Date_TxEnd_PostReport, na.rm = TRUE)) %>%
    dplyr::mutate(coalesced_death_proxy = pmax(Date_LastContact, Date_LastClinicalNote, Date_TxEnd_PostReport, na.rm = TRUE)) %>%
    dplyr::mutate(death_date = dplyr::case_when(
      !is.na(coalesced_death_proxy) & !is.na(parsemonthofdeath) &
        format(coalesced_death_proxy, "%Y-%m") == format(parsemonthofdeath, "%Y-%m") ~ coalesced_death_proxy,
      !is.na(parsemonthofdeath) ~ parsemonthofdeath,
      TRUE ~ as.Date(NA)
    )) %>%
    dplyr::select(patientid, Date_LastFollowUp, death_date)

  drug_episodes <- chemo_and_next %>%
    dplyr::left_join(last_contact_data, by = "patientid", relationship = "many-to-many")

  # Load and process progression data
  prog_data <- haven::read_dta("H:\\PREDiCText\\ek\\Flatiron_exploratory\\_#_Round4\\Progression\\Progression_allDataSets_round1and4.dta") %>%
    dplyr::mutate(progressiondate = as.Date(progressiondate)) %>%
    dplyr::select(patientid, progressiondate)

  first_line_dates <- drug_episodes %>%
    dplyr::group_by(patientid) %>%
    dplyr::filter(linenumber == chemo_line) %>%
    dplyr::summarise(first_line_start = linestartdate)

  prog_data_filtered <- prog_data %>%
    dplyr::left_join(first_line_dates, by = "patientid", relationship = "many-to-many") %>%
    dplyr::filter(progressiondate > first_line_start) %>%
    dplyr::group_by(patientid) %>%
    dplyr::summarise(progression_date = min(progressiondate))

  drug_episodes <- drug_episodes %>%
    dplyr::left_join(prog_data_filtered, by = "patientid", relationship = "many-to-many")

  # Remove patients with lines before 2010
  patients_pre_2010 <- drug_episodes %>%
    dplyr::filter(lubridate::year(linestartdate) < 2010) %>%
    dplyr::pull(patientid)

  drug_episodes <- drug_episodes %>%
    dplyr::filter(!patientid %in% patients_pre_2010)
  n_2010 <- length(unique(patients_pre_2010))
  message(glue::glue("The total number of {tumour} patients with treatment start pre-2010 that were filtered out is {n_2010}"))

  # Add next line and gap to next line
  drug_episodes <- drug_episodes %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(
      next_linestart = dplyr::lead(linestartdate),
      days_to_next_line = as.numeric(next_linestart - lineenddate)
    ) %>%
    dplyr::ungroup()

  return(drug_episodes)
}
