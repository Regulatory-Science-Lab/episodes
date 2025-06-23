#################################### Helper functions for drug episodes construction ###################################
#' Helper Functions and Main Constructor for State Episodes
#'
#' These functions help create cleaned and structured multi-state transition datasets
#' from real-world treatment episode data, suitable for survival analysis or cost-effectiveness modeling.
#'
#' @param drug_episodes A data.frame of drug episode-level information with survival events.
#' @param prog_on_cutoff Integer. Maximum negative offset allowed between progression and lineenddate.
#' @param death_cutoff,progression_cutoff,next_line_cutoff Numeric thresholds used to define transitions.
#'
#' @return A data.frame of structured drug transitions.
#' @importFrom dplyr filter mutate select transmute arrange group_by ungroup lead pull case_when bind_rows if_else summarise count
#' @importFrom lubridate year month
#' @importFrom knitr kable
#' @export
construct_state_episodes <- function(drug_episodes) {
  cutoffs <- generate_cutoffs(drug_episodes)
  drug_episodes <- adjust_dates(drug_episodes)
  lines <- initialize_lines(drug_episodes)
  lines <- fix_death_before_progression(lines)

  on_treatment_list <- build_on_treatment(
    lines,
    death_cutoff = cutoffs$death_cutoff,
    progression_cutoff = cutoffs$progression_cutoff,
    next_line_cutoff = cutoffs$next_line_cutoff
  )

  on_treatment <- on_treatment_list$on_treatment
  progression <- build_progression_from_on(on_treatment, lines)
  off_treatment <- build_off_from_on(on_treatment, lines)

  drug_transitions <- dplyr::bind_rows(
    on_treatment_list$on_treatment_clean,
    progression,
    off_treatment
  )

  drug_transitions <- fix_line_gaps(drug_transitions)

  drug_transitions <- drug_transitions %>%
    dplyr::arrange(patientid, linestartdate, lineenddate) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(event = dplyr::if_else(dplyr::row_number() == dplyr::n(), 0L, 1L)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(lineenddate >= linestartdate)

  transition_counts <- drug_transitions %>%
    dplyr::arrange(patientid, linestartdate, lineenddate) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(next_state = dplyr::lead(state)) %>%
    dplyr::ungroup() %>%
    dplyr::count(state, next_state, sort = TRUE) %>%
    dplyr::arrange(state, next_state)

  print(knitr::kable(transition_counts, caption = "Filtered Transition Counts Table"))

  return(list(drug_transitions = drug_transitions, lines = lines))
}

#' Adjust Dates in Drug Episode Data
adjust_dates <- function(drug_episodes) {
  drug_episodes %>%
    dplyr::mutate(
      progression_date = dplyr::if_else(
        as.numeric(progression_date - linestartdate) < 7 &
          as.numeric(lineenddate - progression_date) > 14,
        as.Date(NA),
        progression_date
      ),
      death_date = dplyr::if_else(
        !is.na(death_date) &
          !is.na(lineenddate) &
          death_date < lineenddate &
          lubridate::year(death_date) == lubridate::year(lineenddate) &
          lubridate::month(death_date) == lubridate::month(lineenddate),
        lineenddate,
        death_date
      )
    )
}

#' Initialize Line-Level Features
initialize_lines <- function(drug_episodes) {
  drug_episodes %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(
      next_start = dplyr::lead(linestartdate),
      base_time = dplyr::first(linestartdate),
      chemo_line = chemo_line
    ) %>%
    dplyr::ungroup()
}

#' Fix Deaths Before Progression
fix_death_before_progression <- function(lines) {
  lines %>%
    dplyr::mutate(
      death_date = dplyr::if_else(
        !is.na(death_date) & !is.na(progression_date) & death_date < progression_date,
        pmax(progression_date, Date_LastFollowUp, na.rm = TRUE),
        death_date
      )
    )
}

#' Build On-Treatment State
build_on_treatment <- function(lines, death_cutoff, progression_cutoff, next_line_cutoff) {
  on_treatment <- lines %>%
    dplyr::filter(linenumber == chemo_line) %>%
    dplyr::mutate(
      state = "On_Treatment_Target_Line",
      next_state = dplyr::case_when(
        !is.na(death_date) & is.na(progression_date) & is.na(next_linestart) &
          difftime(death_date, lineenddate) <= death_cutoff ~ "Death",
        !is.na(progression_date) &
          difftime(progression_date, lineenddate, units = "days") <= progression_cutoff ~ "progression",
        !is.na(next_start) &
          difftime(next_start, lineenddate, units = "days") <= next_line_cutoff ~ "progression",
        TRUE ~ "off_treatment"
      ),
      lineenddate = dplyr::case_when(
        next_state == "progression" & progression_date < next_linestart ~ progression_date - 1,
        next_state == "progression" & !is.na(next_linestart) &
          (next_linestart < progression_date | is.na(progression_date)) ~ next_start - 1,
        next_state == "Death" ~ death_date - 1,
        TRUE ~ lineenddate
      )
    )

  death_patients <- on_treatment %>%
    dplyr::filter(next_state == "Death") %>%
    dplyr::pull(patientid)

  death_target_line <- lines %>%
    dplyr::filter(patientid %in% death_patients) %>%
    dplyr::filter(linenumber == chemo_line) %>%
    dplyr::transmute(
      patientid,
      linestartdate = death_date,
      lineenddate = death_date,
      state = "Death", chemo_line)


  on_treatment_clean <- on_treatment %>%
    dplyr::select(patientid, linestartdate, lineenddate, state, chemo_line) %>%
    dplyr::bind_rows(death_target_line)
  return(list(on_treatment = on_treatment, on_treatment_clean = on_treatment_clean))
}


#' Build progression from on-treatment state
build_progression_from_on <- function(on_treatment, lines) {
patients_progress <- on_treatment %>%
  dplyr::filter(next_state == "progression") %>%
  dplyr::pull(patientid)

progression <- lines %>%
  dplyr::filter(patientid %in% patients_progress) %>%
  dplyr::filter(linenumber == chemo_line) %>%
  dplyr::group_by(patientid) %>%
  dplyr::mutate(
    seg_end = case_when(
      !is.na(death_date) & is.na(next_start)  ~ death_date -1,
      !is.na(next_linestart) ~ next_start - 1,
      TRUE ~ Date_LastFollowUp
    ),

    progression_state = dplyr::case_when(
      !is.na(death_date) & is.na(next_start)  ~ "Death",
      !is.na(next_linestart) ~ "next_line",
      TRUE ~ "censored"
    )
  ) %>%
  dplyr::filter(!is.na(seg_end)) %>%
  dplyr::transmute(
    patientid,
    linestartdate = pmin(progression_date, next_linestart, na.rm = TRUE),
    lineenddate = seg_end,
    state = "progression",
    progression_state,
    chemo_line
  )

progressed_patients <- progression %>%
  dplyr::filter(progression_state == "next_line") %>%
  dplyr::pull(patientid)

progression_next_line <- lines %>%
  dplyr::filter(linenumber != chemo_line) %>%
  dplyr::filter(patientid %in% progressed_patients) %>%
  dplyr::transmute(
    patientid,
    linestartdate,
    lineenddate = case_when(
      !is.na(death_date) ~ death_date - 1,
      TRUE ~ Date_LastFollowUp
    ),
    progression_state = dplyr::case_when(
      !is.na(death_date) ~ "Death",
      TRUE ~ "censored"
    ),
    chemo_line, state = "progression",
  )


progression <- progression %>%
  dplyr::filter(progression_state != "next_line")

progression <- progression %>%
  dplyr::bind_rows(progression_next_line)

progress_death <- progression %>%
  dplyr::filter(progression_state == "Death") %>%
  dplyr::pull(patientid)


death_progress <- lines %>%
  dplyr::filter(patientid %in% progress_death) %>%
  dplyr::filter(linenumber == chemo_line) %>%
  dplyr::transmute(
    patientid,
    linestartdate = death_date,
    lineenddate = death_date,
    state = "Death", chemo_line)

progression <- progression %>%
  dplyr::select(-progression_state)

progression <- dplyr::bind_rows(progression, death_progress)

return(progression)
}

#' Build off-treatment from on-treatment
build_off_from_on <- function(on_treatment, lines) {

  patients_off <- on_treatment %>%
    dplyr::filter(next_state == "off_treatment") %>%
    dplyr::pull(patientid)

  off_treatment <- lines %>%
    dplyr::filter(patientid %in% patients_off) %>%
    dplyr::filter(linenumber == chemo_line) %>%
    dplyr::mutate(
      linestartdate = lineenddate + 1,

      next_state = dplyr::case_when(
        !is.na(next_linestart) & progression_date < next_linestart | !is.na(progression_date) & is.na(next_linestart) ~ "progression",
        !is.na(death_date) & is.na(progression_date) & is.na(next_linestart) ~ "Death",
        !is.na(progression_date) & progression_date > next_linestart | is.na(progression_date) & !is.na(next_linestart) ~ "progression",
        TRUE ~ "censored"
      ),

      lineenddate_off = dplyr::case_when(
        next_state == "progression" & !is.na(progression_date) & progression_date > next_linestart | is.na(progression_date) & !is.na(next_linestart) ~ next_linestart - 1,
        next_state == "progression" & !is.na(next_linestart) & progression_date < next_linestart | !is.na(progression_date) & is.na(next_linestart)  ~ progression_date-1,
        next_state == "Death" ~ death_date - 1,
        TRUE ~ Date_LastFollowUp
      )
    ) %>%
    dplyr::transmute(
      patientid,
      linestartdate,
      lineenddate = lineenddate_off,
      state = "off_treatment",
      next_state,
      chemo_line
    )

  off_treatment <- off_treatment %>%
    dplyr::filter(lineenddate > linestartdate)

  off_treatment_progressed <- off_treatment %>%
    dplyr::filter(next_state == "progression") %>%
    dplyr::pull(patientid)

  progression_off <- lines %>%
    dplyr::filter(patientid %in% off_treatment_progressed) %>%
    dplyr::filter(linenumber == chemo_line) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(
      seg_end = dplyr::case_when(
        !is.na(death_date) & is.na(next_start)  ~ death_date -1,
        !is.na(next_linestart) ~ next_start - 1,
        TRUE ~ pmax(Date_LastFollowUp, lineenddate)
      ),

      progression_state = dplyr::case_when(
        !is.na(death_date) & is.na(next_start)  ~ "Death",
        !is.na(next_linestart) ~ "next_line",
        TRUE ~ "censored"
      )
    ) %>%
    dplyr::filter(!is.na(seg_end)) %>%
    dplyr::transmute(
      patientid,
      linestartdate = pmax(progression_date, next_linestart, na.rm = TRUE),
      lineenddate = seg_end,
      state = "progression",
      progression_state,
      chemo_line
    )

  progressed_patients <- progression_off %>%
    dplyr::filter(progression_state == "next_line") %>%
    dplyr::pull(patientid)


  progression_next_line_off <- lines %>%
    dplyr::filter(linenumber != chemo_line) %>%
    dplyr::filter(patientid %in% progressed_patients) %>%
    dplyr::transmute(
      patientid,
      linestartdate,
      lineenddate = case_when(
        !is.na(death_date) ~ death_date - 1,
        TRUE ~ pmax(Date_LastFollowUp, lineenddate)
      ),
      progression_state = dplyr::case_when(
        !is.na(death_date) ~ "Death",
        TRUE ~ "censored"
      ),
      chemo_line, state = "progression",
    )


  progression_off <- progression_off %>%
    dplyr::filter(progression_state != "next_line")

  progression_off <- progression_off %>%
    dplyr::bind_rows(progression_next_line_off)

  progress_death <- progression_off %>%
    dplyr::filter(progression_state == "Death") %>%
    dplyr::pull(patientid)


  death_progress_off <- lines %>%
    dplyr::filter(patientid %in% progress_death) %>%
    dplyr::filter(linenumber == chemo_line) %>%
    dplyr::transmute(
      patientid,
      linestartdate = death_date,
      lineenddate = death_date,
      state = "Death", chemo_line)

  off_treatment_death <- off_treatment %>%
    dplyr::filter(next_state == "Death") %>%
    dplyr::select(-next_state)

  death_rows <- off_treatment_death %>%
    dplyr::left_join(lines %>% select(patientid, death_date), by = "patientid") %>%
    dplyr::filter(!is.na(death_date)) %>%
    dplyr::mutate(
      linestartdate = death_date,
      lineenddate = death_date,
      state = "Death"
    ) %>%
    dplyr::select(patientid, linestartdate, lineenddate, state, chemo_line)

  off_treatment <- dplyr::bind_rows(off_treatment, progression_off, death_progress_off, death_rows)

  off_treatment <- off_treatment %>%
    dplyr::select(-next_state, -progression_state)

  return(off_treatment)
}

#' Make sure there are no gaps between line dates
fix_line_gaps <- function(drug_transitions) {
  drug_transitions <- drug_transitions %>%
    dplyr::arrange(patientid, linestartdate) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(
      next_linestart = dplyr::lead(linestartdate),
      lineenddate = dplyr::if_else(
        !is.na(next_linestart) & lineenddate != (next_linestart - 1),
        next_linestart - 1,
        lineenddate
      )
    ) %>%
    dplyr::select(-next_linestart) %>%
    dplyr::ungroup()
  return(drug_transitions)
}


