library(targets)
library(tarchetypes)

utils::globalVariables(c(
  "Date_LastClinicalNote", "Date_LastContact", "Date_LastFollowUp",
  "Date_TxEnd_PostReport", "after_death", "chemo_line", "days_to_death",
  "days_to_progression", "death_date", "death_flag", "death_month",
  "detaileddrugcategory", "dist", "drop_line", "episodedate", "event",
  "first_line_start", "gender", "group_id", "is_single_day", "label",
  "lineenddate", "lineenddate_off", "linename", "linenumber",
  "linesetting", "linestartdate", "merge_with_next", "merge_with_prev",
  "month_relative", "new_line", "next_gap", "next_linestart", "next_name",
  "next_start", "next_state", "on_treatment", "overlaps_previous",
  "parsemonthofdeath", "patientid", "prev_end", "prev_gap", "prev_name",
  "prior_lines", "progression_date", "progression_state",
  "progressiondate", "scale_lci", "scale_uci", "seg_end", "shape",
  "shape_lci", "shape_uci", "state", "state_transition", "target_end",
  "target_line", "target_name", "target_start", "time", "."))

targets::tar_option_set(
  packages = c(
    "dplyr", "tidyr", "haven", "openxlsx", "glue",
    "lubridate", "purrr", "survival", "flexsurv", "episodes"
  ),
  imports = "episodes"
)

list(

  # Step 0: Track the .rds file itself
  tar_target(
    tumour_path,
    Sys.getenv("TUMOUR_PATH", unset = NA),
    format = "file"
  ),

  # Step 1: Load tumour_defs from the tracked file path
  targets::tar_target(
    tumour_defs,
    {
      if (is.na(tumour_path) || !file.exists(tumour_path)) {
        stop("TUMOUR_PATH env var not set or file does not exist.")
      }
      tumour_data <- readRDS(tumour_path)
      stopifnot(is.data.frame(tumour_data))
      tumour_data
    }
  ),

  # Step 2: Dynamic branching over rows
  targets::tar_target(
    tumour_row,
    tumour_defs,
    pattern = map(tumour_defs),
    iteration = "list"
  ),

  # Step 3: Process each row
  targets::tar_target(
    tumour_outputs,
    {
      drug_episodes <- prep_episode_data(
        tumour = tumour_row$tumour,
        treatment = tumour_row$treatment,
        drug_separator = ",",
        overlap_threshold = 1,
        exact = tumour_row$exact
      )

      state_data <- construct_state_episodes(drug_episodes)

      weibull_params <- state_exit_weibull_estimates(
        state_data$drug_transitions,
        state_data$lines
      )

      state_summary  <- state_numbers_summary(state_data$drug_transitions)
      death_table    <- death_table(state_data$drug_transitions)

      list(
        tumour = tumour_row$tumour,
        treatment = tumour_row$treatment,
        weibull_params = weibull_params,
        state_summary = state_summary,
        death_table = death_table
      )
    },
    pattern = map(tumour_row),
    iteration = "list"
  ),

  # Step 4: Final output
  targets::tar_target(
    combined_excel_output,
    {
      wb <- openxlsx::loadWorkbook("H:/PREDiCText/nirupama/weibull_estimates/01_public_parameters_updatedNT.xlsx")

      for (output in tumour_outputs) {
        wb <- write_shape_scale(wb, output$weibull_params, output$tumour, output$treatment)
        wb <- write_summary(wb, output$state_summary, output$tumour, output$treatment)
        wb <- write_death_table(wb, output$death_table, output$tumour, output$treatment)
      }

      timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
      out_path <- glue::glue("H:/PREDiCText/nirupama/weibull_estimates/01_public_parameters_updatedNTGC_{timestamp}.xlsx")
      openxlsx::saveWorkbook(wb, out_path, overwrite = FALSE)
      out_path
    },
    cue = targets::tar_cue(mode = "always")
  )
)
