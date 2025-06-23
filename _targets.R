library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    "dplyr", "tidyr", "haven", "openxlsx", "glue",
    "lubridate", "purrr", "survival", "flexsurv", "episodes"
  )
)

list(

  tar_target(
    tumour_defs,
    tibble::tibble(
      tumour = c("nsclc", "breast"),
      treatment = c("cisplatin|carboplatin", "capecitabine"),
      exact = c(FALSE, TRUE)
    )
  ),

  tar_target(
    tumour_outputs,
    {
      tumour <- tumour_defs$tumour
      treatment <- tumour_defs$treatment
      exact <- tumour_defs$exact

      drug_episodes <- prep_episode_data(
        tumour = tumour,
        treatment = treatment,
        drug_separator = ",",
        overlap_threshold = 1,
        exact = exact
      )

      state_data <- construct_state_episodes(drug_episodes)

      weibull_params <- state_exit_weibull_estimates(state_data$drug_transitions, state_data$lines)
      state_summary  <- state_numbers_summary(state_data$drug_transitions)
      death_table    <- death_table(state_data$drug_transitions)

      print(death_table)

      list(
        tumour = tumour,
        treatment = treatment,
        weibull_params = weibull_params,
        state_summary = state_summary,
        death_table = death_table
      )
    },
    pattern = map(tumour_defs),
    iteration = "list"
  ),

  # Final Excel output step
  tar_target(
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
      message(glue::glue("Saved output to: {out_path}"))
      out_path
    }
  )
)
