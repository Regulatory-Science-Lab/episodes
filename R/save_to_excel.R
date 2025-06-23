#' Write Weibull Shape and Scale Parameters to Workbook
#'
#' This function appends new Weibull scale and shape parameter estimates
#' to an existing Excel workbook. It modifies the workbook in memory
#' and returns the updated workbook object (without saving yet).
#'
#' @param wb An `openxlsx` workbook object loaded into R memory.
#' @param params_results A data frame containing parameter estimates,
#' including `state_transition`, `scale`, `shape`, `scale_lci`, `scale_uci`,
#' `shape_lci`, and `shape_uci` columns from the `fit_weibull_by_state` function
#' @param tumour A string indicating the tumour type (e.g., "lung", "breast").
#' @param treatment A string indicating the treatment name (e.g., "cisplatin").
#'
#' @return The modified workbook (`wb`) with new entries appended to the appropriate sheets.
#'
#' @details
#' - Appends `scale` parameters to the "1.7_Weibull_Scale_SoC" sheet.
#' - Appends `shape` parameters to the "1.6_Weibull_Shape_SoC" sheet.
#' - Adds a horizontal border (line) after each new block of inserted rows.
#'
#' Note: You must save the workbook separately after all modifications using `saveWorkbook()`.
#' @export
write_shape_scale <- function(wb, params_results, tumour = "lung", treatment = "cisplatin") {

  scale_dat <- params_results %>%
    dplyr::mutate(tumour = tumour,
                  treatment = treatment,
                  dist = "weibull") %>%
    dplyr::select(tumour, treatment, state = state_transition, value = scale, dist, par1 = scale_lci, par2 = scale_uci)

  shape_dat <- params_results %>%
    dplyr::mutate(tumour = tumour,
                  treatment = treatment,
                  dist = "weibull") %>%
    dplyr::select(tumour, treatment, state = state_transition, value = shape, dist, par1 = shape_lci, par2 = shape_uci)

  wb <- write_params_helper(wb, param_dat = scale_dat, sheet_name = "1.7_Weibull_Scale_SoC")

  wb <- write_params_helper(wb, param_dat = shape_dat, sheet_name = "1.6_Weibull_Shape_SoC")

  return(wb)
}


#' Helper function to append parameters to a specific workbook sheet
#'
#' This helper function finds the next empty row in a specified Excel sheet,
#' writes a block of data, adds a thin horizontal border line below the block,
#' and returns the modified workbook.
#'
#' @param wb An `openxlsx` workbook object.
#' @param param_dat A data frame containing the block of parameters to append.
#' @param sheet_name A string specifying which sheet to write to.
#'
#' @return The modified workbook (`wb`) after writing and styling.
#'
#' @details
#' - Automatically detects the next available row based on existing data.
#' - Does not overwrite previous entries.
#' - Adds a visual separator line after each data block.
#'
#' @export
#' @export
write_params_helper <- function(wb, param_dat, sheet_name = "1.6_Weibull_Shape_SoC") {

  existing_dat <- tryCatch({
    openxlsx::readWorkbook(wb, sheet = sheet_name)
  }, error = function(e) NULL)

  if (!is.null(existing_dat) && nrow(existing_dat) > 0) {
    start_row <- nrow(existing_dat) + 2
  } else {
    start_row <- 2
  }

  message(glue::glue("Writing to sheet {sheet_name} starting at row {start_row}"))

  openxlsx::writeData(wb, sheet = sheet_name, x = param_dat, startRow = start_row, colNames = FALSE)

  openxlsx::addStyle(wb,
                     sheet = sheet_name,
                     style = openxlsx::createStyle(border = "top", borderStyle = "thin"),
                     rows = start_row + nrow(param_dat),
                     cols = 1:7,
                     gridExpand = TRUE
  )

  return(wb)
}



#'  Write sample sizes in states to Workbook
#' This function appends new Weibull scale and shape parameter estimates
#' to an existing Excel workbook. It modifies the workbook in memory
#' and returns the updated workbook object (without saving yet).
#'
#' @param wb An `openxlsx` workbook object loaded into R memory.
#' @param summary_table a table of summary counts from the `state_numbers_summary` function
#' @param tumour A string indicating the tumour type (e.g., "lung", "breast").
#' @param treatment A string indicating the treatment name (e.g., "cisplatin").
#'
#' @return The modified workbook (`wb`) with new entries appended to the appropriate sheets.
#'
#' @details
#' - Appends state number sample size counts to the "1.20_State_Numbers" sheet.
#' - Adds a horizontal border (line) after each new block of inserted rows.
#'
#' Note: You must save the workbook separately after all modifications using `saveWorkbook()`.
#' @export
write_summary <- function(wb, summary_table, tumour = "lung", treatment = "cisplatin") {

  summary_table <- summary_table %>%
    dplyr::mutate(tumour = tumour, treatment = treatment) %>%
    dplyr::select(tumour, treatment, everything())


  wb <- write_params_helper(wb, param_dat = summary_table, sheet_name = "1.20_State_Numbers")

  return(wb)
}



#'  Write death table to Workbook
#' This function appends new Weibull scale and shape parameter estimates
#' to an existing Excel workbook. It modifies the workbook in memory
#' and returns the updated workbook object (without saving yet).
#'
#' @param wb An `openxlsx` workbook object loaded into R memory.
#' @param death_table a table of death and state counts by month from the `death_table` function
#' @param tumour A string indicating the tumour type (e.g., "lung", "breast").
#' @param treatment A string indicating the treatment name (e.g., "cisplatin").
#'
#' @return The modified workbook (`wb`) with new entries appended to the appropriate sheets.
#'
#' @details
#' - Appends state number sample size counts to the "1.20_State_Numbers" sheet.
#' - Adds a horizontal border (line) after each new block of inserted rows.
#'
#' Note: You must save the workbook separately after all modifications using `saveWorkbook()`.
#' @export
write_death_table <- function(wb, death_table, tumour = "lung", treatment = "cisplatin") {

  death_table <- death_table %>%
    dplyr::mutate(tumour = tumour, treatment = treatment) %>%
    dplyr::select(tumour, treatment, everything())

  wb <- write_params_helper(wb, param_dat = death_table, sheet_name = "1.19_Death_Table")

  return(wb)
}

