#' Collapse overlapping treatment lines for each patient
#'
#' Identifies and merges treatment lines that overlap in time for the same patient.
#' Overlapping lines are grouped and merged into one with:
#' - Earliest `linestartdate`
#' - Latest `lineenddate`
#' - Combined `linename` (de-duplicated and cleaned)
#'
#' @param ep_data A data frame with columns `patientid`, `linestartdate`, `lineenddate`, and `linename`
#'
#' @return A cleaned data frame with one row per merged line and columns:
#'   `patientid`, `linestartdate`, `lineenddate`, `linename`, `linenumber`
#'
#' @importFrom dplyr arrange group_by mutate ungroup summarise row_number
#' @importFrom stringr str_replace_all str_squish
#' @export
#' @examples
#' \dontrun{
#'  collapse_overlapping_lines(ep_data)
#' }


collapse_overlapping_lines <- function(ep_data) {
 ep_overlap_flagged <- ep_data %>%
    dplyr::arrange(patientid, linestartdate, linenumber) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(
      prev_end = dplyr::lag(lineenddate),
      overlaps_previous = dplyr::if_else(!is.na(prev_end) & linestartdate < prev_end, TRUE, FALSE)
    ) %>%
    dplyr::ungroup()
  
  ep_grouped <- ep_overlap_flagged %>%
    dplyr::arrange(patientid, linestartdate) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(
      group_id = cumsum(!overlaps_previous | is.na(overlaps_previous))
    ) %>%
    dplyr::ungroup()
  
  ep_cleaned <- ep_grouped %>%
    dplyr::group_by(patientid, group_id) %>%
    dplyr::summarise(
      linestartdate = min(linestartdate),
      lineenddate   = max(lineenddate),
      linename      = paste(sort(unique(linename)), collapse = ","),
      .groups = "drop"
    ) %>%
    dplyr::arrange(patientid, linestartdate) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(linenumber = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # Clean up repeated drug names and standardize formatting
  ep_cleaned <- ep_cleaned %>%
    dplyr::mutate(
      linename = stringr::str_replace_all(linename, "\\b(\\w+)( \\, \\1)+\\b", "\\1"),
      linename = stringr::str_squish(gsub("\\s*\\,\\s*", ", ", linename))
    )
  
  return(ep_cleaned)
}

