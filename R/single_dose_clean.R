#' Clean and collapse single-dose cancer treatment lines
#'
#' This function identifies and cleans single-day treatment lines in a longitudinal
#' drug exposure dataset. It attempts to merge single-day lines with adjacent lines
#' if they share overlapping drugs, and drops lines that are isolated and non-mergeable.
#'
#' It supports threshold-based partial overlap merging and reconstructs treatment
#' episodes by recomputing line numbers after merging or dropping.
#'
#' @param ep_data data frame containing line-level treatment data. Must include
#'   columns: `patientid`, `linestartdate`, `lineenddate`, `linenumber`, `linename`.
#' @param drug_separator A character string used to split multiple drugs in `linename`.
#'   Default is ",".
#' @param overlap_threshold An integer indicating the minimum number of overlapping
#'   drugs between line names required to consider lines similar enough to merge.
#'   Default is 1.
#'
#' @return A cleaned `data.frame` with the same columns (`patientid`, `linestartdate`,
#'   `lineenddate`, `linename`, `linenumber`), but with collapsed or dropped single-day lines.
#'
#' @examples
#' \dontrun{
#' cleaned_lines <- single_dose_clean(ep_data)
#' }
#' @import dplyr
#' @importFrom stringr str_squish
#' @export
single_dose_clean <- function(ep_data, drug_separator = ",", overlap_threshold = 1) {

  # Helper: split drugs into a sorted set
  split_drugs <- function(name) sort(trimws(unlist(strsplit(name, drug_separator))))

  # Helper: check if two line names have overlapping drugs (at least `threshold`)
  has_overlap <- function(name1, name2, threshold = overlap_threshold) {
    if (is.na(name1) | is.na(name2)) return(FALSE)
    length(intersect(split_drugs(name1), split_drugs(name2))) >= threshold
  }

  # Flag single-day lines and gather adjacent line info
  ep_flagged <- ep_data %>%
    dplyr::arrange(patientid, linestartdate, linenumber) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(
      is_single_day = as.numeric(lineenddate - linestartdate) == 0,
      prev_end      = dplyr::lag(lineenddate),
      prev_name     = dplyr::lag(linename),
      prev_gap      = as.numeric(linestartdate - prev_end),
      next_start    = dplyr::lead(linestartdate),
      next_name     = dplyr::lead(linename),
      next_gap      = as.numeric(next_start - lineenddate)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      merge_with_prev = mapply(function(x, y, gap, single) {
        single & !is.na(gap) & has_overlap(x, y)
      }, linename, prev_name, prev_gap, is_single_day),

      merge_with_next = mapply(function(x, y, gap, single) {
        single & !is.na(gap) & has_overlap(x, y)
      }, linename, next_name, next_gap, is_single_day),

      drop_line = is_single_day & !merge_with_prev & !merge_with_next
    )

  # Identify which lines will be merged and to whom
  to_merge <- dplyr::bind_rows(
    ep_flagged %>% dplyr::filter(merge_with_prev) %>% dplyr::mutate(target_line = linenumber - 1),
    ep_flagged %>% dplyr::filter(merge_with_next) %>% dplyr::mutate(target_line = linenumber + 1)
  )

  # Prepare target lines to merge with
  target_lines <- ep_flagged %>%
    dplyr::select(patientid, linenumber, linestartdate, lineenddate, linename) %>%
    dplyr::rename(
      target_line   = linenumber,
      target_start  = linestartdate,
      target_end    = lineenddate,
      target_name   = linename
    )

  to_merge_full <- dplyr::left_join(to_merge, target_lines, by = c("patientid", "target_line"))

  # Collapse line info for merged lines
  collapsed_lines <- to_merge_full %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      linestartdate = min(linestartdate, target_start, na.rm = TRUE),
      lineenddate   = max(lineenddate, target_end, na.rm = TRUE),
      linename = paste(sort(unique(c(
        split_drugs(linename),
        split_drugs(target_name)
      ))), collapse = ",")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(patientid, linestartdate, lineenddate, linename) %>%
    dplyr::distinct()

  # Remove merged or dropped lines from original
  lines_to_remove <- dplyr::bind_rows(
    to_merge %>% dplyr::select(patientid, linenumber),
    to_merge %>% dplyr::select(patientid, linenumber = target_line),
    ep_flagged %>% dplyr::filter(drop_line) %>% dplyr::select(patientid, linenumber)
  )

  ep_retained <- ep_flagged %>%
    dplyr::anti_join(lines_to_remove, by = c("patientid", "linenumber")) %>%
    dplyr::select(patientid, linestartdate, lineenddate, linename)

  # Combine retained + collapsed lines and renumber
  ep_final <- dplyr::bind_rows(ep_retained, collapsed_lines) %>%
    dplyr::arrange(patientid, linestartdate) %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(linenumber = dplyr::row_number()) %>%
    dplyr::ungroup()

  return(ep_final)
}

