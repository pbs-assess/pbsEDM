#' Create vector of named lags for making the summary plot from MVE
#'
#' Converts a list of lags (as created by [create_subset_lags()]) into a
#'   character vector of unique lag names.
#'
#' @param lags [list()] where each element is a list with variable names and
#'   their associated lag vectors (e.g., `list(age11 = c(1, 2))`).
#'
#' @return A character vector of unique strings in the format
#'   `{variable_name}_lag_{lag_value}`, where each lag value gets its own
#'   string (e.g., `"age11_lag_1"` and `"age11_lag_2"`). Duplicates are
#'   removed.
#'
#' @export
#'
#' @examples
#' lags_list <- list(
#'   list(age11 = c(1)),
#'   list(age11 = c(1, 2))
#' )
#' create_named_lags(lags_list)
create_named_lags <- function(lags) {
  result <- character()

  for (element in lags) {
    for (var_name in names(element)) {
      lag_values <- element[[var_name]]
      for (lag_value in lag_values) {
        named_lag <- paste0(var_name, "_lag_", lag_value)
        result <- c(result, named_lag)
      }
    }
  }

  # Return unique values
  unique(result)
}
