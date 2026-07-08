#' Create vector of named lags for making the summary plot from MVE
#'
#' Converts a list of lags (as created by [create_subset_lags()]) into a
#'   character vector of unique lag names, optionally ordered by variable.
#'
#' @param lags [list()] where each element is a list with variable names and
#'   their associated lag vectors (e.g., `list(age11 = c(1, 2))`).
#' @param order [character()] or `NULL`. Specifies the order of variable names
#'   in the output. Variables listed here appear first in the returned vector,
#'   in the order specified. Variables not in `order` appear after, in
#'   alphabetical order. If `NULL` (default), all variables are in alphabetical
#'   order. Within each variable, lag values are always in numerical order.
#'
#' @return A character vector of unique strings in the format
#'   `{variable_name}_lag_{lag_value}`, where each lag value gets its own
#'   string (e.g., `"age11_lag_1"` and `"age11_lag_2"`). Duplicates are
#'   removed. The output is ordered according to the `order` parameter.
#'
#' @export
#'
#' @examples
#' lags_list <- list(
#'   list(age11 = c(1)),
#'   list(age11 = c(1, 2))
#' )
#' create_named_lags(lags_list)
#'
#' # With specified order
#' create_named_lags(lags_list, order = c("age11"))
create_named_lags <- function(lags, order = NULL) {
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

  # Remove duplicates
  result <- unique(result)

  # Parse strings to extract variable names and lag values
  var_lag_pairs <- strsplit(result, "_lag_")
  vars <- sapply(var_lag_pairs, function(x) x[1])
  lags_extracted <- sapply(var_lag_pairs, function(x) as.numeric(x[2]))

  # Get unique variable names
  unique_vars <- unique(vars)

  # Sort variable names according to order parameter
  if (is.null(order)) {
    sorted_vars <- sort(unique_vars)
  } else {
    # Variables in order appear first (in specified order)
    vars_in_order <- order[order %in% unique_vars]
    # Remaining variables in alphabetical order
    vars_not_in_order <- sort(unique_vars[!(unique_vars %in% order)])
    sorted_vars <- c(vars_in_order, vars_not_in_order)
  }

  # Build final result with ordered variables and numerically sorted lags
  final_result <- character()
  for (var in sorted_vars) {
    # Get all lag values for this variable
    var_lags <- lags_extracted[vars == var]
    # Sort lags numerically
    sorted_var_lags <- sort(unique(var_lags))
    # Create strings in order
    for (lag in sorted_var_lags) {
      final_result <- c(final_result, paste0(var, "_lag_", lag))
    }
  }

  final_result
}
