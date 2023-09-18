#' Lag Superset Columns
#'
#' @param data [matrix()] or [data.frame()] with named [numeric()] columns
#' @param lags [list()] of a named vector of lags for each explanatory
#'   variable.
#' @param superset [list()] superset of lags corresponding to the parent state
#'   space reconstruction
#' @param beyond [logical()]
#'
#' @return [tibble::tibble()]
#' @export
#'
superset_columns <- function (data,
                              lags,
                              superset = NULL,
                              beyond = FALSE) {
  # Define superset columns
  if (is.null(superset)) superset <- lags
  lag_sizes <- unlist(lags, use.names = FALSE)
  lag_names <- rep(names(lags), lengths(lags))
  lag_cols <- paste0(lag_names, "_", lag_sizes)
  sup_sizes <- unlist(superset, use.names = FALSE)
  sup_names <- rep(names(superset), lengths(superset))
  if (beyond) {
    n_rows <- nrow(data) + 1
  } else {
    n_rows <- nrow(data)
  }
  lag_mat <- matrix(0, nrow = n_rows, ncol = length(sup_sizes))
  colnames(lag_mat) <- paste0(sup_names, "_", sup_sizes)
  lag_mat[, lag_cols] <- 1
  # Return
  return(tibble::as_tibble(lag_mat))
}
