#' Untransform State Space Forecast Values
#'
#' @param x [numeric()][vector()] observed values
#' @param y [numeric()][vector()] state space forecast values
#'
#' @return [numeric()][vector()]
#' @export
#'
untransform_forecasts <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return untransformed forecast
  return(mean(x, na.rm = TRUE) + y * stats::sd(x, na.rm = TRUE))
}
