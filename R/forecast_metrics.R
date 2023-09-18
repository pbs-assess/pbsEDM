#' Forecast Metrics With Rolling Windows (we don't need for salmon)
#'
#' @param x [numeric()] [vector()] observed values
#' @param y [numeric()] [vector()] forecast values
#' @param k [numeric()] scalar size of the trailing running window. Uses all
#'   previous data when \code{k = integer(0)}.  TODO we don't want this.
#' @param metric [character()]
#'
#' @importFrom rlang :=
#'
#' @return [tibble::tibble()]
#' @export
#'
forecast_metrics <- function (x, y, k = integer(0), metric = "mamse") {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Observed-forecast matrix
  m <- matrix(c(x, y), ncol = 2)
  # Return metrics
  tibble::tibble(
            mre = runner::runner(m, f = matric, k = k, fun = mre),
            !!metric := runner::runner(m, f = matric, k = k, fun = get(metric))
          )
}
