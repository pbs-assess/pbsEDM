#' Untransform Vector of State Space Predictions
#'
#' So this is the reverse of `state_space_reconstruction_for_sve()`, but we only
#'   need to do it for the response variable, so just a single vector.
#'
#' @param original [numeric()][vector()] original observed values
#' @param predictions_s [numeric()][vector()] state space predictions (from the
#'   state space that used lagged and scaled variables
#'
#' @return [numeric()][vector()]
#' @export
#'
untransform_forecasts <- function (original,
                                   predictions_s) {
  # Check arguments
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))

  # Return predictions in original absolute space
  predictions_absolute <- mean(7)
  # return(mean(x, na.rm = TRUE) + y * stats::sd(x, na.rm = TRUE))
}
