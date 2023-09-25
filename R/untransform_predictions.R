#' Untransform Vector of State Space Predictions
#'
#' So this is the reverse of `state_space_reconstruction_for_sve()`, but we only
#'   need to do it for the response variable, so just a single vector. Using the
#'   notation of in Appendix of first manuscript.
#'
#' @param N_observed [numeric()][vector()] original observed values
#' @param Y_predicted [numeric()][vector()] state space predictions (from the
#'   state space that used lagged and scaled variables), though just the
#'   response column
#' @param max_lag maximum lag used
#' @return [numeric()][vector()]
#' @export
#'
untransform_predictions <- function (N_observed,
                                   Y_predicted,
                                   max_lag) {
  # Check arguments
  checkmate::assert_true(length(N_observed) == length(Y_predicted))

  # First Y_predicted that is not NA (see notes)
  #  should be Y_predicted[max_lag + 2], and we have up to
  #  Y_predicted[T]. This will likely change when don't do first-differencing.

  min_t_Y_predicted <- max_lag + 2

  stopifnot(all(is.na(Y_predicted[1:(min_t_Y_predicted - 1)])))
  stopifnot(!any(is.na(Y_predicted[(min_t_Y_predicted):length(Y_predicted)])))

  Z_observed <- c(diff(N_observed), NA)

  Z_mean <- mean(Z_observed,
                na.rm = TRUE)

  Z_sd <- sd(Z_observed,
             na.rm = TRUE)

  Z_predicted <- Z_mean + Z_sd * Y_predicted

  # Then something like (double check the indexing)
  N_predicted <- rep(NA,
                     length(N_observed) + 1)
  N_predicted[min_t_Y_predicted] <- N_observed[min_t_Y_predicted]  # Only for
                                        # first-differencing, have to set first one
  for(i in (min_t_Y_predicted):length(N_observed)){
     N_predicted[i + 1] <- Z_mean + Z_sd * Y_predicted[i] + N_predicted[i]
  }
  N_predicted[min_t_Y_predicted] <- NA   # For first-differencing, haven't
                                        # predictd this one so set to NA

  # Return a prediction vector of original N_t, but also with N_{T+1},
  #  which is what we are ultimately after. So in original absolute space
  return(N_predicted)
}
