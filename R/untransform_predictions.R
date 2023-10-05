#' Untransform Vector of State Space Predictions
#'
#' So this is the reverse of `state_space_reconstruction_for_sve()`, but we only
#'   need to do it for the response variable, so just a single vector. Was using the
#'   notation of in Appendix of first manuscript, but only one vector so
#'   changing to `response_` to match the call from `single_view_embedding_for_sve()`.
#'
#' @param response_observed [numeric()][vector()] original observed values
#' @param response_s_predicted [numeric()][vector()] state space predictions (from the
#'   state space that used lagged and scaled variables), though just the
#'   response column
#' @param max_lag maximum lag used
#' @param positive_response_only logical, if TRUE then set any negative
#'   predictions of the response variable to be the minimum observed value
#' @return [numeric()][vector()] a prediction vector of original `N_t` in
#'   absolute space, so
#'   `response_predicted`, but also with a value for `T+1`, which is what we are
#'   ultimately after.

#' @export
#'
untransform_predictions <- function(response_observed,
                                    response_s_predicted,
                                    max_lag,
                                    positive_response_only = TRUE){
  # Check arguments
  checkmate::assert_true(length(response_observed) == length(response_s_predicted))

  # First response_s_predicted that is not NA (see notes)
  #  should be response_s_predicted[max_lag + 2], and we have up to
  #  response_s_predicted[T]. This will likely change when don't do first-differencing.

  min_t_response_s_predicted <- max_lag + 2

  stopifnot(all(is.na(response_s_predicted[1:(min_t_response_s_predicted - 1)])))
  stopifnot(!any(is.na(response_s_predicted[(min_t_response_s_predicted):length(response_s_predicted)])))    # These detected some errors, so keep in.

  Z_observed <- c(diff(response_observed), NA)

  Z_mean <- mean(Z_observed,
                na.rm = TRUE)

  Z_sd <- sd(Z_observed,
             na.rm = TRUE)

  Z_predicted <- Z_mean + Z_sd * response_s_predicted

  # response_predicted_{t+1} = response_observed_{t} + Z_predicted_{t}
  #  I'd originaly used response_predicted_{t} which compounded any errors,
  #  leading to best rho of 0.59 for simualted_4 and lags_lots in
  #  mve_understanding.Rmd.
  #  Think the NA's should just flow through.
  response_predicted <- c(NA,                     # deals with the t+1
                          response_observed + Z_predicted)

  if(positive_response_only){
    min_response_observed <- min(response_observed)
    response_predicted[response_predicted < 0] = min_response_observed
  }

  return(response_predicted)
}
