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
#' @return [numeric()][vector()] a prediction vector of original `N_t` in
#'   absolute space, so
#'   `response_predicted`, but also with a value for `T+1`, which is what we are
#'   ultimately after.

#' @export
#'
untransform_predictions <- function(response_observed,
                                    response_s_predicted,
                                    max_lag){
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

  # Useful, but not needed as repeated below in loop
  # Z_predicted <- Z_mean + Z_sd * response_s_predicted

# HERE TODO pretty sure need to double check this; but seems okay. Thought wrong
# because not great answers when only doing simple view, so see if needs
# changing once multiview working.
  # Then something like (double check the indexing)
  response_predicted <- rep(NA,
                     length(response_observed) + 1)
  response_predicted[min_t_response_s_predicted] <- response_observed[min_t_response_s_predicted]  # Only for
                                        # first-differencing, have to set first one
  for(i in (min_t_response_s_predicted):length(response_observed)){
     response_predicted[i + 1] <- Z_mean + Z_sd * response_s_predicted[i] + response_predicted[i]
  }
  response_predicted[min_t_response_s_predicted] <- NA   # For first-differencing, haven't
                                        # predicted this one so set to NA

  return(response_predicted)
}
