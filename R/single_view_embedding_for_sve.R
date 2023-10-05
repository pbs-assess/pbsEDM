#' Single View Embedding (updated version)
#'
#' Andy taken Luke's original and adapting it. Cumbersome function name for now,
#'   but being consistent with other new function names. Call various functions
#'   to do a single view embedding for a given set of lags as per Hao and Sugihara.
#'
#' @param data [matrix()] or [data.frame()] with named [numeric()] columns
#' @param response [character()] column name of the response variable in
#'   \code{data}
#' @param lags [list()] of a named vector of lags for each explanatory
#'   variable.
#' @param metric [character()]
#' @param ... currently just `max_allowed_correlation` to pass onto `state_space_reconstruction_for_sve()`
#'
#' @author Andrew M. Edwards and Luke A. Rogers
#'
#' @return [tibble()] with columns `response`, `response_predicted` (so both in
#'  original absolute numbers, and then `response_s` and `response_s_predicted`
#'  both in scaled (first differenced and then scaled) co-ordinates, so we can then
#'  calculate both types of correlation coefficient (or whatever else) to rank each
#'  individual single view embedding. OR a single NA if the lagged variables are
#'   highly correlated (see `state_space_reconstruction_for_sve()`.
#' @examples
#' \dontrun{
#' h_simulated <- 0.1095 + sample(1:180) * 0.001 # has mean of 0.2
#' simulated_4 <- EDMsimulate::salmon_sim(h = h_simulated)
#' res <- single_view_embedding_for_sve(data = simulated_4,
#'                              response = "R_t",
#'                              lags = create_subset_lags(list(R_t = 0:4,
#'                                                             S_t = 0:8
#'                                                             ))[[16000]]) # picking a specific subset of
#'                                                                          # potential lags
#' res %>% as.data.frame()
#' # Shows that can have R_t_predicted bigger than any original R_t, e.g. line
#' #  73, because R_t_s[72] was the largest possible, and previous R_t was not
#' #  very small. This may change with different seeds, as hadn't set, but idea
#' #  should hold.
#' }
#'
#' @export
single_view_embedding_for_sve <- function(data,
                                          response,
                                          lags,
                                          ...){

  # Define the state space reconstruction, using scaled variables --------------

  ssr <- state_space_reconstruction_for_sve(data,
                                            lags,
                                            ...)

  if(all(is.na(ssr))){      # Because lagged variables are highly correlated, so
                       # don't want to consider those reconstructions
    return(NA)         # Use this in multiview_embedding()
  }

  # Also need scaled response variable
  response_as_lag <- list(xxx = 0)
  names(response_as_lag) <- response

  response_s <- state_space_reconstruction_for_sve(data,
                                                   lags = response_as_lag,
                                                   response_only = TRUE)

  # Compute state space distances between points -------------------------------

  # Matrix of allowed neighbour distances, rows corresponding to
  #   focal point times, columns to neighbour point time. The value represents the
  #   distance from the focal point to the neighbour point. Disallowed focal
  #   point and neighbour combinations have value NA, but more are added in
  #   state_space_forecasts_for_sve() where each $t^*$ is considered.

  distances <- state_space_distances_for_sve(ssr)

  # Make predictions for all allowable focal times $t^*$, taking into account the
  # candidate nearest neighbours,

  # Define lags of the response variable
  if(response %in% names(lags)){
    lags_of_response_variable <- lags[[which(names(lags) == response)]]
  } else {
    lags_of_response_variable <- NULL
  }

# browser()

  prediction_indices <- state_space_forecasts_for_sve(ssr,
                                                      distances,
                                                      lags_of_response_variable =
                                                      lags_of_response_variable,
                                                      response_s = response_s)

  # Take off the final one because we are shifting (predicting t_star+1 from the
  # prediction indices for each t_star); need NA at beginning.
  response_s_predicted <- c(NA,
                             response_s[prediction_indices[-length(prediction_indices)]])

  # Define observed ------------------------------------------------------------
  response_observed <- dplyr::pull(data,
                                   response)

# TODO THIS SHOULD MAYBE BE LAGS_OF_RESPONSE_VARIABLE
  # Back transform and unlag to give predictions for original response variable
  response_abs_predicted <- untransform_predictions(response_observed = response_observed,
                                                    response_s_predicted = response_s_predicted,
                                                    max_lag = max(unlist(lags,
                                                      use.names = FALSE)))
  # What to return. Lags is an input so don't need that, can just have a tibble
  #  of things that are calculated within this function. Can calculate rho etc outside.
  to_return <- tibble::tibble(
                         response_obs = c(response_observed, NA),     # Since
                                        # predicted will have T+1 prediction
                         response_pred = response_abs_predicted,
                         response_scaled = c(response_s, NA),
                         response_scaled_predicted = c(response_s_predicted, NA))

  names(to_return) <- c(response,
                        paste0(response, "_predicted"),
                        paste0(response, "_s"),
                        paste0(response, "_s_predicted"))

  to_return        # Or NA earlier if highly correlated dimensions for this
                   # state space reconstruction
}
