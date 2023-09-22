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
#' @param superset [list()] superset of lags corresponding to parent SSR TODO
#'   NEEDED?? No, it will be kept track of within sve function and from what was
#'   called, no need to repeat.
#'
#' @author Andrew M. Edwards and Luke A. Rogers
#'
#' @return [data.frame()]
#' @export
#'
single_view_embedding_for_sve <- function(data,
                                          response,
                                          lags,
                                          metric = "rmse",
                                          superset = NULL) {

  # Define the state space reconstruction, using scaled variables --------------

  ssr <- state_space_reconstruction_for_sve(data,
                                            lags)

  # Also need scaled response variable
  response_as_lag <- list(xxx = 0)
  names(response_as_lag) <- response

  response_s <- state_space_reconstruction_for_sve(data,
                                                   response = response,
                                                   lags = response_as_lag,
                                                   response_only = TRUE)


  # Compute state space distances between points -------------------------------

  # Matrix of allowed neighbour distances, rows corresponding to
  #   focal point times, columns to neighbour point time. The value represents the
  #   distance from the focal point to the neighbour point. Disallowed focal
  #   point and neighbour combinations have value NA.

  distances <- state_space_distances_for_sve(ssr)

  # Make predictions for all allowable focal times $t^*$, taking into account the
  # candidate nearest neighbours

  prediction_indices <- state_space_forecasts_for_sve(ssr,
                                                      distances,
                                                      max_lag = max(unlist(lags,
                                                                           use.names = FALSE)))


  # Take off the final one because we are shifting (predicting t_star+1 from the
  # prediction indices for each t_star; need NA at beginning.
  response_s_prediction <- c(NA,
                             response_s[prediction_indices[-length(prediction_indices)]])

  # Define observed ------------------------------------------------------------
  # observed <- c(dplyr::pull(data, response), NA)[seq_along(ssr_forecasts)]    # NA seems to get added then removed, at least in mve_understanding.Rmd example



  # HERE # see vignette - need to untransofrm

  # Compute forecast -----------------------------------------------------------

  forecast <- untransform_forecasts(observed, ssr_forecasts)

  # Return results -------------------------------------------------------------

  rows <- seq_along(forecast)

  tibble::tibble(
            set = rep(0:1, c(index - 1, nrow(data) - index + 2))[rows],
            time = seq_len(nrow(data) + 1L)[rows],
            points = c(0, as.vector(apply(distances, 1, function (x) sum(!is.na(x)))))[rows],
            dim = rep(ncol(ssr), nrow(data) + 1L)[rows],
            observed = observed,
            forecast = forecast,
            forecast_metrics(observed, forecast, window, metric),
            superset_columns(data, lags, superset)
          )
}
