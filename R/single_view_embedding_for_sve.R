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
#' @param superset [list()] superset of lags corresponding to parent SSR TODO NEEDED??
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
                                            response,
                                            lags)

  # Compute state space distances between points -------------------------------

  # Matrix of allowed neighbour distances, rows corresponding to
  #   focal point times, columns to neighbour point time. The value represents the
  #   distance from the focal point to the neighbour point. Disallowed focal
  #   point and neighbour combinations have value NA.

  distances <- state_space_distances_for_sve(ssr)

  # Make forecasts for all allowable focal times $t^*$, taking into account the
  # candidate nearest neighbours

  ssr_forecasts <- state_space_forecasts_for_sve(ssr,
                                                 distances)

  # TODO ANdy got to HERE, now editing state_space_forecasts_for_sve() to deal
  # with valid focal times.
  # Define observed ------------------------------------------------------------

  observed <- c(dplyr::pull(data, response), NA)[seq_along(ssr_forecasts)]    # NA seems to get added then removed, at least in mve_understanding.Rmd example

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