#' Single View Embedding
#'
#' @param data [matrix()] or [data.frame()] with named [numeric()] columns
#' @param response [character()] column name of the response variable in
#'   \code{data}
#' @param lags [list()] of a named vector of lags for each explanatory
#'   variable.
#' @param index [integer()]
#' @param buffer [integer()] number of forecasts prior to \code{index}
#' @param window [integer()] forecast metric moving window width
#' @param metric [character()]
#' @param beyond [logical()]
#' @param superset [list()] superset of lags corresponding to parent SSR
#'
#' @author Luke A. Rogers
#'
#' @return [data.frame()]
#' @export
#'
single_view_embedding <- function (data,
                                   response,
                                   lags,
                                   index = 50L,
                                   buffer = 10L,
                                   window = integer(0),
                                   metric = "rmse",
                                   beyond = FALSE,
                                   superset = NULL) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_integerish(index, lower = 20, upper = nrow(data), len = 1)
  checkmate::assert_integerish(buffer, lower = 1, upper = 10, len = 1)

  # Define the state space reconstruction --------------------------------------

  ssr <- state_space_reconstruction(data, response, lags)

  # Compute state space distances between points -------------------------------

  # - Rows in X are points in the SSR
  # - Each row in X_distance corresponds to a focal point in the SSR
  # - Each column in X_distance corresponds to a potential neighbour in the SSR
  # - Elements of X_distance correspond to distances to neighbours
  # - NA elements indicate disallowed neighbours for a given focal point

  distances <- state_space_distances(ssr, index, buffer)

  # Compute centred and scaled forecasts ---------------------------------------

  # - Create neighbour index matrix
  # - Create neighbour matrices
  # - Project neighbour matrices
  # - Compute X_forecast vector

  ssr_forecasts <- state_space_forecasts(ssr, distances, beyond)

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
