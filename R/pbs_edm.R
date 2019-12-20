#' Use EDM to Perform Out-of-Sample Forecasting 
#'
#' @param data_frame A data frame with named columns (list of numeric vectors)
#' @param lags Named list of numeric vectors giving lags to use for each
#'     column. List names must match column names. (list of numeric vectors)
#' @param from_user Rows to predict from (numeric vector)
#' @param into_user Rows to predict into (numeric vector)
#' @param forecast_distance Forecast distance (numeric scalar)
#' @param symmetric_exclusion Symmetric exclusion radius? (logical scalar)
#' @param include_stats Return forecast stats? (logical scalar)
#' @param include_forecasts Return forecast and observations? (logical scalar)
#' @param include_neighbours Return neighbours? (logical scalar)
#'
#' @return A tibble (work in progress)
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#'
#' @examples
#'   data_frame <- data.frame(x = simple_ts, y = simple_ts, z = simple_ts)
#'   lags <- list(x = 0:3, y = 0:1, z = 0)
#'   pbs_edm(data_frame, lags)
#' 
pbs_edm <- function(data_frame,
                    lags,
                    from_user = seq_len(nrow(data_frame)), 
                    into_user = seq_len(nrow(data_frame)), 
                    forecast_distance = 1L,
                    symmetric_exclusion = FALSE,
                    include_stats = TRUE,
                    include_forecasts = TRUE,
                    include_neighbours = TRUE) {
  
  # Check arguments
  stopifnot(
    is.data.frame(data_frame),
    is.list(lags),
    all(is.element(names(lags), names(data_frame))),
    is.numeric(from_user),
    is.vector(from_user),
    is.numeric(into_user),
    is.vector(into_user),
    is.numeric(forecast_distance),
    is.logical(symmetric_exclusion),
    is.logical(include_stats),
    is.logical(include_forecasts),
    is.logical(include_neighbours)
  )
  
  # Make lag matrix
  lag_tibble <- make_lag_tibble(data_frame, lags)
  
  # Calculate Euclidean distances
  distance_tibble <- make_dist_tibble(lag_tibble)
  
  # Specify global indices
  global_indices <- make_global_indices(lag_tibble, 
                                        from_user, 
                                        forecast_distance)
  from_global <- dplyr::pull(global_indices, from)
  
  # Make forecasts
  forecast_tibble <- mapply(FUN = make_simplex_forecast,
                            from_index = from_global,
                            MoreArgs = list(
                              from_global = from_global,
                              lags = lags,
                              lag_tibble = lag_tibble,
                              distance_tibble = distance_tibble,
                              forecast_distance = forecast_distance,
                              symmetric_exclusion = symmetric_exclusion
                            ),
                            SIMPLIFY = FALSE) %>%
    dplyr::bind_rows() %>%
    dplyr::right_join(tibble::tibble(index = seq_len(nrow(data_frame))),
                      by = "index")
  
  # Calculate statistics
  correlation <- cor(x = dplyr::pull(forecast_tibble, observation),
                     y = dplyr::pull(forecast_tibble, forecast),
                     use = "pairwise.complete.obs")
  rmse <- sqrt(mean((dplyr::pull(forecast_tibble, observation) - 
                       dplyr::pull(forecast_tibble, forecast))^2, na.rm = TRUE))
  
  # Extract forecasts and neighbours
  forecasts <- list(dplyr::select(forecast_tibble, -neighbours))
  neighbours <- list(dplyr::select(forecast_tibble, neighbours))
  
  # Make summary tibble
  summary_tibble <- tibble::tibble(forecast_distance = forecast_distance)
  
  if (include_stats) {
    summary_tibble <- summary_tibble %>% 
      dplyr::mutate(correlation = correlation, rmse = rmse)
  }
  if (include_forecasts) {
    summary_tibble <- summary_tibble %>% dplyr::mutate(forecasts = forecasts)
  }
  if (include_neighbours) {
    summary_tibble <- summary_tibble %>% dplyr::mutate(neighbours = neighbours)
  }
  
  # Return summary tibble
  summary_tibble
}













