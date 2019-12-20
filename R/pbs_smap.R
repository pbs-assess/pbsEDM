

pbs_smap <- function(data_frame,
                     lags,
                     local_weight = 0,
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
    is.numeric(local_weight),
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
  forecast_tibble <- mapply(FUN = make_smap_forecast,
                            from_index = from_global,
                            MoreArgs = list(
                              from_global = from_global,
                              lags = lags,
                              local_weight = local_weight,
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

