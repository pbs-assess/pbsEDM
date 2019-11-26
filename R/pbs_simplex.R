#' Perform out-of-sample forecasting via simplex projection
#'
#' @param time_series A numeric vector with possible NA values
#' @param embed_dim An integer embedding dimension
#' @param libr_ind An integer vector of indices to predict from
#' @param pred_ind An integer vector of indicies to predict to
#' @param forecast_dist The number of time steps forward for each forecast
#' @param lag_size The number of time steps separating successive lags
#'
#' @return
#' @export
#'
#' @examples
pbs_simplex <- function (time_series,
                         embed_dim, 
                         libr_ind = seq_along(time_series), 
                         pred_ind = libr_ind,
                         lag_size = 1, 
                         forecast_dist = 1) {
  # Check arguments
  # Make lagged matrix
  lag_mat <- pbs_make_lags(time_series, embed_dim, lag_size)
  # Calculate Euclidean distances
  lag_dist <- dist(lag_mat, upper = TRUE) %>% tidy() %>% drop_na()
  # Instantiate prediction vector
  pred_vec <- rep(NA, length(time_series))
  # Iterate over prediction set
  for (time_ind in pred_ind) {
    # Identify allowable indices
    rel_libr_ind <- setdiff(
      libr_ind, 
      c(
        seq_len((embed_dim - 1) * lag_size),
        (time_ind - embed_dim):(time_ind + embed_dim),
        (length(time_series) - forecast_dist + 1):(length(time_series))
      )
    )
    # Identify nearest neighbours
    rel_lag_dist <- lag_dist %>% 
      filter(item2 %in% time_ind, item1 %in% rel_libr_ind) %>%
      arrange(distance) %>%
      mutate(n = row_number()) %>%
      filter(n <= embed_dim + 1)
    # Are there enough points?
    if (nrow(rel_lag_dist) > embed_dim) {
      # Calculate the weights
      omega_weights <- exp(-rel_lag_dist$distance/rel_lag_dist$distance[1])
      # Identify the iterated indicies
      iterated_ind <- rel_lag_dist$item1 + forecast_dist
      # Predict the value
      pred_vec[time_ind] <- sum(
        time_series[iterated_ind] * omega_weights / sum(omega_weights)
      )
    } else {
      pred_vec[time_ind] <- NA
    }
  }
  pred_obs_cor <- cor(pred_vec, time_series, use = "pairwise.complete.obs")
  # Return
  list(
    embed_dim = embed_dim,
    lag_size = lag_size,
    forecast_dist = forecast_dist,
    rho = pred_obs_cor,
    pred_vec = pred_vec
  )
}



