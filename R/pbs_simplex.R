#' Perform out-of-sample forecasting via simplex projection
#'
#' @param time_series A time series with possible NA values (numeric vector).
#' @param embed_dim The embedding dimension for the lagged-coordinate 
#'   embedding of the time series (integer scalar)
#' @param lag_size The number of time steps separating successive lags
#'   (integer scalar)
#' @param forecast_dist The number of time steps forward to forecast
#'   (integer scalar)
#' @param use_indices Time series indices to predict from (integer vector).
#'   Indices will be removed if the associated lagged coordinate vector 
#'   or its projection forward by forecast_distance (1) has one or more
#'   NA values, (2) contains one or more indices that appear in the focal
#'   vector, or (3) falls outside the time series.
#' @param predict_indices Time series indices to predict to (integer vector).
#'   Indices will be removed if the associated lagged coordinate vector 
#'   or its projection forward by forecast_distance (1) has one or more
#'   NA values, (2) contains one or more indices that appear in the focal
#'   vector, or (3) falls outside the time series.
#' @param return_value One of "stats", "predictions", or "both".
#' 
#' @return Either a tibble of statistics summarizing the forecast and
#'   forecast accuracy, a vector of predicted values, or a list of both.
#' @export
#'
#' @examples pbs_simplex(rnorm(30, 0, 1))
pbs_simplex <- function (time_series,
                         embed_dim = 2,
                         lag_size = 1,
                         forecast_dist = 1,
                         use_indices = seq_len(length(time_series)),
                         predict_indices = seq_len(length(time_series)),
                         return_value = "stats") {
  # Check arguments

  # Make a matrix of time series lags
  lag_mat <- pbs_make_lags(time_series, embed_dim, lag_size)

  # Calculate Euclidean distances among row vectors
  lag_dist <- dist(lag_mat, diag = F, upper = TRUE) %>%
    broom::tidy() %>%
    tidyr::drop_na()

  # Instantiate prediction vector
  pred_vec <- rep(NA, length(time_series))
  predict_indices <- setdiff(predict_indices, length(predict_indices))

  # Iterate over prediction set
  for (time_ind in predict_indices) {
    # Identify allowable indices
    rel_libr_ind <- setdiff(
      use_indices,
      c(
        seq_len((embed_dim - 1) * lag_size),
        (time_ind - embed_dim):(time_ind + embed_dim),
        (length(time_series) - forecast_dist + 1):(length(time_series))
      )
    )
    # Identify nearest neighbours
    rel_lag_dist <- lag_dist %>%
      dplyr::filter(item2 %in% time_ind,
                    item1 %in% rel_libr_ind) %>%
      dplyr::arrange(distance) %>%
      dplyr::mutate(n = row_number()) %>%
      dplyr::filter(n <= embed_dim + 1)
    # Are there enough points?
    if (nrow(rel_lag_dist) > embed_dim) {
      # Calculate the weights
      omega_weights <- exp(-rel_lag_dist$distance / rel_lag_dist$distance[1])
      # Identify the iterated indicies
      iterated_ind <- rel_lag_dist$item1 + forecast_dist
      # Predict the value
      pred_vec[time_ind + forecast_dist] <- sum(
        time_series[iterated_ind] * omega_weights
      )  / sum(omega_weights)
    } else {
      pred_vec[time_ind] <- NA
    }
  }
  pred_obs_cor <- cor(pred_vec,
                      time_series,
                      use = "pairwise.complete.obs")
  # Return
  list(
    embed_dim = embed_dim,
    lag_size = lag_size,
    forecast_dist = forecast_dist,
    rho = pred_obs_cor,
    pred_vec = pred_vec
  )
}
