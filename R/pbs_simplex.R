#' Use Simplex Projection to Perform Out-of-Sample Forecasting
#'
#' @param data A time series with possible NA values (numeric vector or
#'   tibble with a named column).
#' @param col_name If `data` is a tibble, the time series column name 
#'   (character scalar).
#' @param embed_dim The embedding dimension for the lagged-coordinate 
#'   embedding of the time series (integer scalar)
#' @param lag_size The number of time steps separating successive lags
#'   (integer scalar)
#' @param pred_dist The number of time steps forward to forecast
#'   (integer scalar)
#' @param lib_ind Time series indices to predict from (integer vector).
#'   Indices will be removed if the associated lagged coordinate vector 
#'   or its projection forward by pred_dist (1) has one or more
#'   NA values, (2) contains one or more indices that appear in the focal
#'   vector, or (3) falls outside the time series.
#' @param pred_ind Time series indices to predict to (integer vector).
#'   Indices will be removed if the associated lagged coordinate vector 
#'   or its projection forward by pred_dist (1) has one or more
#'   NA values, (2) contains one or more indices that appear in the focal
#'   vector, or (3) falls outside the time series.
#' 
#' @return A list of tibbles
#'   
#' @importFrom magrittr %>%
#' @export
#'
#' @examples pbs_simplex(tibble::tibble(x = 1:100), "x")
#' 
pbs_simplex <- function (data,
                         col_name = NULL,
                         embed_dim = 2L,
                         lag_size = 1L,
                         pred_dist = 1L,
                         lib_ind = seq_len(NROW(data)),
                         pred_ind = seq_len(NROW(data))) {
  # Check arguments
  stopifnot(
    is.vector(data) | tibble::is_tibble(data),
    is.numeric(data) | is.character(col_name),
    is.numeric(embed_dim),
    is.numeric(lag_size),
    is.numeric(pred_dist),
    is.numeric(lib_ind),
    is.numeric(pred_ind),
    round(embed_dim) == embed_dim,
    round(lag_size) == lag_size,
    round(pred_dist) == pred_dist,
    round(lib_ind) == lib_ind,
    round(pred_ind) == pred_ind,
    NROW(data) > (embed_dim - 1) * lag_size
  )
  
  # If data is a vector, make into a tibble
  if (is.numeric(data)) {
    data <- tibble::tibble(x = data)
    col_name <- "x"
  }
  
  # Make a tibble of time series lags
  lag_tbl <- pbs_make_lags(data, col_name, embed_dim, lag_size)

  # Calculate Euclidean distances among row vectors
  lag_dist <- lag_tbl %>%
    as.matrix() %>%
    stats::dist(diag = FALSE, upper = TRUE) %>%
    broom::tidy() %>%
    tidyr::drop_na()

  # Initialize indices
  data_ind <- data %>% nrow() %>% seq_len()
  non_na_ind <- which(!is.na(rowSums(lag_tbl)))
  proj_ok_ind <- data_ind %>% intersect(non_na_ind - pred_dist)
  lib_ind <- lib_ind %>% intersect(non_na_ind) %>% intersect(proj_ok_ind)
  pred_to_ind <- pred_ind %>% intersect(non_na_ind) # Possibly from NAs
  pred_from_ind <- (pred_to_ind - pred_dist) %>% intersect(non_na_ind)
  
  # Instantiate prediction vector
  pred_vec <- rep(NA, NROW(data))

  # Iterate over prediction set
  for (time_ind in pred_from_ind) {
    exclude_ind <- seq(
      from = time_ind + pred_dist,
      by = lag_size,
      length.out = embed_dim
    )
    rel_lib_ind <- lib_ind %>% setdiff(time_ind) %>% setdiff(exclude_ind)
    
    # Identify nearest neighbours
    rel_lag_dist <- lag_dist %>%
      dplyr::filter(item2 %in% time_ind,
                    item1 %in% rel_lib_ind) %>%
      dplyr::arrange(distance) %>%
      dplyr::mutate(n = row_number()) %>%
      dplyr::filter(n <= embed_dim + 1)
    # Are there enough points?
    if (nrow(rel_lag_dist) > embed_dim) {
      # Calculate the weights
      omega_weights <- exp(-rel_lag_dist$distance / rel_lag_dist$distance[1])
      # Identify the iterated indicies
      iterated_ind <- rel_lag_dist$item1 + pred_dist
      # Predict the value
      pred_vec[time_ind + pred_dist] <- sum(
        dplyr::pull(data, col_name)[iterated_ind] * omega_weights
      )  / sum(omega_weights)
    } else {
      pred_vec[time_ind + pred_dist] <- NA
    }
  }
  pred_obs_cor <- stats::cor(pred_vec,
                             dplyr::pull(data, col_name),
                             use = "pairwise.complete.obs")
  # Return
  list(
    stats_tbl = tibble::tibble( # Needs variance etc.
      embed_dim = embed_dim,
      lag_size = lag_size,
      pred_dist = pred_dist,
      pred_obs_cor = pred_obs_cor
    ),
    pred_tbl = tibble::tibble(
      obs = data[col_name],
      pred = pred_vec
    ),
    nbr_list = NULL # List of tibbles of neighbour distances
  )
}
