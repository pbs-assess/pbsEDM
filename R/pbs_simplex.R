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
                         col_name = "x",
                         embed_dim = 2L,
                         lag_size = 1L,
                         pred_dist = 1L,
                         lib_ind = seq_len(NROW(data)),
                         pred_ind = seq_len(NROW(data))) {
  
  # Create tibble 'data_tbl'
  data_tbl <- pbs_make_tibble(data, col_name)
  
  # Check arguments
  assertthat::assert_that(
    assertthat::is.count(embed_dim),
    assertthat::is.count(lag_size),
    assertthat::is.count(pred_dist),
    is.numeric(lib_ind),
    is.vector(lib_ind),
    is.numeric(pred_ind),
    is.vector(pred_ind),
    nrow(data_tbl) > (embed_dim - 1) * lag_size
  )

  # Create a tibble of time series lags
  lag_tbl <- pbs_make_lags(data_tbl, col_name, embed_dim, lag_size)

  # Calculate Euclidean distances among row vectors
  lag_dist <- pbs_calc_dist(lag_tbl)

  # Initialize indices
  data_ind <- data_tbl %>% nrow() %>% seq_len()
  non_na_ind <- which(!is.na(rowSums(lag_tbl)))
  proj_ok_ind <- data_ind %>% intersect(non_na_ind - pred_dist)
  lib_ind <- lib_ind %>% intersect(non_na_ind) %>% intersect(proj_ok_ind)
  pred_to_ind <- pred_ind %>% intersect(non_na_ind) # Possibly from NAs
  pred_from_ind <- (pred_to_ind - pred_dist) %>% intersect(non_na_ind)
  
  # Instantiate prediction vector and neighbour list
  pred_vec <- rep(NA, nrow(data_tbl))
  nbr_list <- list()

  # Iterate over prediction set
  for (time_ind in pred_from_ind) {
    
    # Exclude invalid indices
    exclude_ind <- seq(
      from = time_ind + pred_dist,
      by = lag_size,
      length.out = embed_dim
    )
    
    # Specify the valid indices
    rel_lib_ind <- lib_ind %>% setdiff(time_ind) %>% setdiff(exclude_ind)
    
    # Identify nearest neighbours
    rel_lag_dist <- pbs_make_nbrs(lag_dist, 
                                  time_ind, 
                                  rel_lib_ind, 
                                  embed_dim, 
                                  pred_dist)
    
    # Concatenate to the neighbour list
    nbr_list[[time_ind]] <- rel_lag_dist
    
    # Are there enough points?
    if (nrow(rel_lag_dist) > embed_dim) {
      
      # Pull vectors
      proj_ind <- dplyr::pull(rel_lag_dist, proj_nbr_ind)
      proj_vals <- dplyr::pull(data_tbl, col_name)[proj_ind]
      weight <- dplyr::pull(rel_lag_dist, weight)

      # Predict the value
      pred_vec[time_ind + pred_dist] <- sum(proj_vals * weight) / sum(weight) 

    } else {
      
      # Prediction is NA
      pred_vec[time_ind + pred_dist] <- NA
      
    }
  }
  
  # Caluculate the prediction statistics
  obs_vec <- dplyr::pull(data_tbl, col_name)
  pred_obs_cor <- stats::cor(obs_vec, pred_vec, use = "pairwise.complete.obs")
  pred_obs_rmse <- sqrt(mean((obs_vec - pred_vec)^2, na.rm = TRUE))

  # Return a list
  list(
    stats_tbl = tibble::tibble( # Needs variance etc.
      embed_dim = embed_dim,
      lag_size = lag_size,
      pred_dist = pred_dist,
      pred_obs_cor = pred_obs_cor,
      rmse = pred_obs_rmse
    ),
    pred_tbl = tibble::tibble(
      obs = data_tbl[col_name],
      pred = pred_vec
    ),
    nbr_list = nbr_list # List of tibbles of neighbour distances
  )
}
