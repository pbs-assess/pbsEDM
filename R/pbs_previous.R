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
#' @param one_sided Should only the focal value be excluded? (logical)
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
#' @examples pbs_simplex_1(simple_ts, "x")
#' 
pbs_simplex_1 <- function (data,
                         col_name = "x",
                         embed_dim = 2L,
                         lag_size = 1L,
                         pred_dist = 1L,
                         one_sided = TRUE,
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
    exclude_ind <- util_exclude_indices(time_ind, 
                                        pred_dist, 
                                        lag_size, 
                                        embed_dim, 
                                        one_sided = one_sided)
    
    # Specify the valid indices
    rel_ind <- lib_ind %>% setdiff(time_ind) %>% setdiff(exclude_ind)
    
    # Identify nearest neighbours
    nbr_dist <- pbs_make_nbrs(lag_dist, 
                              time_ind, 
                              rel_ind, 
                              embed_dim, 
                              pred_dist)
    
    # Augment by weights
    nbr_dist <- nbr_dist %>% 
      dplyr::mutate(weight = exp(-distance / distance[1]))
    
    # Concatenate to the neighbour list
    nbr_list[[time_ind]] <- nbr_dist
    
    # Store the number of neighbours
    num_nbrs <- nrow(nbr_dist)
    
    # Are there enough neighbours?
    if (num_nbrs > embed_dim) {
      
      # Pull vectors
      proj_ind <- dplyr::pull(nbr_dist, nbr_proj)
      proj_vals <- dplyr::pull(data_tbl, col_name)[proj_ind]
      weight <- dplyr::pull(nbr_dist, weight)
      
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

#' Use An S-Map to Perform Out-of-Sample Forecasting
#'
#' @param data A time series with possible NA values (numeric vector or
#'   tibble with a named column).
#' @param col_name If `data` is a tibble, the time series column name 
#'   (character scalar).
#' @param embed_dim The embedding dimension for the lagged-coordinate 
#'   embedding of the time series (integer scalar)
#' @param local_weight The exponential local weighting parameter
#'   (integer scalar)
#' @param lag_size The number of time steps separating successive lags
#'   (integer scalar)
#' @param pred_dist The number of time steps forward to forecast
#'   (integer scalar)
#' @param one_sided Should only the focal value be excluded? (logical)
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
#' @examples pbs_s_map_1(simple_ts, "x", local_weight = 2)
#' 
pbs_s_map_1 <- function(data,
                      col_name = "x",
                      embed_dim = 2L,
                      local_weight = 0L,
                      lag_size = 1L,
                      pred_dist = 1L,
                      one_sided = TRUE,
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
    exclude_ind <- util_exclude_indices(time_ind, 
                                        pred_dist, 
                                        lag_size, 
                                        embed_dim, 
                                        one_sided = one_sided)
    
    # Specify the valid indices
    rel_ind <- lib_ind %>% setdiff(time_ind) %>% setdiff(exclude_ind)
    
    # Identify nearest neighbours
    nbr_dist <- pbs_make_nbrs(lag_dist, 
                              time_ind, 
                              rel_ind, 
                              embed_dim, 
                              pred_dist,
                              max_nbrs = nrow(data_tbl))
    
    # Augment by weights
    nbr_dist <- nbr_dist %>% 
      dplyr::mutate(
        weight = exp(-local_weight * distance / mean(distance, na.rm = TRUE))
      )
    
    # Concatenate to the neighbour list
    nbr_list[[time_ind]] <- nbr_dist
    
    # Store the number of neighbours
    num_nbrs <- nrow(nbr_dist)
    
    # Are there enough neighbours?
    if (num_nbrs > embed_dim) {
      
      # Pull vectors
      proj_ind <- dplyr::pull(nbr_dist, nbr_proj)
      proj_vals <- dplyr::pull(data_tbl, col_name)[proj_ind]
      
      # Calculate the B vector and A matrix
      a_mat <- matrix(NA_real_, nrow = num_nbrs, ncol = embed_dim)
      b_vec <- rep(NA_real_, num_nbrs)
      for (i in seq_len(num_nbrs)) {
        b_vec[i] <- nbr_dist$weight[i] * proj_vals[i]
        for (j in seq_len(embed_dim)) {
          a_mat[i, j] <- nbr_dist$weight[i] * lag_tbl[[nbr_dist$nbr_ind[i], j]]
        }
      }
      
      # Solve for C using Singluar Value Decomposition
      a_svd <- svd(a_mat) # m x n
      d_inv <- diag(1 / a_svd$d) # diag n x n
      u_mat <- a_svd$u # m x n
      v_mat <- a_svd$v # n x n
      c_vec <- v_mat %*% d_inv %*% t(u_mat) %*% b_vec
      
      # Predict the value
      pred_vec[time_ind + pred_dist] <- sum(c_vec * lag_tbl[time_ind, ])
      
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
      local_weight = local_weight,
      lag_size = lag_size,
      pred_dist = pred_dist,
      pred_obs_cor = pred_obs_cor,
      rmse = pred_obs_rmse
    ),
    pred_tbl = tibble::tibble(
      obs = obs_vec,
      pred = pred_vec
    ),
    nbr_list = nbr_list # List of tibbles of neighbour distances
  )  
}


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
#'   pbs_edm_1(data_frame, lags)
#' 
pbs_edm_1 <- function(data_frame,
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


#' Use S-Mapping to Perform Out-of-Sample Forecasting
#'
#' @param data_frame A data frame with named columns (list of numeric vectors)
#' @param lags Named list of numeric vectors giving lags to use for each
#'     column. List names must match column names. (list of numeric vectors)
#' @param local_weight State-dependence parameter (numeric scalar)
#' @param from_user Rows to predict from (numeric vector)
#' @param into_user Rows to predict into (numeric vector)
#' @param forecast_distance Forecast distance (numeric scalar)
#' @param symmetric_exclusion Symmetric exclusion radius? (logical scalar)
#' @param include_stats Return forecast stats? (logical scalar)
#' @param include_forecasts Return forecast and observations? (logical scalar)
#' @param include_neighbours Return neighbours? (logical scalar)
#'
#' @return A tibble (work in progress)
#' @export
#'
#' @examples 
#'   data_frame <- data.frame(x = simple_ts)
#'   lags <- list(x = 0:3)
#'   local_weight <- 0
#'   pbs_smap_1(data_frame, lags, local_weight)
#' 
pbs_smap_1 <- function(data_frame,
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













