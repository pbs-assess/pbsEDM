#' Use S-Mapping to Perform Out-of-Sample Forecasting 
#'
#' @param data_frame A data frame with named columns (list of numeric vectors)
#' @param lags Named list of numeric vectors giving lags to use for each
#'     column. List names must match column names. (list of numeric vectors)
#' @param local_weight State-dependence parameter (numeric scalar)
#' @param forecast_distance Forecast distance (numeric scalar)
#' @param symmetric_exclusion Symmetric exclusion radius? (logical scalar)
#' @param include_stats Return forecast stats? (logical scalar)
#' @param include_forecasts Return forecast and observations? (logical scalar)
#' @param include_neighbours Return neighbours? (logical scalar)
#'
#' @return A tibble with results, forecasts, neighbours, distances, and weights
#' 
#' @export
#'
#' @examples
#' data_frame <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' pbs_smap(data_frame, lags)
#' 
pbs_smap <- function(data_frame,
                     lags,
                     local_weight = 0,
                     # from_user = seq_len(nrow(data_frame)), 
                     # into_user = seq_len(nrow(data_frame)), 
                     forecast_distance = 1L,
                     symmetric_exclusion = FALSE) {
  
  # Check arguments
  stopifnot(
    is.data.frame(data_frame),
    is.list(lags),
    is.numeric(local_weight),
    all(is.element(names(lags), names(data_frame))),
    is.numeric(forecast_distance),
    is.logical(symmetric_exclusion)
  )
  
  # Make lags matrix
  lag_sizes_vector <- unlist(lags, use.names = FALSE)
  col_names_vector <- rep(names(lags), lengths(lags))
  cols_list <- mapply(FUN = dplyr::pull,
                      var = col_names_vector,
                      MoreArgs = list(.data = data_frame),
                      SIMPLIFY = FALSE)
  lags_matrix <- mapply(FUN = dplyr::lag,
                        x = cols_list,
                        n = lag_sizes_vector)
  observations <- lags_matrix[, 1]
  
  # Make distance matrix
  lags_matrix[which(is.na(rowSums(lags_matrix))), ] <- NA # For distance
  distance_matrix <- as.matrix(dist(lags_matrix, diag = TRUE, upper = TRUE))
  diag(distance_matrix) <- NA
  
  # Exclude indices that project beyond the time series
  threshold <- ncol(distance_matrix) - forecast_distance
  distance_matrix[, which(seq_len(ncol(distance_matrix)) > threshold)] <- NA
  distance_matrix[which(seq_len(ncol(distance_matrix)) > threshold), ] <- NA
  
  # Exclude indices that contain the projection of the focal value
  # TODO: Check whether lags_unique is correct
  indices <- seq_len(nrow(distance_matrix))
  lags_unique <- unique(lag_sizes_vector)
  focal_indices <- rep(indices, each = length(lags_unique))
  exclude_indices <- focal_indices + forecast_distance + lags_unique
  within_range <- which(exclude_indices %in% indices)
  exclude_matrix <- as.matrix(data.frame(x = focal_indices[within_range],
                                         y = exclude_indices[within_range]))
  distance_matrix[exclude_matrix] <- NA
  
  # Exclude indices that project to a row that contains NAs
  na_indices <- which(is.na(rowSums(lags_matrix)))
  exclude_indices <- na_indices - forecast_distance
  within_range <- which(exclude_indices %in% indices)
  focal_indices <- rep(indices, each = length(within_range))
  exclude_matrix <- as.matrix(data.frame(x = focal_indices,
                                         y = rep(exclude_indices[within_range],
                                                 length(focal_indices))))
  distance_matrix[exclude_matrix] <- NA

  # Exclude indices symmetrically around the forecast index
  if (symmetric_exclusion == TRUE) {
    symm_na_indices <- focal_indices + forecast_distance - lags_unique
    within_range <- which(symm_na_indices %in% indices)
    symm_na_matrix <- as.matrix(data.frame(x = focal_indices[within_range],
                                           y = symm_na_indices[within_range]))
    distance_matrix[symm_na_matrix] <- NA
  }
  
  # Make neighbour distance matrix
  nbr_dst_matrix <- t(apply(distance_matrix, 1, sort, na.last = TRUE))
  
  # Make neighbour index matrix
  nbr_ind_matrix <- t(apply(distance_matrix, 1, order))
  nbr_ind_matrix[which(is.na(nbr_dst_matrix))] <- NA
  
  # Calculate weight matrix: row gives focal index, column gives ordered nbr
  weight_matrix <- t(apply(nbr_dst_matrix, 1, 
                           function(x, y) exp(-y * x / mean(x, na.rm = TRUE)),
                           y = local_weight))
  
  # Make projected neighbour index matrix
  lag_ind_matrix <- apply(nbr_ind_matrix, 2, dplyr::lag, n = forecast_distance)
  pro_ind_matrix <- lag_ind_matrix + forecast_distance
  
  # Make projected neighbour value matrix: row gives focal, column gives ord nbr
  pro_val_matrix <- t(apply(pro_ind_matrix, 1, function(x, y) y[x, 1],
                            y = lags_matrix)) # Works because NAs are NA_real_
  
  # Make projected distances matrix
  pro_dst_matrix <- apply(nbr_dst_matrix, 2, dplyr::lag, n = forecast_distance)
  
  # Make projected neighbour weights matrix
  pro_wts_matrix <- apply(weight_matrix, 2, dplyr::lag, n = forecast_distance)
  
  # Make projected lag matrix
  pro_lag_matrix <- apply(lags_matrix, 2, dplyr::lag, n = forecast_distance)
  
  # Calculate the B matrix of column vectors. Dimensions correspond to:
  # - Focal index
  # - Nearest neighbours ordered relative to focal index
  b_matrix <- pro_wts_matrix * pro_val_matrix
  
  # Replace NAs by zeros
  b_matrix[which(is.na(b_matrix))] <- 0
  
  # Compute the W array of matrices. Dimensions correspond to:
  # - Nearest neighbours ordered relative to focal index
  # - Lagged row vector index
  # - Focal index
  w_array <- sapply(X = seq_len(nrow(data_frame)), 
                    FUN = function(X, w, y) w[X, ] %*% t(rep(1, y)), 
                    w = pro_wts_matrix,
                    y = length(lags_unique),
                    simplify = "array")
  
  # Compute the L array of lagged row vectors. Dimensions correspond to:
  # - Nearest neighbours ordered relative to focal index
  # - Lagged row vector index
  # - Focal index
  l_array <- sapply(X = seq_len(nrow(data_frame)),
                    FUN = function(X, l, m) l[m[X, ], ],
                    l = lags_matrix,
                    m = lag_ind_matrix, # Double check
                    simplify = "array")
  
  # Compute the A array of matrices. Dimensions correspond to:
  # - Nearest neighbours ordered relative to focal index
  # - Lagged row vector index
  # - Focal index
  a_array <- w_array * l_array
  
  # Replace NAs by zeros
  a_array[which(is.na(a_array))] <- 0
  
  # Decompose A matrices by SVD
  svd_list <- apply(a_array, 3, svd)
  
  # Simplify
  vdu_array <- sapply(
    X = seq_len(nrow(data_frame)),
    FUN = function(X, s) s[[X]]$v %*% diag(1/s[[X]]$d) %*% t(s[[X]]$u),
    s = svd_list,
    simplify = "array")
  
  # Solve for C matrix
  c_matrix <- sapply(X = seq_len(nrow(data_frame)),
                     FUN = function(X, a, b) a[,, X] %*% b[X, ],
                     a = vdu_array,
                     b = b_matrix)

  # Make forecasts
  forecasts <- sapply(X = seq_len(nrow(data_frame)),
                      FUN = function(X, l, m) sum(m[, X] * l[X, ]),
                      m = c_matrix,
                      l = pro_lag_matrix)

  # Compute statistics
  rho <- cor(observations, forecasts, use = "pairwise.complete.obs")
  rmse <- sqrt(mean((observations - forecasts)^2, na.rm = TRUE))
  
  # Return tibble
  tibble::tibble(E = length(unlist(lags)),
                 theta = local_weight, 
                 rho = rho, 
                 rmse = rmse,
                 lags = lags,
                 observations = list(observations = observations),
                 forecasts = list(forecasts = forecasts),
                 neighbours = list(neighbours = pro_ind_matrix), 
                 distances = list(distances = pro_dst_matrix),
                 weights = list(weights = pro_wts_matrix))
}
