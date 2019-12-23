#' Title
#'
#' @param data_frame 
#' @param lags 
#' @param interpret_lags 
#' @param from_user 
#' @param into_user 
#' @param forecast_distance 
#' @param symmetric_exclusion 
#' @param include_stats 
#' @param include_forecasts 
#' @param include_neighbours 
#'
#' @return
#' @export
#'
#' @examples
#' data_frame <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' pbs_edm(data_frame, lags)
#' 
pbs_edm <- function(data_frame,
                    lags,
                    # from_user = seq_len(nrow(data_frame)), 
                    # into_user = seq_len(nrow(data_frame)), 
                    forecast_distance = 1L,
                    symmetric_exclusion = FALSE,
                    include_stats = TRUE,
                    include_forecasts = TRUE,
                    include_neighbours = TRUE) {
  
  # Check arguments
  
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
  
  # Make distance matrix
  lags_matrix[which(is.na(rowSums(lags_matrix))), ] <- NA # For correct distance
  distance_matrix <- as.matrix(dist(lags_matrix, diag = TRUE, upper = TRUE))
  diag(distance_matrix) <- NA
  
  # Exclude indices that project beyond the time series
  threshold <- ncol(distance_matrix) - forecast_distance
  distance_matrix[, which(seq_len(ncol(distance_matrix)) > threshold)] <- NA
  distance_matrix[which(seq_len(ncol(distance_matrix)) > threshold), ] <- NA
  
  # Exclude indices that contain the projection of the focal value
  indices <- seq_len(nrow(distance_matrix))
  lags_unique <- unique(lag_sizes_vector)
  focal_indices <- rep(indices, each = length(lags_unique))
  exclude_indices <- focal_indices + forecast_distance + lags_unique
  within_range <- which(exclude_indices %in% indices)
  exclude_matrix <- as.matrix(data.frame(x = focal_indices[within_range],
                                         y = exclude_indices[within_range]))
  distance_matrix[exclude_matrix] <- NA
  
  # Exclude indices symmetrically around the forecast index
  if (symmetric_exclusion == TRUE) {
    symm_na_indices <- focal_indices + forecast_distance - lags_unique
    within_range <- which(symm_na_indices %in% indices)
    symm_na_matrix <- as.matrix(data.frame(x = focal_indices[within_range],
                                           y = symm_na_indices[within_range]))
    distance_matrix[symm_na_matrix] <- NA
  }

  # Make neighbour index matrix
  num_nbrs <- length(lags_unique) + 1
  nbr_ind_matrix <- t(apply(distance_matrix, 1, order))[, seq_len(num_nbrs)]
  nbr_ind_matrix[which(rowSums(!is.na(distance_matrix)) < num_nbrs), ] <- NA
  
  # Make neighbour value matrix
  nbr_val_matrix <- t(apply(nbr_ind_matrix, 1, function(x, y) y[x, 1],
                            y = lags_matrix)) # Works because NAs are NA_real_
  
  # Calculate weight matrix
  weight_matrix <- t(apply(nbr_val_matrix, 1, function(x) exp(-x / x[1])))
  
  # Make projected neighbour index matrix
  pro_ind_matrix <- apply(nbr_ind_matrix, 2, dplyr::lag, n = forecast_distance)
  pro_ind_matrix <- pro_ind_matrix + forecast_distance
  
  # Make projected neighbour value matrix
  pro_val_matrix <- t(apply(pro_ind_matrix, 1, function(x, y) y[x, 1],
                            y = lags_matrix)) # Works because NAs are NA_real_
  
  # Make projected neighbour weight matrix
  pro_wt_matrix <- apply(weight_matrix, 2, dplyr::lag, n = forecast_distance)

  # Make forecasts
  forecasts <- rowSums(pro_val_matrix * pro_wt_matrix) / rowSums(pro_wt_matrix)
  
  # Compute statistics
  rho <- cor(lags_matrix[, 1], forecasts, use = "pairwise.complete.obs")
  rmse <- sqrt(mean((lags_matrix[, 1] - forecasts)^2, na.rm = TRUE))
  
  # Return tibble
  tibble::tibble(lags = lags, rho = rho, rmse = rmse,
                 forecasts = list(
                   tibble::tibble(observation = lags_matrix[, 1],
                                  forecast = forecasts)),
                 neighbours = list(nbr_ind_matrix), 
                 distances = list(distance_matrix))
}


