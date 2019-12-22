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
#' data_frame <- data.frame(x = 1:30)
#' lags <- list(x = 0:3)
#' pbs_edm(data_frame, lags)
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
  neighbours_matrix <- t(apply(distance_matrix, 1, order))[, seq_len(num_nbrs)]

  # Calculate weights matrix
  
  # Make projected neighbour index matrix
  projected_matrix <- neighbours_matrix + forecast_distance
  
  # Make projected neighbour value matrix

  # Make forecasts
  forecasts <- rowSums(value_matrix * weight_matrix) / rowSums(weight_matrix)
  
  # Compute statistics
  
  # Return tibble
  
}


