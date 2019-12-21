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
#' 
pbs_edm <- function(data_frame,
                    lags,
                    interpret_lags = NULL,
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
  
  # Make neighbour index matrix
  
  # Calculate weights vector
  
  # Make projected neighbour index matrix
  
  # Make projected neighbour value matrix
  
  # Make forecasts
  forecasts <- rowSums(value_matrix * weight_matrix) / rowSums(weight_matrix)
  
  # Compute statistics
  
  # Return tibble
  
}


