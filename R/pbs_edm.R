#' Use EDM to Perform Out-of-Sample Forecasting 
#'
#' @param data A data frame with named columns (list of numeric vectors)
#' @param names Names of numeric columns in data (character vector)
#' @param dims Number of embedding dimensions for each column (numeric vector)
#' @param lag Number of time steps separating successive lags (numeric scalar)
#' @param dist Forecast distance (numeric scalar)
#' @param symm Symmetric exclusion radius? (logical scalar)
#' @param from Rows to predict from (numeric vector)
#' @param into Rows to predict into (numeric vector)
#' @param stats Return only stats and omit neighbour lists? (logical scalar)
#'
#' @return A data frame
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#'
#' @examples
#' 
pbs_edm <- function(data,
                    names,
                    dims,
                    lag = 1L,
                    dist = 1L,
                    symm = FALSE,
                    from = seq_len(nrow(mat)), 
                    into = seq_len(nrow(mat)), 
                    stats = TRUE) {
  
  # Check arguments
  stopifnot(
    is.data.frame(data),
    is.character(names),
    is.numeric(dims),
    is.vector(dims),
    is.numeric(lag),
    is.numeric(dist),
    is.logical(symm),
    is.numeric(from),
    is.vector(from),
    is.numeric(into),
    is.vector(into),
    is.logical(stats)
  )
  # Make lag matrix
  lag_matrix <- as.matrix(combine_lag_tibbles(data, names, dims, lag))
  
  # Calculate Euclidean distances
  dist_tibble <- make_dist_tibble(lag_matrix)
  
  # Specify global indices
  global_indices <- make_global_indices(lag_matrix, from, into, dist)
  
  # Generate forecasts
    
    # Specify local indices
  
    # Identify nearest neighbours
  
    # Calculate weights
  
    # Forecast value
  
  # Calculate statistics
  
  # Return data frame
  
}













