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
                    lag,
                    dist,
                    symm,
                    from,
                    into,
                    stats) {
  
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
  # Make lag matrices
  
  # Calculate Euclidean distances
  
  # Specify global indices
  
  # Generate predictions
  
    # Specify local indices
  
    # Identify nearest neighbours
  
    # Calculate weights
  
    # Forecast value
  
  # Calculate statistics
  
  # Return data frame
  
}













