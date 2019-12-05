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
#' @examples pbs_s_map(tibble::tibble(x = 1:100), "x")
#' 
pbs_s_map <- function(data,
											col_name = NULL,
											embed_dim = 2L,
											local_weight = 0L,
											lag_size = 1L,
											pred_dist = 1L,
											lib_ind = seq_len(NROW(data)),
											pred_ind = seq_len(NROW(data))) {
	
	# Check arguments
	
	# Make sure data is a tibble
	
	# Create a tibble of time series lags
	lag_tbl <- pbs_make_lags(data, col_name, embed_dim, lag_size)
	
	# Calculate Euclidean distances among row vectors
	lag_dist <- pbs_calc_dist(lag_tbl)
	
	# 
}