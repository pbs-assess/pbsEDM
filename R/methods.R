#' Perform Out-of-Sample Forecasting via Empirical Dynamic Modelling
#'
#' @param nt [data.frame()] A data frame with named columns
#' @param lags [list()] A named list of integer vectors
#' @param forecast_distance [integer(1)] The forecast distance
#' @param first_difference [logical(1)] 
#' @param centre_and_scale [logical(1)]
#' @param show_calculations [logical(1)]
#' 
#' @details Only lags and columns explicitly named in `lags` are used. For
#' example, an unlagged time series in a column named `Salinity` could be 
#' included via `lags = list(..., Salinity = c(0))`.
#'
#' @return A list
#' @export
#'
#' @examples
#' xt <- matrix(rep(1:10, 5), ncol = 5)
#' colnames(xt) <- c("A", "B", "C", "D", "E")
#' lags <- list(A = c(0), B = c(0, 1), C = c(0, 1, 2))
#' 
#' 
pbsEDM <- function (nt,
										lags,
										forecast_distance,
										first_difference,
										centre_and_scale,
										show_calculations) {
	
	#----------------- Check arguments ------------------------------------------#
	
	
	
	#----------------- Transform time series ------------------------------------#

	xt <- as.matrix(nt[, names(lags)])
	# First difference and buffer with NAs
	if (first_difference) {
		xt <- rbind(apply(xt, 2, diff), NA_real_)
	}
	# Centre and scale
	if (centre_and_scale) {
		xt_means <- apply(xt, 2, mean, na.rm = TRUE)
		xt_sds <- apply(xt, 2, sd, na.rm = TRUE)
		xt <- t((t(xt) - xt_means) / xt_sds)
	}
	
	#----------------- Create lagged matrix -------------------------------------#
	
	lags_size <- unlist(lags, use.names = FALSE)
	lags_name <- rep(names(lags), lengths(lags))
	xt_lags <- pbsLAG(xt[, lags_name], lags_size)
	colnames(xt_lags) <- paste0(lags_name, "_", lags_size)
	
	#----------------- Create computation objects -------------------------------#
	
	
	
	#----------------- Prepare return values ------------------------------------#
	
	
	
	#----------------- Return a list --------------------------------------------#
	
}