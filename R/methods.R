#' Perform Out-of-Sample Forecasting via Empirical Dynamic Modelling
#'
#' @param nt [data.frame()] A data frame with named columns for the raw (unlagged) variables
#' @param lags [list()] A named list of integer vectors specifying the lags for each variable
#' @param forecast_distance [integer(1)] The forecast distance
#' @param first_difference [logical(1)] First difference each variable before lagging?
#' @param centre_and_scale [logical(1)] Centre and scale each variable before lagging?
#'
#' @details Only lags of variables explicitly named in \code{lags} are used. The
#' observed time series whose values are forecast must be specified first. For
#' example, given a data frame with variables \code{Predator}, \code{Prey} and \code{Salinity}, 
#' the unlagged version of the variable \code{Predator} can be specified as the 
#' observed time series by 
#' 
#' \code{lags = list(Predator = c(0, ...), ...)}. 
#' 
#' The unlagged version of \code{Predator} can be forecast from the first and second lags 
#' of \code{Predator}, unlagged and second lag of \code{Prey}, and unlagged \code{Salinity} by
#' specifying:
#' 
#' \code{lags = list(Predator = c(0:2), Prey = c(0, 2), Salinity = c(0))}.
#' 
#' @return A list containing:
#' 
#' \itemize{
#' 	 \item nt_results = NULL,
#'   \item nt_observed = NULL,
#'   \item nt_forecast = NULL,
#'   \item xt_results [data.frame()] A summary of results
#'   \item xt_observed [numeric()] Observed (possibly transformed) with \code{length = nrow(nt)}
#'   \item xt_forecast [numeric()] Forecast (of possibly transformed) with \code{length = nrow(nt)}
#'   \item xt_lags [matrix()] A column matrix of lagged variables with \code{nrow = nrow(nt)}
#'   \item xt_distance [matrix()] Focal index (row), neighbour (column), distance (value)
#'   \item xt_nbr_index [matrix()] Focal index (row), rank (column), and neighbour index (value)
#'   \item xt_nbr_value [matrix()] Focal index (row), rank (column), and neighbour value (value)
#'   \item xt_nbr_distance [matrix()] Focal index (row), rank (column), and neighbour distance (value)
#'   \item xt_nbr_weight [matrix()] Focal index (row), rank (column), and neighbour weight (value)
#'   \item xt_prj_index [matrix()] Projected index (row), rank (column), and neighbour index (value)
#'   \item xt_prj_value [matrix()] Projected index (row), rank (column), and neighbour index (value)
#'   \item xt_prj_weight  [matrix()] Projected index (row), rank (column), and neighbour weight (value)
#'   \item forecast_distance = [integer(1)]
#'   \item first_difference = [logical(1)]
#'   \item centre_and_scale = [logical(1)]
#' }
#' 
#'
#' @author Luke A. Rogers
#' @export
#'
#' @examples
#' nt <- matrix(rep(1:30, 5), ncol = 5)
#' colnames(nt) <- c("A", "B", "C", "D", "E")
#' lags <- list(A = c(0, 1, 2), B = c(0, 1), C = c(0, 1, 2))
#' m1 <- pbsEDM(nt, lags)
#'
#' nt <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m2 <- pbsEDM(nt, lags)
#'
#' nt <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m3 <- pbsEDM(nt, lags, first_difference = TRUE)
#'
pbsEDM <- function (nt,
										lags,
										forecast_distance = 1L,
										first_difference = FALSE,
										centre_and_scale = FALSE) {

	#----------------- Check arguments ------------------------------------------#

	stopifnot(
		is.numeric(nt) || is.matrix(nt) || is.data.frame(nt),
		is.list(lags),
		is.integer(forecast_distance) && length(forecast_distance) == 1L,
		is.logical(first_difference) && length(first_difference) == 1L,
		is.logical(centre_and_scale) && length(centre_and_scale) == 1L
	)

	#----------------- Transform time series ------------------------------------#

	xt <- as.matrix(nt[, names(lags)])
	colnames(xt) <- names(lags)
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

	# xt_lags is a matrix of named lagged column vectors
	lags_size <- unlist(lags, use.names = FALSE)
	lags_name <- rep(names(lags), lengths(lags))
	xt_lags <- pbsLAG(xt[, lags_name], lags_size)
	colnames(xt_lags) <- paste0(lags_name, "_", lags_size)

	#----------------- Create distance matrix -----------------------------------#

	# xt_dist[i, j] is the distance between row vectors i and j in xt_lags
	xt_lags_na <- xt_lags
	xt_lags_na[which(is.na(rowSums(xt_lags_na))), ] <- NA # For distance
	xt_dist <- as.matrix(dist(xt_lags_na))

	#----------------- Exclude elements from the distance matrix ----------------#

	# Exclude the xt_lags focal row vector
	diag(xt_dist) <- NA

	# Exclude xt_lags row vectors that contain NAs
	na_rows <- which(is.na(rowSums(xt_lags)))
	xt_dist[na_rows, ] <- NA
	xt_dist[, na_rows] <- NA

	# Exclude xt_lags rows that project beyond xt_lags
	seq_rows <- seq_len(nrow(xt_lags))
	na_rows <- which((seq_rows + forecast_distance) > max(seq_rows))
	xt_dist[na_rows, ] <- NA
	xt_dist[, na_rows] <- NA

	# Exclude xt_lags rows that project to rows that contain NAs
	prj_rows <- (which(is.na(rowSums(xt_lags))) - forecast_distance)
	na_rows <- prj_rows[which(prj_rows > 0)]
	# xt_dist[na_rows, ] <- NA # Commented to allow forecast from here
	xt_dist[, na_rows] <- NA

	# Exclude xt_lags rows that contain a projection of the focal value
	rep_rows <- rep(seq_rows, each = length(lags[[1]]))
	prj_rows <- rep_rows + forecast_distance + lags[[1]]
	na_mat <- matrix(c(rep_rows, prj_rows), ncol = 2)[which(prj_rows <= nrow(xt_lags)),]
	xt_dist[na_mat] <- NA # Specify [row, col] pairs

	#----------------- Create neighbour index matrix ----------------------------#

	# nbr_inds is an nrow(xt_lags) x num_nbrs matrix of xt_lags row indices
	num_nbrs <- length(lags_size) + 1
	seq_nbrs <- seq_len(num_nbrs)
	nbr_inds <- t(apply(xt_dist, 1, order))[, seq_nbrs]
	nbr_inds[which(rowSums(!is.na(xt_dist)) < num_nbrs), ] <- NA

	#----------------- Create neighbour matrices --------------------------------#

	# nbr_vals is a matrix of values from xt_lags[, 1] corresponding to nbr_inds
	nbr_vals <- t(apply(nbr_inds, 1, function(x, y) y[x, 1], y = xt_lags))
	nbr_dist <- t(apply(xt_dist, 1, sort, na.last = T))[, seq_nbrs]
	nbr_wgts <- t(apply(nbr_dist, 1, function(x) exp(-x / x[1])))

	#----------------- Project neighbour matrices -------------------------------#

	prj_inds <- pbsLAG(nbr_inds, forecast_distance) + forecast_distance
	prj_vals <- t(apply(prj_inds, 1, function(x, y) y[x, 1], y = xt_lags))
	prj_wgts <- pbsLAG(nbr_wgts, forecast_distance)

	#----------------- Prepare return values ------------------------------------#

	xt_observed <- xt_lags[, 1]
	xt_forecast <- as.vector(rowSums(prj_vals * prj_wgts) / rowSums(prj_wgts))
	xt_rho <- cor(xt_observed, xt_forecast, use = "pairwise.complete.obs")
	xt_rmse <- sqrt(mean((xt_observed - xt_forecast)^2, na.rm = TRUE))
	xt_dim <- length(lags_size)
	xt_results <- data.frame(xt_dim = xt_dim,
													 xt_rho = xt_rho,
													 xt_rmse = xt_rmse,
													 stringsAsFactors = FALSE)

	#----------------- Return a list --------------------------------------------#

	structure(
		list(
			nt_results = NULL,
			nt_observed = NULL,
			nt_forecast = NULL,
			xt_results = xt_results,
			xt_observed = xt_observed,
			xt_forecast = xt_forecast,
			xt_lags = xt_lags,
			xt_distance = xt_dist,
			xt_nbr_index = nbr_inds,
			xt_nbr_value = nbr_vals,
			xt_nbr_distance = nbr_dist,
			xt_nbr_weight = nbr_wgts,
			xt_prj_index = prj_inds,
			xt_prj_value = prj_vals,
			xt_prj_weight = prj_wgts,
			forecast_distance = as.integer(forecast_distance),
			first_difference = first_difference,
			centre_and_scale = centre_and_scale
		),
		class = "pbsEDM"
	)
}


#' Perform Out-of-Sample Forecasting via S-Mapping
#'
#' @param nt [data.frame()] A data frame with named columns
#' @param lags [list()] A named list of integer vectors
#' @param local_weight [numeric(1)]
#' @param forecast_distance [integer(1)]
#' @param first_difference [logical(1)]
#' @param centre_and_scale [logical(1)]
#'
#' @details Only lags of variables explicitly named in \code{lags} are used. The
#' observed time series whose values are forecast must be specified first. For
#' example, given a data frame with variables \code{Predator}, \code{Prey} and \code{Salinity}, 
#' the unlagged version of the variable \code{Predator} can be specified as the 
#' observed time series by 
#' 
#' \code{lags = list(Predator = c(0, ...), ...)}. 
#' 
#' The unlagged version of \code{Predator} can be forecast from the first and second lags 
#' of \code{Predator}, unlagged and second lag of \code{Prey}, and unlagged \code{Salinity} by
#' specifying:
#' 
#' \code{lags = list(Predator = c(0:2), Prey = c(0, 2), Salinity = c(0))}.
#' 
#' @return A list containing:
#' 
#' \itemize{
#' 	 \item nt_results = NULL,
#'   \item nt_observed = NULL,
#'   \item nt_forecast = NULL,
#'   \item xt_results [data.frame()] A summary of results
#'   \item xt_observed [numeric()] Observed (possibly transformed) with \code{length = nrow(nt)}
#'   \item xt_forecast [numeric()] Forecast (of possibly transformed) with \code{length = nrow(nt)}
#'   \item xt_lags [matrix()] A column matrix of lagged variables with \code{nrow = nrow(nt)}
#'   \item xt_distance [matrix()] Focal index (row), neighbour (column), distance (value)
#'   \item xt_nbr_index [matrix()] Focal index (row), rank (column), and neighbour index (value)
#'   \item xt_nbr_distance [matrix()] Focal index (row), rank (column), and neighbour distance (value)
#'   \item xt_nbr_weight [matrix()] Focal index (row), rank (column), and neighbour weight (value)
#'   \item xt_prj_index [matrix()] Projected index (row), rank (column), and neighbour index (value)
#'   \item xt_prj_value [matrix()] Projected index (row), rank (column), and neighbour index (value)
#'   \item xt_prj_weight  [matrix()] Projected index (row), rank (column), and neighbour weight (value)
#'   \item forecast_distance = [integer(1)]
#'   \item first_difference = [logical(1)]
#'   \item centre_and_scale = [logical(1)]
#' }
#'
#' @author Luke A. Rogers
#' @export
#'
#' @examples
#' nt <- matrix(rep(1:30, 5), ncol = 5)
#' colnames(nt) <- c("A", "B", "C", "D", "E")
#' lags <- list(A = c(0, 1, 2), B = c(0, 1), C = c(0, 1, 2))
#' m1 <- pbsEDM(nt, lags)
#'
#' nt <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m2 <- pbsSMAP(nt, lags)
#'
pbsSMAP <- function (nt,
										 lags,
										 local_weight = 0,
										 forecast_distance = 1L,
										 first_difference = FALSE,
										 centre_and_scale = FALSE) {

	#----------------- Check arguments ------------------------------------------#

	stopifnot(
		is.numeric(nt) || is.matrix(nt) || is.data.frame(nt),
		is.list(lags),
		is.numeric(local_weight) && length(local_weight) == 1L,
		is.integer(forecast_distance) && length(forecast_distance) == 1L,
		is.logical(first_difference) && length(first_difference) == 1L,
		is.logical(centre_and_scale) && length(centre_and_scale) == 1L
	)

	#----------------- Transform time series ------------------------------------#

	xt <- as.matrix(nt[, names(lags)])
	colnames(xt) <- names(lags)
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

	# xt_lags is a matrix of named lagged column vectors
	lags_size <- unlist(lags, use.names = FALSE)
	lags_name <- rep(names(lags), lengths(lags))
	xt_lags <- pbsLAG(xt[, lags_name], lags_size)
	colnames(xt_lags) <- paste0(lags_name, "_", lags_size)

	#----------------- Create distance matrix -----------------------------------#

	# xt_dist[i, j] is the distance between row vectors i and j in xt_lags
	xt_lags_na <- xt_lags
	xt_lags_na[which(is.na(rowSums(xt_lags_na))), ] <- NA # For distance
	xt_dist <- as.matrix(dist(xt_lags_na))

	#----------------- Exclude elements from the distance matrix ----------------#

	# Exclude the xt_lags focal row vector
	diag(xt_dist) <- NA

	# Exclude xt_lags row vectors that contain NAs
	na_rows <- which(is.na(rowSums(xt_lags)))
	xt_dist[na_rows, ] <- NA
	xt_dist[, na_rows] <- NA

	# Exclude xt_lags rows that project beyond xt_lags
	seq_rows <- seq_len(nrow(xt_lags))
	na_rows <- which((seq_rows + forecast_distance) > max(seq_rows))
	xt_dist[na_rows, ] <- NA
	xt_dist[, na_rows] <- NA

	# Exclude xt_lags rows that project to rows that contain NAs
	prj_rows <- (which(is.na(rowSums(xt_lags))) - forecast_distance)
	na_rows <- prj_rows[which(prj_rows > 0)]
	# xt_dist[na_rows, ] <- NA
	xt_dist[, na_rows] <- NA

	# Exclude xt_lags rows that contain a projection of the focal value
	rep_rows <- rep(seq_rows, each = length(lags[[1]]))
	prj_rows <- rep_rows + forecast_distance + lags[[1]]
	na_mat <- matrix(c(rep_rows, prj_rows), ncol = 2)[which(prj_rows <= nrow(xt_lags)),]
	xt_dist[na_mat] <- NA # Specify [row, col] pairs

	#----------------- Create neighbour index matrix ----------------------------#

	nbr_dist <- t(apply(xt_dist, 1, sort, na.last = TRUE))
	nbr_inds <- t(apply(xt_dist, 1, order))
	nbr_inds[which(is.na(nbr_dist))] <- NA
	# nbr_vals <- t(apply(nbr_inds, 1, function(x, y) y[x, 1], y = xt_lags))
	nbr_wgts <- t(apply(nbr_dist, 1,
											function(x, y) exp(-y * x / mean(x, na.rm = TRUE)),
											y = local_weight))

	#----------------- Compute lag of neighbour index matrix --------------------#

	# TODO: Needed?
	lag_inds <- pbsLAG(nbr_inds, forecast_distance)

	#----------------- Project neighbour matrices -------------------------------#

	prj_inds <- pbsLAG(nbr_inds, forecast_distance) + forecast_distance
	prj_vals <- t(apply(prj_inds, 1, function(x, y) y[x, 1], y = xt_lags))
	prj_wgts <- pbsLAG(nbr_wgts, forecast_distance)

	#----------------- Project xt_lag matrix ------------------------------------#

	prj_lags <- pbsLAG(xt_lags, forecast_distance)

	#----------------- Compute B matrix for SVD ---------------------------------#

	# The row gives the focal index
	# The col gives the nearest neighbours ordered relative to focal index
	b_matrix <- prj_wgts * prj_vals
	b_matrix[which(is.na(b_matrix))] <- 0

	#----------------- Compute W array of matrices for SVD ----------------------#

	# The row (first dimension) gives the nearest neighbours relative to focal
	# The (second dimension) gives the xt_lags row vector index
	# The col (third dimension) gives the focal index
	w_array <- sapply(X = seq_rows,
										FUN = function(X, w, y) w[X, ] %*% t(rep(1, y)),
										w = prj_wgts,
										y = length(lags_size),
										simplify = "array")

	#----------------- Compute L array of lagged row vectors for SVD ------------#

	# The row (first dimension) gives the nearest neighbours relative to focal
	# The (second dimension) gives the xt_lags row vector index
	# The col (third dimension) gives the focal index
	l_array <- sapply(X = seq_rows,
										FUN = function(X, l, m) l[m[X, ], ],
										l = xt_lags,
										m = lag_inds, # Double check
										simplify = "array")

	#----------------- Compute A array of matrices for SVD ----------------------#

	# The row (first dimension) gives the nearest neighbours relative to focal
	# The (second dimension) gives the xt_lags row vector index
	# The col (third dimension) gives the focal index
	a_array <- w_array * l_array
	a_array[which(is.na(a_array))] <- 0

	#----------------- Solve for C matrix via SVD -------------------------------#

	# Decompose A matrices by SVD
	svd_list <- apply(a_array, 3, svd)

	# Simplify
	vdu_array <- sapply(
		X = seq_rows,
		FUN = function(X, s) s[[X]]$v %*% diag(1/s[[X]]$d) %*% t(s[[X]]$u),
		s = svd_list,
		simplify = "array")

	# Solve for C matrix
	c_matrix <- sapply(X = seq_rows,
										 FUN = function(X, a, b) a[,, X] %*% b[X, ],
										 a = vdu_array,
										 b = b_matrix)

	#----------------- Make forecasts -------------------------------------------#

	xt_forecast <- sapply(X = seq_rows,
												FUN = function(X, l, m) sum(m[, X] * l[X, ]),
												m = c_matrix,
												l = prj_lags)
	xt_forecast[is.nan(xt_forecast)] <- NA_real_

	#----------------- Prepare return values ------------------------------------#

	xt_observed <- xt_lags[, 1]
	xt_rho <- cor(xt_observed, xt_forecast, use = "pairwise.complete.obs")
	xt_rmse <- sqrt(mean((xt_observed - xt_forecast)^2, na.rm = TRUE))
	xt_dim <- length(lags_size)
	xt_theta <- local_weight
	xt_results <- data.frame(xt_dim = xt_dim,
													 xt_theta = xt_theta,
													 xt_rho = xt_rho,
													 xt_rmse = xt_rmse,
													 stringsAsFactors = FALSE)

	#----------------- Return a list --------------------------------------------#

	structure(
		list(
			nt_results = NULL,
			nt_observed = NULL,
			nt_forecast = NULL,
			xt_results = xt_results,
			xt_observed = xt_observed,
			xt_forecast = xt_forecast,
			xt_lags = xt_lags,
			xt_distance = xt_dist,
			xt_nbr_index = nbr_inds,
			xt_nbr_distance = nbr_dist,
			xt_nbr_weight = nbr_wgts,
			xt_prj_index = prj_inds,
			xt_prj_value = prj_vals,
			xt_prj_weight = prj_wgts,
			forecast_distance = as.integer(forecast_distance),
			first_difference = first_difference,
			centre_and_scale = centre_and_scale
		),
		class = "pbsEDM"
	)
}
