#' Create N Matrix
#'
#' @param N [matrix()] or [data.frame()] with [numeric()] columns.
#' @param lags [list()] of named integer vectors specifying the lags to use for
#'   each time series in \code{N}.
#' @param p The integer forecast distance.
#'
#' @return [matrix()] N
#' @export
#'
#' @examples
#' N <- data.frame(x = 1:10, y = 11:20)
#' lags <- list(x = c(0, 1, 2), y = c(0, 1))
#' pbsN(N, lags)
#' 
pbsN <- function (N, lags, p = 1L) {
	# Check arguments
	stopifnot(
		is.matrix(N) || is.data.frame(N),
		is.list(lags),
		all(is.element(names(lags), colnames(N))),
		length(unique(names(lags))) == length(names(lags)),
		length(unique(colnames(N))) == length(colnames(N)),
		is.numeric(as.vector(unlist(N[, names(lags)]))),
		is.numeric(as.vector(unlist(lags))),
		lags[[1]][1] == 0L,
		is.integer(p) && length(p) == 1L
	)
	# Compute N
	N <- as.matrix(N[, names(lags)]) # Retain only columns named in lags
	N <- rbind(N, array(NA_real_, dim = c(p, ncol(N)))) # Augment by rows
	colnames(N) <- names(lags) # Assign column names to new N
	return(N)	# Return new N
}

#' Create Z Matrix
#'
#' @param N [matrix()] with [numeric()] columns.
#' @param first_difference [logical()] First difference columns of N?
#'
#' @return [matrix()] Z
#' @export
#'
#' @examples
#' N <- data.frame(x = 1:10, y = 11:20)
#' lags <- list(x = c(0, 1, 2), y = c(0, 1))
#' N <- pbsN(N, lags)
#' Z <- pbsZ(N, first_difference = FALSE)
#' 
pbsZ <- function (N, first_difference) {
	# Check arguments
	stopifnot(
		is.matrix(N),
		is.logical(first_difference) && length(first_difference) == 1L
	)
	# Compute Z
	if (first_difference) {
		diff(N)
	} else {
		N
	}
}

#' Create Y Matrix
#'
#' @param Z [matrix()] with [numeric()] columns.
#' @param centre_and_scale [logical()] Centre and scale columns of Z?
#'
#' @return [matrix()] Y
#' @export
#'
pbsY <- function (Z, centre_and_scale) {
	# Check arguments
	stopifnot(
		is.matrix(Z),
		is.logical(centre_and_scale) && length(centre_and_scale) == 1L
	)
	# Compute Y
	if (centre_and_scale) {
		Z_means <- apply(Z, 2, mean, na.rm = TRUE)
		Z_sds <- apply(Z, 2, sd, na.rm = TRUE)
		t((t(Z) - Z_means) / Z_sds) # TODO: Print warning if divides by zero
	} else {
		Z
	}
}

#' Create X Matrix
#'
#' @param Y [matrix()] with [numeric()] columns.
#' @param lags [list()] of named integer vectors specifying the lags to use for
#'   each time series in \code{N}.
#'
#' @return [matrix()] X
#' @export
#'
pbsX <- function (Y, lags) {
	# Check arguments
	stopifnot(
		is.matrix(Y),
		is.list(lags),
		all(is.element(names(lags), colnames(Y))),
		length(unique(names(lags))) == length(names(lags)),
		length(unique(colnames(Y))) == length(colnames(Y)),
		is.numeric(as.vector(unlist(Y[, names(lags)]))),
		is.numeric(as.vector(unlist(lags))),
		lags[[1]][1] == 0L
	)
	# Compute X
	lags_size <- unlist(lags, use.names = FALSE)
	lags_name <- rep(names(lags), lengths(lags))
	# Rows in X are points in state space
	X <- pbsLag(Y[, lags_name, drop = FALSE], lags_size)
	colnames(X) <- paste0(lags_name, "_", lags_size)
	class(X) <- unique(c(class(X), "pbsSSR"))
	return(X)
}

#' Create State Space Reconstruction Matrix
#'
#' @param N [matrix()] or [data.frame()] with named [numeric()] columns 
#'   for the response variable and covariate time series.
#' @param lags [list()] of named integer vectors specifying the lags to use for
#'   each time series in \code{N}.
#' @param p [integer()] The forecast distance.
#' @param first_difference [logical()] First-difference each time series?
#' @param centre_and_scale [logical()] Centre and scale each time series?
#'
#' @return [matrix()] State space reconstruction.
#' @export
#'
pbsSSR <- function (N,
										lags,
										p = 1L,
										first_difference = FALSE,
										centre_and_scale = FALSE) {
	# Compute X the state space reconstruction (SSR)
	N <- pbsN(N = N, lags = lags, p = p)
	Z <- pbsZ(N = N, first_difference = first_difference)
	Y <- pbsY(Z = Z, centre_and_scale = centre_and_scale)
	X <- pbsX(Y = Y, lags = lags)
	return(X)
}

#' Compute Distances Between Points in the State Space Reconstruction
#'
#' @param X [matrix()] with named [numeric()] columns.
#' @param lags [list()] of named integer vectors specifying the lags to use for
#'   each time series in \code{N}.
#' @param p [integer()] The forecast distance.
#' @param first_difference [logical()] First-difference each time series?
#'
#' @return
#' @export
#'
pbsDist <- function (X,
										 lags,
										 p = 1L,
										 first_difference = FALSE) {
	# Check arguments
	stopifnot(
		is.matrix(X) & is.numeric(X),
		is.element("pbsSSR", class(X)),
		is.list(lags),
		length(unique(names(lags))) == length(names(lags)),
		is.numeric(as.vector(unlist(lags))),
		lags[[1]][1] == 0L,
		is.integer(p) && length(p) == 1L,
		is.logical(first_difference) && length(first_difference) == 1L
	)
	# Avoid partial-component distances
	X[is.na(rowSums(X)), ] <- NA 
	# Compute distance matrix
	X_distance <- as.matrix(dist(X))
	# Exclude the focal point from each row of neighbours
	diag(X_distance) <- NA
	# Index points in X that contain NAs
	na_rows <- which(is.na(rowSums(X)))
	# Exclude all neighbours of focal points that themselves contain NAs
	X_distance[na_rows, ] <- NA
	# Exclude all neighbours that contain NAs
	X_distance[, na_rows] <- NA
	# Index points that project to points that contain NAs
	na_rows <- (which(is.na(rowSums(X))) - p)
	na_rows <- subset(na_rows, na_rows > 0)
	# Exclude neighbours that project to points that contain NAs
	X_distance[, na_rows] <- NA
	# Index points whose elements include a projection of the focal value
	seq_rows <- seq_len(nrow(X))
	fcl_rows <- rep(seq_rows, each = length(lags[[1]]))
	prj_rows <- fcl_rows + p + lags[[1]]
	na_mat <- matrix(c(fcl_rows, prj_rows), ncol = 2)[which(prj_rows <= nrow(X)),]
	# Exclude neighbours whose elements include a projection of the focal value
	X_distance[na_mat] <- NA # Specify [row, col] pairs
	if (first_difference) {
		# Index points whose elements include a first difference of the focal value
		dif_rows <- fcl_rows - 1 + lags[[1]]
		na_mat <- matrix(c(fcl_rows, dif_rows), ncol = 2)[which(dif_rows <= nrow(X)),]
		# Exclude neighbours whose elements include a first difference of the focal value
		X_distance[na_mat] <- NA # Specify [row, col] pairs
	}
	# Set class
	class(X_distance) <- unique(c(class(X_distance), "pbsDist"))
	return(X_distance)
}

#' Return a Lagged Matrix
#'
#' @param x [matrix()] A vector or column matrix
#' @param n [integer()] The lag sizes
#'
#' @return A matrix or vector. If `x` is a vector then returns a vector `length(x)`
#'   with `NA` for the first `n` values then the first `length(x) - n` values of
#'   `x`.
#'   
#' @export
#'
#' @examples
#' pbsLag(1:10, 3)
#' pbsLag(matrix(rep(1:10, 2), nrow = 10), 3)
#' pbsLag(matrix(rep(1:10, 2), nrow = 10), c(3, 5))
#'
pbsLag <- function (x,
                    n = 1) {

	# Check arguments
	stopifnot(
		is.matrix(x) || is.numeric(x),
		is.numeric(n)
	)
	# Lag x
	m <- as.matrix(x)
	if (length(n) == 1) {n <- rep(n, ncol(m))}
	for (i in seq_along(n)) {
		m[, i] <- c(rep(NA_real_, floor(n[i])), m[, i])[seq_along(m[, i])]
	}
	if (is.vector(x)) {
		m <- as.vector(m)
	} else if (is.matrix(x)) {
		m <- as.matrix(m)
	}
	m
}

# As for sizeSpectra package, need these to avoid warnings related to dplyr
# commands (e.g. referring to the column names within dplyr::filter() ).
# Copied from the warning given by check() (that puts them alphabetical, then
# maybe add some more at the end).
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))
if (getRversion() >= "2.15.1") {
	utils::globalVariables(c(
		"Nx_lags_orig",
		"Xt",
		"Xtmin1",
		"d",
		"my.pred",
		"rEDM.pred"
	))
}
