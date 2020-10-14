

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
	X <- pbsLag(Y[, lags_name], lags_size) # Rows in X are points in state space
	colnames(X) <- paste0(lags_name, "_", lags_size)
	class(X) <- unique(c(class(X), "pbsSSR"))
	return(X)
}

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
#' If `x` is a matrix then TODO [Andy thinks it does the
#'   obvious thing].
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
	if (is.vector(x)) {m <- as.vector(m)}
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
