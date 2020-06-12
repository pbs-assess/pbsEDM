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
#' pbsLAG(1:10, 3)
#' pbsLAG(matrix(rep(1:10, 2), nrow = 10), 3)
#' pbsLAG(matrix(rep(1:10, 2), nrow = 10), c(3, 5))
#'
pbsLAG <- function (x,
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
