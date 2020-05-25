#' Return a Lagged Matrix
#'
#' @param x [matrix()] A vector or column matrix
#' @param n [integer()] The lag sizes
#'
#' @return A matrix or vector
#' @export
#'
#' @examples 
#' pbsLAG(1:10, 3)
#' pbsLAG(matrix(rep(1:10, 2), nrow = 10), c(3, 5))
#' 
pbsLAG <- function (x, n) {
	m <- as.matrix(x)
	for (i in seq_along(n)) {
		m[, i] <- c(rep(NA_real_, floor(n[i])), m[, i])[seq_along(m[, i])]
	}
	if (is.vector(x)) {m <- as.vector(m)}
	m
}