#' Create Lags Of A Vector Or Matrix
#'
#' @param x [vector()] or column [matrix()] to lag
#' @param n [integer()] lag sizes
#'
#' @return [vector()] or column [matrix()] of lagged values
#'
#' @author Luke A. Rogers
#'
#' @export
#'
#' @examples
#' create_lags(1:10, 3)
#' create_lags(matrix(rep(1:10, 2), nrow = 10), 3)
#' create_lags(matrix(rep(1:10, 2), nrow = 10), c(3, 5))
#'
#' create_lags(matrix(rep(1:10, 2), nrow = 10), c(0, 1))
#' create_lags(matrix(rep(1:10, 2), nrow = 10), c(0, 0))
#' create_lags(matrix(rep(1:10, 2), nrow = 10), c(0, -1))
#'
create_lags <- function (x, n = 1L) {

  # 0.0 Check arguments --------------------------------------------------------

  stopifnot(
    is.matrix(x) || is.numeric(x),
    is.numeric(n))

  # 1.0 Define m and n ---------------------------------------------------------

  # Coerce to matrix
  m <- as.matrix(x)

  # Define lags
  if (length(n) == 1) {
    n <- rep(n, ncol(m))
  }

  # 2.0 Create lags of m -------------------------------------------------------

  # Create positive or negative lags and buffer by NAs
  for (i in seq_along(n)) {
    if (n[i] >= 0) {
      m[, i] <- c(rep(NA_real_, floor(n[i])), m[, i])[seq_along(m[, i])]
    } else {
      m[, i] <- c(m[, i], rep(NA_real_, floor(-n[i])))[seq_along(m[, i]) - n[i]]
    }
  }

  # 3.0 Coerce m to vector or matrix -------------------------------------------

  if (is.vector(x)) {
    m <- as.vector(m)
  } else if (is.matrix(x)) {
    m <- as.matrix(m)
  }

  # 4.0 Return m ---------------------------------------------------------------

  return(m)
}
