#' Represent An Integer As A Binary Vector
#'
#' @param x [integer()] >= 0
#' @param digits [integer()] length of the output vector
#'
#' @return [integer()] [vector()] binary representation of x
#' @export
#'
#' @examples
#' binary(10, digits = 8)
#'
#' for (i in 0:32) print(binary(i, digits = 6))
#'
binary <- function (x, digits = NULL) {
  # Check arguments

  # Compute binary vector
  v <- c()
  while (x > 0) {
    r <- x %% 2
    x <- x %/% 2
    v <- c(r, v)
  }
  if (!is.null(digits)) {
    if (digits < length(v)) {
      stop("condition digits >= length(v) must hold")
    } else {
      v <- c(rep(0, digits - length(v)), v)
    }
  }
  return(v)
}
