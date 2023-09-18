#' Forecast via Simplex Projection
#'
#' @param N [data.frame()] or [numeric()] Response variable time series as
#'   first column or as a vector.
#' @param E [integer()] Vector of embedding dimensions.
#' @param p [integer()] Scalar forecast distance.
#' @param first_difference [logical()] First-difference the time series
#' @param centre_and_scale [logical()] Centre and scale the time series?
#' @param verbose [logical()] Print progress?
#'
#' @return A list of class \code{pbsSimplex}
#' @author Luke A. Rogers
#' @export
#'
#' @examples
#' N <- data.frame(x = 1:30)
#' s1 <- pbsSimplex(N)
#'
#' N <- data.frame(x = simple_ts)
#' s2 <- pbsSimplex(N)
#'
pbsSimplex <- function (N,
                        E = 1:10,
                        p = 1L,
                        first_difference = FALSE,
                        centre_and_scale = FALSE,
                        verbose = FALSE) {

  # Check arguments ------------------------------------------------------------

  stopifnot(
    is.data.frame(N) | (is.vector(N) & is.numeric(N)),
    is.numeric(E) & floor(E) == E & E > 0L,
    is.integer(p) && length(p) == 1L,
    is.logical(first_difference) && length(first_difference) == 1L,
    is.logical(centre_and_scale) && length(centre_and_scale) == 1L,
    is.logical(verbose) && length(verbose) == 1L
  )

  # Define N -------------------------------------------------------------------

  if (is.numeric(N)) {
    N <- data.frame(Obs = N)
  } else {
    colnames(N)[1] <- "Obs"
  }

  # Compute --------------------------------------------------------------------

  results_list <- list()
  results <- data.frame()

  # TODO: Parallelize this using sockets for compatibility with Windows
  for (i in seq_along(E)) {
    # Define lags
    lags <- list(Obs = seq(0L, E[i] - 1))
    # Store value
    results_list[[i]] <- pbsEDM(N,
                                lags,
                                p,
                                first_difference,
                                centre_and_scale,
                                verbose = FALSE)
    results <- rbind(results, results_list[[i]]$results)
  }

  # Return value ---------------------------------------------------------------

  return(structure(list(
    results_list = results_list,
    results = results),
    class = c("pbsSimplex")))
}
