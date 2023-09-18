##' Do the S-map calculations for a vector of `theta` values
##'
##' Given a vector of theta values and data, do the S-map calculation for each theta, using the
##' fast `smap_efficient()`. If values are not already first-differenced then set
##' `first_difference = TRUE` so the differencing is done within
##' `smap_efficient()`.
##'
##' @param N A data frame with named columns for the response variable and
##'   covariate time series.
##' @param lags A list of named integer vectors specifying the lags to use for
##'   each time series in \code{N}.
##' @param theta_vec Vector of theta values to use.
##' @param ... Further options to pass to `smap_efficient()`. In particular
##'   `first_difference` may need to be TRUE, the default is FALSE.
##'
##' @return Vector of values of `rho` corresponding to each value of `theta_vec`.
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'  N <- data.frame(x = simple_ts)
##'  lags <- list(x = 0:1)
##'  res <- smap_thetavec(N, lags, theta_vec = seq(0, 2, by=0.1), first_difference = TRUE)
##'  plot(res)
##' }
##' # And see pbsSmap vignette.
smap_thetavec <- function(N,
                          lags,
                          theta_vec = seq(0, 1, by=0.1),
                          ...){

  rho_vec = rep(NA,
                length(theta_vec))

  for(i in 1:length(theta_vec)){
    rho_vec[i] <- smap_efficient(N,
                                 lags = lags,
                                 theta = theta_vec[i],
                                 ...)
  }

  return(rho_vec)
}
