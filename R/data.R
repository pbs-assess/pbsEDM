## Descriptions of saved data objects.

##' Simple example time series
##'
##' Example simulation results originally from Carrie's salmonTraj() function,
##' that end up some discrepancies in results between `rEDM` and our own EDM
##' code. Contains first difference values, $X(t) = N(t+1) - N(t)$ where $N(t)$
##' are population numbers from the from a simple Ricker-like model. `X` is not
##' standardised (which is okay as univariate time series).
##' @format Vector of length 99 (since original `N(t)` was for 100, with each
##' element representing $X(t) = N(t+1) - N(t)$.
##'
##' @source Generated from running `data-raw/simple_ts.R`.
"simple_ts"
