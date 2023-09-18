#' SMap That Returns Only The Forecast Accuracy
#'
#' @param N A data frame with named columns for the response variable and
#'   covariate time series.
#' @param lags A list of named integer vectors specifying the lags to use for
#'   each time series in \code{N}.
#' @param theta The numeric local weighting parameter.
#' @param p The integer forecast distance.
#' @param first_difference Logical. First-difference each time series?
#' @param centre_and_scale Logical. Centre and scale each time series?
#' @param exclusion_radius Number of points around ${\bf x}_t^*$ to exclude as
#'   candidate nearest neighbours; either default of `half` as used for our manuscript
#'   (see equation (6)), or a number to match
#'   the `exclusionRadius` setting in `rEDM::Simplex()`. See `?pbsDist` for more details.
#'
#' @return [numeric()] The forecast accuracy rho
#' @export
#'
#' @examples
##' \donttest{
##'  N <- data.frame(x = simple_ts)
##'  lags <- list(x = 0:1)
##'  m4 <- smap_efficient(N, lags, first_difference = TRUE)
##' }
##' # And see pbsSmap vignette.
smap_efficient <- function (N,
                            lags = NULL,
                            theta = 0,
                            p = 1L,
                            first_difference = FALSE,
                            centre_and_scale = FALSE,
                            exclusion_radius = "half"){

  #----------------- Define X -------------------------------------------------#

  N <- pbsN(N = N, lags = lags, p = p)
  Z <- pbsZ(N = N, first_difference = first_difference)
  Y <- pbsY(Z = Z, centre_and_scale = centre_and_scale)
  X <- pbsX(Y = Y, lags = lags)

  #----------------- Exclude disallowed neighbours ----------------------------#

  # Distances between points in state space (row vectors in X)
  X_distance <- pbsDist(X,
                        lags,
                        p,
                        first_difference,
                        exclusion_radius = exclusion_radius)

  #----------------- Create neighbour index matrix ----------------------------#
  # TODO: Continue from here (notation and algorithm)

  nbr_dist <- t(apply(X_distance, 1, sort, na.last = TRUE))
  nbr_inds <- t(apply(X_distance, 1, order))
  nbr_inds[which(is.na(nbr_dist))] <- NA
  # nbr_vals <- t(apply(nbr_inds, 1, function(x, y) y[x, 1], y = X))
  nbr_wgts <- t(apply(nbr_dist,
                      1,
                      function(x, y) exp(-y * x / mean(x, na.rm = TRUE)),
                      y = theta))

  #----------------- Compute lag of neighbour index matrix --------------------#

  # TODO: Needed?
  lag_inds <- pbsLag(nbr_inds, p)

  #----------------- Project neighbour matrices -------------------------------#

  prj_inds <- pbsLag(nbr_inds, p) + p
  prj_vals <- t(apply(prj_inds,
                      1,
                      function(x, y) y[x, 1],
                      y = X))
  prj_wgts <- pbsLag(nbr_wgts, p)

  #----------------- Project xt_lag matrix ------------------------------------#

  prj_lags <- pbsLag(X, p)

  #----------------- Compute B matrix for SVD ---------------------------------#

  # The row gives the focal index
  # The col gives the nearest neighbours ordered relative to focal index
  b_matrix <- prj_wgts * prj_vals
  b_matrix[which(is.na(b_matrix))] <- 0

  #----------------- Compute W array of matrices for SVD ----------------------#

  # The row (first dimension) gives the nearest neighbours relative to focal
  # The (second dimension) gives the X row vector index
  # The col (third dimension) gives the focal index
  seq_rows <- seq_len(nrow(X))
  lags_size <- unlist(lags, use.names = FALSE)
  w_array <- sapply(X = seq_rows,
                    FUN = function(X, w, y) w[X, ] %*% t(rep(1, y)),
                    w = prj_wgts,
                    y = length(lags_size),
                    simplify = "array")

  #----------------- Compute L array of lagged row vectors for SVD ------------#

  # The row (first dimension) gives the nearest neighbours relative to focal
  # The (second dimension) gives the X row vector index
  # The col (third dimension) gives the focal index
  l_array <- sapply(X = seq_rows,
                    FUN = function(X, l, m) l[m[X, ], ],
                    l = X,
                    m = lag_inds, # Double check
                    simplify = "array")

  #----------------- Compute A array of matrices for SVD ----------------------#

  # The row (first dimension) gives the nearest neighbours relative to focal
  # The (second dimension) gives the X row vector index
  # The col (third dimension) gives the focal index
  a_array <- w_array * l_array
  a_array[which(is.na(a_array))] <- 0

  #----------------- Solve for C matrix via SVD -------------------------------#

  # Decompose A matrices by SVD
  svd_list <- apply(a_array, 3, svd)

  # Simplify
  vdu_array <- sapply(X = seq_rows,
                      FUN = function(X, s) s[[X]]$v %*% diag(1/s[[X]]$d) %*%
                        t(s[[X]]$u),
                      s = svd_list,
                      simplify = "array")

  # Solve for C matrix
  c_matrix <- sapply(X = seq_rows,
                     FUN = function(X, a, b) a[,, X] %*% b[X, ],
                     a = vdu_array,
                     b = b_matrix)

  #----------------- Store observed as vectors --------------------------------#

  X_observed <- X[, 1]
  N_observed <- N[, 1]  # Check that response var. is in first column

  #----------------- Compute X_forecast ---------------------------------------#

  X_forecast <- sapply(X = seq_rows,
                       FUN = function(X, l, m) sum(m[, X] * l[X, ]),
                       m = c_matrix,
                       l = prj_lags)
  X_forecast[is.nan(X_forecast)] <- NA_real_

  #----------------- Prepare return values ------------------------------------#

  X_rho <- stats::cor(X_observed, X_forecast, use = "pairwise.complete.obs")

  # Return rho -----------------------------------------------------------------

  return(X_rho)
}
