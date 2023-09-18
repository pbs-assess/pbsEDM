#' Forecast via S-Mapping
#'
#' @description Perform short-term nonlinear forecasting via S-mapping..
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
#' @param verbose Logical. Print progress?
#'
#' @details The name of the first element in \code{lags} must match the name of
#'   the response variable in \code{N}. Unlagged time series, including the
#'   response variable, must be specified by a zero in the corresponding named
#'   vector in \code{lags}. For example, given a \code{data.frame} with named
#'   columns \code{Predator}, \code{Prey} and \code{Temperature},
#'   \code{Predator} can be specified as the unlagged response variable by
#'
#'   \code{lags = list(Predator = c(0, ...), ...)}.
#'
#'   This places the unlagged time series of \code{Predator} abundance (or its
#'   optionally first-differenced and/or centred and scaled counterpart) along
#'   the first axis of the reconstructed state space. To predict \code{Predator}
#'   abundance from its first two lags, and from the unlagged and first lags of
#'   \code{Prey} and \code{Temperature}, \code{lags} can be specified as
#'
#'   \code{lags = list(Predator = c(0:2), Prey = c(0:1), Temperature = c(0:1))}.
#'
#'   This example generalizes to arbitrary (possibly non-consecutive) lags of
#'   arbitrarily many covariates (up to limitations of time series length).
#'
#' @return A list of class \code{pbsEDM} containing:
#'
#'   \itemize{
#' \item \code{N} [matrix()] Response variable and unlagged covariates as
#'   columns
#'
#' \item \code{N_observed} [vector()] Response variable time series
#'
#' \item \code{N_forecast} [vector()] Forecast of response variable time series
#'
#' \item \code{X} [matrix()] Unlagged and lagged state variables as columns
#'
#' \item \code{X_observed} [vector()] Transformed response variable time series
#'
#' \item \code{X_forecast} [vector()] Forecast of transformed response variable
#'
#' \item \code{X_distance} [matrix()] Square distance \code{matrix} between
#'   pairs of points in state space (pairs of rows in \code{X})
#'
#' \item \code{neighbour_distance} [matrix()] Distance by focal time (row) and
#'   rank (column)
#'
#' \item \code{neighbour_index} [matrix()]
#'   Neighbour index by focal time (row) and distance rank (column)
#'
#' \item
#'   \code{neighbour_value} [matrix()] Neighbour value by focal time (row) and
#'   distance rank (column)
#'
#' \item \code{neighbour_weight} [matrix()] Neighbour weight by focal time (row)
#'   and distance rank (column)
#'
#' \item
#'   \code{projected_index} [matrix()] Projected neighbour index by projected
#'   time (row) and neighbour distance rank (column)
#'
#' \item
#'   \code{projected_value} [matrix()] Projected neighbour value by projected
#'   time (row) and neighbour distance rank (column)
#'
#' \item
#'   \code{projected_weight} [matrix()] Projected neighbour weight by projected
#'   time (row) and neighbour distance rank (column)
#'
#' \item \code{lags} [list()]
#'   A named list of integer vectors specifying the lags to use for each time
#'   series in \code{N}
#'
#' \item \code{theta} [numeric()] Local weighting parameter
#'
#' \item \code{p} [integer()] The forecast distance
#'
#' \item \code{first_difference} [logical()] First difference each time series?
#'
#' \item \code{centre_and_scale} [logical()] Centre and scale each time series?
#'
#' \item \code{results} [data.frame()] A summary of forecast accuracy
#' }
#'
#'
#' @author Luke A. Rogers
#' @export
#'
#' @examples
#' N <- matrix(rep(1:30, 5), ncol = 5)
#' colnames(N) <- c("A", "B", "C", "D", "E")
#' lags <- list(A = c(0, 1, 2), B = c(0, 1), C = c(0, 1, 2))
#' m1 <- pbsSmap(N, lags)
#'
#' N <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m2 <- pbsSmap(N, lags)
#'
#' N <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m3 <- pbsSmap(N, lags, theta = 2)
#'
#' N <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m4 <- pbsSmap(N, lags, first_difference = TRUE)
#'
#'
#' # And see pbsSmap vignette.
pbsSmap <- function (N,
                     lags = NULL,
                     theta = 0,
                     p = 1L,
                     first_difference = FALSE,
                     centre_and_scale = FALSE,
                     exclusion_radius = "half",
                     verbose = FALSE) {

  tictoc::tic("pbsEDM")

  #----------------- Check arguments ------------------------------------------#

  stopifnot(
    is.matrix(N) || is.data.frame(N),
    is.list(lags),
    all(is.element(names(lags), colnames(N))),
    length(unique(names(lags))) == length(names(lags)),
    length(unique(colnames(N))) == length(colnames(N)),
    is.numeric(as.vector(unlist(N[, names(lags)]))),
    is.numeric(as.vector(unlist(lags))),
    lags[[1]][1] == 0L,
    is.numeric(theta) && length(theta) == 1L,
    is.integer(p) && length(p) == 1L,
    is.logical(first_difference) && length(first_difference) == 1L,
    is.logical(centre_and_scale) && length(centre_and_scale) == 1L,
    is.logical(verbose) && length(verbose) == 1L
  )

  #----------------- Define X -------------------------------------------------#

  if (verbose) cat("\ndefining state space X\n")
  N <- pbsN(N = N, lags = lags, p = p)
  Z <- pbsZ(N = N, first_difference = first_difference)
  Y <- pbsY(Z = Z, centre_and_scale = centre_and_scale)
  X <- pbsX(Y = Y, lags = lags)

  #----------------- Exclude disallowed neighbours ----------------------------#

  if (verbose) cat("defining neighbour distances\n")
  # Distances between points in state space (row vectors in X)
  X_distance <- pbsDist(X,
                        lags,
                        p,
                        first_difference,
                        exclusion_radius = exclusion_radius)

  #----------------- Create neighbour index matrix ----------------------------#
  # TODO: Continue from here (notation and algorithm)

  if (verbose) cat("defining neighbours\n")
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

  if (verbose) cat("projecting nearest neighbours forward by p\n")
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

  if (verbose) cat("computing forecasts of X\n")
  X_forecast <- sapply(X = seq_rows,
                       FUN = function(X, l, m) sum(m[, X] * l[X, ]),
                       m = c_matrix,
                       l = prj_lags)
  X_forecast[is.nan(X_forecast)] <- NA_real_

  #----------------- Compute Z_forecast ---------------------------------------#

  if (centre_and_scale) {
    Z_means <- apply(Z, 2, mean, na.rm = TRUE)
    Z_sds <- apply(Z, 2, stats::sd, na.rm = TRUE)
    Z_forecast <- Z_means[1] + X_forecast * Z_sds[1] # X_forecast == Y_forecast
  } else {
    Z_forecast <- X_forecast
  }

  #----------------- Compute N_forecast ---------------------------------------#

  if (verbose) cat("computing forecasts of N\n")
  if (first_difference) {
    # N_forecast_{t} = N_observed_{t-1} + Z_forecast_{t-1}
    N_forecast <- pbsLag(N_observed) + pbsLag(c(Z_forecast, NA_real_))
  } else {
    N_forecast <- Z_forecast
  }

  # Prepare p-value ------------------------------------------------------------

  X_pval <- smap_surrogates(N,
                            lags,
                            theta,
                            p,
                            first_difference,
                            centre_and_scale)

  #----------------- Prepare return values ------------------------------------#

  N_rho <- stats::cor(N_observed, N_forecast, use = "pairwise.complete.obs")
  N_rmse <- sqrt(mean((N_observed - N_forecast)^2, na.rm = TRUE))
  X_rho <- stats::cor(X_observed, X_forecast, use = "pairwise.complete.obs")
  X_rmse <- sqrt(mean((X_observed - X_forecast)^2, na.rm = TRUE))
  E <- length(lags_size)
  results <- data.frame(E = E,
                        theta = theta,
                        N_rho = N_rho,
                        N_rmse = N_rmse,
                        X_rho = X_rho,
                        X_rmse = X_rmse,
                        X_pval = X_pval,
                        stringsAsFactors = FALSE)

  #----------------- Return a list --------------------------------------------#

  tictoc::toc()

  if (verbose) cat("returning list of class pbsEDM\n")
  structure(
    list(
      N = N,
      N_observed = N_observed,
      N_forecast = N_forecast,
      Z = Z,
      Z_observed = Z[, 1],
      Z_forecast = Z_forecast,
      Y = Y,
      Y_observed = Y[, 1],
      Y_forecast = X_forecast,
      X = X,
      X_observed = X_observed,
      X_forecast = X_forecast,
      X_distance = X_distance,
      neighbour_distance = nbr_dist,
      neighbour_index = nbr_inds,
      neighbour_value = NULL,
      neighbour_weight = nbr_wgts,
      projected_index = prj_inds,
      projected_value = prj_vals,
      projected_weight = prj_wgts,
      lags = lags,
      theta = theta,
      p = as.integer(p),
      first_difference = first_difference,
      centre_and_scale = centre_and_scale,
      exclusion_radius = exclusion_radius,
      results = results
    ),
    class = c("pbsEDM", "pbsSmap")
  )
}
