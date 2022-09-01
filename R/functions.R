#' Forecast via Empirical Dynamic Modelling
#'
#' @description Perform short-term nonlinear forecasting via Empirical Dynamic
#'   Modelling.
#'
#' @param N A data frame with named columns for the response variable and
#'   covariate time series.
#' @param lags A list of named integer vectors specifying the lags to use for
#'   each time series in \code{N}.
#' @param p The integer forecast distance.
#' @param first_difference Logical. First-difference each time series?
#' @param centre_and_scale Logical. Centre and scale each time series?
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
#' \item \code{neighbour_index} [matrix()] Neighbour index by focal time (row)
#'   and distance rank (column)
#'
#' \item \code{neighbour_value} [matrix()] Neighbour value by focal time (row)
#'   and distance rank (column)
#'
#' \item \code{neighbour_weight} [matrix()] Neighbour weight by focal time (row)
#'   and distance rank (column)
#'
#' \item \code{projected_index} [matrix()] Projected neighbour index by
#'   projected time (row) and neighbour distance rank (column)
#'
#' \item \code{projected_value} [matrix()] Projected neighbour value by
#'   projected time (row) and neighbour distance rank (column)
#'
#' \item \code{projected_weight} [matrix()] Projected neighbour weight by
#'   projected time (row) and neighbour distance rank (column)
#'
#'
#' \item \code{lags} [list()] A named list of integer vectors specifying the
#'   lags to use for each time  series in \code{N}
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
#' @author Luke A. Rogers
#' @export
#'
#' @examples
#' N <- matrix(rep(1:30, 5), ncol = 5)
#' colnames(N) <- c("A", "B", "C", "D", "E")
#' lags <- list(A = c(0, 1, 2), B = c(0, 1), C = c(0, 1, 2))
#' m1 <- pbsEDM(N, lags, verbose = TRUE)
#'
#' N <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m2 <- pbsEDM(N, lags, verbose = TRUE)
#'
#' N <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m3 <- pbsEDM(N, lags, first_difference = TRUE, verbose = TRUE)
#'
pbsEDM <- function (N,
                    lags,
                    p = 1L,
                    first_difference = FALSE,
                    centre_and_scale = FALSE,
                    verbose = FALSE) {

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
  X_distance <- pbsDist(X, lags, p, first_difference)

  #----------------- Create neighbour index matrix ----------------------------#
  # TODO: Continue from here (notation and algorithm)

  if (verbose) cat("defining nearest neighbours\n")
  # nbr_inds is an nrow(X) x num_nbrs matrix of X row indices
  lags_size <- unlist(lags, use.names = FALSE)
  lags_name <- rep(names(lags), lengths(lags))
  num_nbrs <- length(lags_size) + 1
  seq_nbrs <- seq_len(num_nbrs)
  nbr_inds <- t(apply(X_distance, 1, order))[, seq_nbrs]
  nbr_inds[which(rowSums(!is.na(X_distance)) < num_nbrs), ] <- NA

  #----------------- Create neighbour matrices --------------------------------#

  # nbr_vals is a matrix of values from X[, 1] corresponding to nbr_inds
  nbr_vals <- t(apply(nbr_inds, 1, function(x, y) y[x, 1], y = X))
  nbr_dist <- t(apply(X_distance, 1, sort, na.last = T))[, seq_nbrs]
  nbr_wgts <- t(apply(nbr_dist, 1, function(x) exp(-x / x[1])))

  #----------------- Project neighbour matrices -------------------------------#

  if (verbose) cat("projecting nearest neighbours forward by p\n")
  prj_inds <- pbsLag(nbr_inds, p) + p
  prj_vals <- t(apply(prj_inds,
                      1,
                      function(x, y) y[x, 1],
                      y = X))
  prj_wgts <- pbsLag(nbr_wgts, p)

  #----------------- Store observed as vectors --------------------------------#

  X_observed <- X[, 1]
  N_observed <- N[, 1]  # Check that response var. is in first column

  #----------------- Compute X_forecast ---------------------------------------#

  if (verbose) cat("computing forecasts of X\n")
  X_forecast <- as.vector(rowSums(prj_vals * prj_wgts) / rowSums(prj_wgts))

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

  #----------------- Prepare return values ------------------------------------#

  N_rho <- stats::cor(N_observed, N_forecast, use = "pairwise.complete.obs")
  N_rmse <- sqrt(mean((N_observed - N_forecast)^2, na.rm = TRUE))
  X_rho <- stats::cor(X_observed, X_forecast, use = "pairwise.complete.obs")
  X_rmse <- sqrt(mean((X_observed - X_forecast)^2, na.rm = TRUE))
  E <- length(lags_size)
  results <- data.frame(E = E,
                        N_rho = N_rho,
                        N_rmse = N_rmse,
                        X_rho = X_rho,
                        X_rmse = X_rmse,
                        stringsAsFactors = FALSE)

  #----------------- Return a list of class pbsEDM ----------------------------#

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
      neighbour_value = nbr_vals,
      neighbour_weight = nbr_wgts,
      projected_index = prj_inds,
      projected_value = prj_vals,
      projected_weight = prj_wgts,
      lags = lags,
      p = as.integer(p),
      first_difference = first_difference,
      centre_and_scale = centre_and_scale,
      results = results
    ),
    class = "pbsEDM"
  )
}

##' Do the `pbsEDM()` calculation for vector of `E` values
##'
##' Only works for one variable, non-differenced $N_t$, for now, as I haven't thought how to incorporate
##' the lags variable properly yet.
##'
##' @param N_t [vector] Vector of non-differenced values $N_t$, with time assumed
##'   to be `1:length(N_t)`.
##' @param E_vec The vector of embedding dimensions to try. ACTUALLY lags, see
##'   Issue 28.
##' @param ... Further options to pass to `pbsEDM()`.
##' @return List of `pbsEDM` lists, each main component corresponds to a value
##'   of `E`, given by `results$E`
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'   aa <- pbsEDM_Evec(NY_lags_example$N_t)
##' }
pbsEDM_Evec <- function(N_t,
                        E_vec = 1:10,
                        ...){
  E_res <- list()
  for(i in 1:length(E_vec)){
    E_res[[i]] <- pbsEDM(data.frame(N_t = N_t),
                         lags = list(N_t = 0:E_vec[i]),
                         first_difference = TRUE,
                         ...)
  }
  return(E_res)    # could make it class pbsEDM_Evec to automate plot.pbsEDM_Evec
}


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
  X_distance <- pbsDist(X, lags, p, first_difference)

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
      results = results
    ),
    class = c("pbsEDM", "pbsSmap")
  )
}

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
                            centre_and_scale = FALSE){

  #----------------- Define X -------------------------------------------------#

  N <- pbsN(N = N, lags = lags, p = p)
  Z <- pbsZ(N = N, first_difference = first_difference)
  Y <- pbsY(Z = Z, centre_and_scale = centre_and_scale)
  X <- pbsX(Y = Y, lags = lags)

  #----------------- Exclude disallowed neighbours ----------------------------#

  # Distances between points in state space (row vectors in X)
  X_distance <- pbsDist(X, lags, p, first_difference)

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


#' Perform Surrogate Test For Nonlinearity
#'
#' @param N A data frame with named columns for the response variable and
#'   covariate time series.
#' @param lags A list of named integer vectors specifying the lags to use for
#'   each time series in \code{N}.
#' @param theta The numeric local weighting parameter.
#' @param p The integer forecast distance.
#' @param first_difference Logical. First-difference each time series?
#' @param centre_and_scale Logical. Centre and scale each time series?
#'
#' @return [numeric()] Quantile for empirical delta rho among surrogates
#' @export
#' # See pbsSmap vignette for example use.
#'
smap_surrogates <- function (N,
                             lags,
                             theta,
                             p = 1L,
                             first_difference = FALSE,
                             centre_and_scale = FALSE){

  # Compute the empirical delta rho --------------------------------------------

  # Rho at theta
  rho_at_theta <- smap_efficient(N,
                                 lags,
                                 theta,
                                 p,
                                 first_difference,
                                 centre_and_scale)
  # Rho at zero
  rho_at_zero <- smap_efficient(N,
                                lags,
                                0L,
                                p,
                                first_difference,
                                centre_and_scale)
  # Delta rho
  empirical_delta_rho <- rho_at_theta - rho_at_zero

  # Compute the surrogate delta rho --------------------------------------------

  # Number of iterations
  n_surrogates <- 100L
  # Initialize
  surrogate_delta_rho <- numeric(n_surrogates)
  # Iterate
  for (i in seq_len(n_surrogates)) {
    # Permute N
    row_permute <- c(sample(seq_len(nrow(N) - 1L)), nrow(N))
    P <- N[row_permute, , drop = FALSE]
    # Rho at theta
    rho_at_theta <- smap_efficient(P,
                                   lags,
                                   theta,
                                   p,
                                   first_difference,
                                   centre_and_scale)
    # Rho at zero
    rho_at_zero <- smap_efficient(P,
                                  lags,
                                  0L,
                                  p,
                                  first_difference,
                                  centre_and_scale)
    # Delta rho
    surrogate_delta_rho[i] <- rho_at_theta - rho_at_zero
  }

  # Compute the quantile -------------------------------------------------------

  p_val <- 1 - stats::ecdf(surrogate_delta_rho)(empirical_delta_rho)
  if (p_val == 0) p_val <- NA_real_

  # Return the p-value ---------------------------------------------------------

  return(p_val)
}

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
