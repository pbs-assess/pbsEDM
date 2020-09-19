#' Forecast via Empirical Dynamic Modelling
#'
#' @param N [data.frame()] A data frame with named columns for the response 
#' variable and covariate time series
#' @param lags [list()] A named list of integer vectors specifying the lags to
#' use for each time series in \code{N}
#' @param p [integer()] The forecast distance
#' @param first_difference [logical()] First difference each time series?
#' @param centre_and_scale [logical()] Centre and scale each time series?
#'
#' @details The name of the first element in \code{lags} must match the name of 
#' the response variable in \code{N}. Unlagged time series, including the 
#' response variable, must be specified by a zero in the corresponding named
#' vector in \code{lags}. For example, given a \code{data.frame} with named
#' columns \code{Predator}, \code{Prey} and \code{Temperature}, \code{Predator} 
#' can be specified as the unlagged response variable by 
#' 
#' \code{lags = list(Predator = c(0, ...), ...)}. 
#' 
#' This places the unlagged time series of \code{Predator} abundance along the
#' first axis of the reconstructed state space. To predict \code{Predator}
#' abundance from its first two lags, and from the unlagged and first lags of
#' \code{Prey} and \code{Temperature}, \code{lags} can be specified as
#' 
#' \code{lags = list(Predator = c(0:2), Prey = c(0:1), Temperature = c(0:1))}. 
#' 
#'
#' @return A list of class \code{pbsEDM} containing:
#'
#' \itemize{
#'   \item \code{N} [matrix()] Response variable and unlagged covariates as columns
#'   \item \code{N_observed} [vector()] Response variable time series
#'   \item \code{N_forecast} [vector()] Forecast of response variable time series
#'   \item \code{X} [matrix()] Unlagged and lagged state variables as columns
#'   \item \code{X_observed} [vector()] Transformed response variable time series
#'   \item \code{X_forecast} [vector()] Forecast of transformed response variable
#'   \item \code{X_distance} [matrix()] Square distance \code{matrix} between pairs of 
#'   points in state space (pairs of rows in \code{X})
#'   \item \code{neighbour_distance} [matrix()] Distance by focal time (row) and rank
#'   (column)
#'   \item \code{neighbour_index} [matrix()] Neighbour index by focal time (row) and 
#'   distance rank (column)
#'   \item \code{neighbour_value} [matrix()] Neighbour value by focal time (row) and 
#'   distance rank (column)
#'   \item \code{neighbour_weight} [matrix()] Neighbour weight by focal time (row) and 
#'   distance rank (column)
#'   \item \code{projected_index} [matrix()] Projected neighbour index by projected 
#'   time (row) and neighbour distance rank (column)
#'   \item \code{projected_value} [matrix()] Projected neighbour value by projected 
#'   time (row) and neighbour distance rank (column)
#'   \item \code{projected_weight} [matrix()] Projected neighbour weight by projected 
#'   time (row) and neighbour distance rank (column)
#'   \item \code{lags} [list()] A named list of integer vectors specifying the lags to
#'   use for each time series in \code{N}
#'   \item \code{p} [integer()] The forecast distance
#'   \item \code{first_difference} [logical()] First difference each time series?
#'   \item \code{centre_and_scale} [logical()] Centre and scale each time series?
#'   \item \code{results} [data.frame()] A summary of forecast accuracy
#' }
#'
#' @author Luke A. Rogers
#' @export
#'
#' @examples
#' N <- matrix(rep(1:30, 5), ncol = 5)
#' colnames(N) <- c("A", "B", "C", "D", "E")
#' lags <- list(A = c(0, 1, 2), B = c(0, 1), C = c(0, 1, 2))
#' m1 <- pbsEDM(N, lags)
#'
#' N <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m2 <- pbsEDM(N, lags)
#'
#' N <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m3 <- pbsEDM(N, lags, first_difference = TRUE)
#'
pbsEDM <- function (N,
                    lags,
                    p = 1L,
                    first_difference = FALSE,
                    centre_and_scale = FALSE) {

  #----------------- Check arguments ------------------------------------------#

  stopifnot(
    is.numeric(N) || is.matrix(N) || is.data.frame(N),
    is.list(lags),
    is.integer(p) && length(p) == 1L,
    is.logical(first_difference) && length(first_difference) == 1L,
    is.logical(centre_and_scale) && length(centre_and_scale) == 1L
  )

  #----------------- Define nt_observed ---------------------------------------#

  #	nt_observed <- pbsLag(as.vector(N[, names(lags)[1]]), lags[[1]][1])
  # ANDY commenting out since breaks Travis (vignette). Think we just want it
  # like I have in the next if statement.

  #----------------- Transform time series ------------------------------------#

  xt <- as.matrix(N[, names(lags)])
  colnames(xt) <- names(lags)
  # First difference and buffer with NAs
  if (first_difference) {
    nt_observed = xt           # original data before first differencing
    xt <- rbind(apply(xt, 2, diff), NA_real_)
  } else {
    nt_observed = NULL
  }

  if (centre_and_scale) {
    xt_means <- apply(xt, 2, mean, na.rm = TRUE)
    xt_sds <- apply(xt, 2, sd, na.rm = TRUE)
    xt <- t((t(xt) - xt_means) / xt_sds)
  }

  #----------------- Create lagged matrix -------------------------------------#

  # X is a matrix of named lagged column vectors
  lags_size <- unlist(lags, use.names = FALSE)
  lags_name <- rep(names(lags), lengths(lags))
  X <- pbsLag(xt[, lags_name], lags_size)
  colnames(X) <- paste0(lags_name, "_", lags_size)

  #----------------- Create distance matrix -----------------------------------#

  # xt_dist[i, j] is the distance between row vectors i and j in X
  xt_lags_na <- X
  xt_lags_na[which(is.na(rowSums(xt_lags_na))), ] <- NA # For distance
  xt_dist <- as.matrix(dist(xt_lags_na))

  #----------------- Exclude elements from the distance matrix ----------------#

  # Exclude the X focal row vector
  diag(xt_dist) <- NA

  # Exclude X row vectors that contain NAs
  na_rows <- which(is.na(rowSums(X)))
  xt_dist[na_rows, ] <- NA
  xt_dist[, na_rows] <- NA

  # Exclude X rows that project beyond X
  seq_rows <- seq_len(nrow(X))
  na_rows <- which((seq_rows + p) > max(seq_rows))
  xt_dist[na_rows, ] <- NA
  xt_dist[, na_rows] <- NA

  # Exclude X rows that project to rows that contain NAs
  prj_rows <- (which(is.na(rowSums(X))) - p)
  na_rows <- prj_rows[which(prj_rows > 0)]
  # xt_dist[na_rows, ] <- NA # Commented to allow forecast from here
  xt_dist[, na_rows] <- NA

  # Exclude X rows that contain a projection of the focal value
  rep_rows <- rep(seq_rows, each = length(lags[[1]]))
  prj_rows <- rep_rows + p + lags[[1]]
  na_mat <- matrix(c(rep_rows, prj_rows), ncol = 2)[which(prj_rows <= nrow(X)),]
  xt_dist[na_mat] <- NA # Specify [row, col] pairs

  #----------------- Create neighbour index matrix ----------------------------#

  # nbr_inds is an nrow(X) x num_nbrs matrix of X row indices
  num_nbrs <- length(lags_size) + 1
  seq_nbrs <- seq_len(num_nbrs)
  nbr_inds <- t(apply(xt_dist, 1, order))[, seq_nbrs]
  nbr_inds[which(rowSums(!is.na(xt_dist)) < num_nbrs), ] <- NA

  #----------------- Create neighbour matrices --------------------------------#

  # nbr_vals is a matrix of values from X[, 1] corresponding to nbr_inds
  nbr_vals <- t(apply(nbr_inds, 1, function(x, y) y[x, 1], y = X))
  nbr_dist <- t(apply(xt_dist, 1, sort, na.last = T))[, seq_nbrs]
  nbr_wgts <- t(apply(nbr_dist, 1, function(x) exp(-x / x[1])))

  #----------------- Project neighbour matrices -------------------------------#

  prj_inds <- pbsLag(nbr_inds, p) + p
  prj_vals <- t(apply(prj_inds,
                      1,
                      function(x, y) y[x, 1],
                      y = X))
  prj_wgts <- pbsLag(nbr_wgts, p)

  #----------------- Prepare return values ------------------------------------#

  X_observed <- X[, 1]
  X_forecast <- as.vector(rowSums(prj_vals * prj_wgts) / rowSums(prj_wgts))
  rho <- cor(X_observed, X_forecast, use = "pairwise.complete.obs")
  rmse <- sqrt(mean((X_observed - X_forecast)^2, na.rm = TRUE))
  E <- length(lags_size)
  results <- data.frame(E = E,
                        rho = rho,
                        rmse = rmse,
                        stringsAsFactors = FALSE)

  #----------------- Return a list --------------------------------------------#

  structure(
    list(
      N = N,
      N_observed = N[,1], # Check that response var. is in first column
      N_forecast = NULL, # TODO: recover raw forecast from centred & scaled
      X = X,
      X_observed = X_observed,
      X_forecast = X_forecast,
      X_distance = xt_dist,
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
##' Only works for one variable, non-differenced N(t), for now, as I haven't thought how to incorporate
##' the lags variable properly yet.
##'
##' @return List of `pbsEDM` lists, each main component corresponds to a value
##'   of `E`, given by `results$E`
##' @export
##' @author Andrew Edwards
##' @param Nt [vector] Vector of non-differenced values N(t), with time assumed
##'   to be `1:length(Nt)`.
##' @param E_vec The vector of embedding dimensions to try.
##' @param ... Further options to pass to `pbsEDM()`.
##' @examples
##' \donttest{
##'   aa <- pbsEDM_Evec(Nx_lags_orig$Nt)
##' }
pbsEDM_Evec <- function(Nt,
                     E_vec = 1:10,
                     ...){
  E_res <- list()
  for(i in 1:length(E_vec)){
    E_res[[i]] <- pbsEDM(data.frame(Nt = Nt),
                         lags = list(Nt = 0:E_vec[i]),
                         first_difference = TRUE,
                         ...)
  }
  return(E_res)    # could make it class pbsEDM_Evec to automate plot.pbsEDM_Evec
}


#' Forecast via S-Mapping
#' 
#' @param N [data.frame()] A data frame with named columns for the response 
#' variable and covariate time series
#' @param lags [list()] A named list of integer vectors specifying the lags to
#' use for each time series in \code{N}
#' @param theta [numeric()] Local weighting parameter 
#' @param p [integer()] The forecast distance
#' @param first_difference [logical()] First difference each time series?
#' @param centre_and_scale [logical()] Centre and scale each time series?
#'
#' @details The name of the first element in \code{lags} must match the name of 
#' the response variable in \code{N}. Unlagged time series, including the 
#' response variable, must be specified by a zero in the corresponding named
#' vector in \code{lags}. For example, given a \code{data.frame} with named
#' columns \code{Predator}, \code{Prey} and \code{Temperature}, \code{Predator} 
#' can be specified as the unlagged response variable by 
#' 
#' \code{lags = list(Predator = c(0, ...), ...)}. 
#' 
#' This places the unlagged time series of \code{Predator} abundance along the
#' first axis of the reconstructed state space. To predict \code{Predator}
#' abundance from its first two lags, and from the unlagged and first lags of
#' \code{Prey} and \code{Temperature}, \code{lags} can be specified as
#' 
#' \code{lags = list(Predator = c(0:2), Prey = c(0:1), Temperature = c(0:1))}. 
#'
#' @return A list of class \code{pbsEDM} containing:
#'
#' \itemize{
#'   \item \code{N} [matrix()] Response variable and unlagged covariates as columns
#'   \item \code{N_observed} [vector()] Response variable time series
#'   \item \code{N_forecast} [vector()] Forecast of response variable time series
#'   \item \code{X} [matrix()] Unlagged and lagged state variables as columns
#'   \item \code{X_observed} [vector()] Transformed response variable time series
#'   \item \code{X_forecast} [vector()] Forecast of transformed response variable
#'   \item \code{X_distance} [matrix()] Square distance \code{matrix} between pairs of 
#'   points in state space (pairs of rows in \code{X})
#'   \item \code{neighbour_distance} [matrix()] Distance by focal time (row) and rank
#'   (column)
#'   \item \code{neighbour_index} [matrix()] Neighbour index by focal time (row) and 
#'   distance rank (column)
#'   \item \code{neighbour_value} [matrix()] Neighbour value by focal time (row) and 
#'   distance rank (column)
#'   \item \code{neighbour_weight} [matrix()] Neighbour weight by focal time (row) and 
#'   distance rank (column)
#'   \item \code{projected_index} [matrix()] Projected neighbour index by projected 
#'   time (row) and neighbour distance rank (column)
#'   \item \code{projected_value} [matrix()] Projected neighbour value by projected 
#'   time (row) and neighbour distance rank (column)
#'   \item \code{projected_weight} [matrix()] Projected neighbour weight by projected 
#'   time (row) and neighbour distance rank (column)
#'   \item \code{lags} [list()] A named list of integer vectors specifying the lags to
#'   use for each time series in \code{N}
#'   \item \code{theta} [numeric()] Local weighting parameter 
#'   \item \code{p} [integer()] The forecast distance
#'   \item \code{first_difference} [logical()] First difference each time series?
#'   \item \code{centre_and_scale} [logical()] Centre and scale each time series?
#'   \item \code{results} [data.frame()] A summary of forecast accuracy
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
#' m1 <- pbsEDM(N, lags)
#'
#' N <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' m2 <- pbsSMAP(N, lags)
#'
pbsSmap <- function (N,
                     lags,
                     theta = 0,
                     p = 1L,
                     first_difference = FALSE,
                     centre_and_scale = FALSE) {

  #----------------- Check arguments ------------------------------------------#

  stopifnot(
    is.numeric(N) || is.matrix(N) || is.data.frame(N),
    is.list(lags),
    is.numeric(theta) && length(theta) == 1L,
    is.integer(p) && length(p) == 1L,
    is.logical(first_difference) && length(first_difference) == 1L,
    is.logical(centre_and_scale) && length(centre_and_scale) == 1L
  )

  #----------------- Define nt_observed ---------------------------------------#

  nt_observed <- pbsLag(as.vector(N[, names(lags)[1]]), lags[[1]][1])

  #----------------- Transform time series ------------------------------------#

  xt <- as.matrix(N[, names(lags)])
  colnames(xt) <- names(lags)
  # First difference and buffer with NAs
  if (first_difference) {
    xt <- rbind(apply(xt, 2, diff), NA_real_)
  }
  # Centre and scale
  if (centre_and_scale) {
    xt_means <- apply(xt, 2, mean, na.rm = TRUE)
    xt_sds <- apply(xt, 2, sd, na.rm = TRUE)
    xt <- t((t(xt) - xt_means) / xt_sds)
  }

  #----------------- Create lagged matrix -------------------------------------#

  # X is a matrix of named lagged column vectors
  lags_size <- unlist(lags, use.names = FALSE)
  lags_name <- rep(names(lags), lengths(lags))
  X <- pbsLag(xt[, lags_name], lags_size)
  colnames(X) <- paste0(lags_name, "_", lags_size)

  #----------------- Create distance matrix -----------------------------------#

  # xt_dist[i, j] is the distance between row vectors i and j in X
  xt_lags_na <- X
  xt_lags_na[which(is.na(rowSums(xt_lags_na))), ] <- NA # For distance
  xt_dist <- as.matrix(dist(xt_lags_na))

  #----------------- Exclude elements from the distance matrix ----------------#

  # Exclude the X focal row vector
  diag(xt_dist) <- NA

  # Exclude X row vectors that contain NAs
  na_rows <- which(is.na(rowSums(X)))
  xt_dist[na_rows, ] <- NA
  xt_dist[, na_rows] <- NA

  # Exclude X rows that project beyond X
  seq_rows <- seq_len(nrow(X))
  na_rows <- which((seq_rows + p) > max(seq_rows))
  xt_dist[na_rows, ] <- NA
  xt_dist[, na_rows] <- NA

  # Exclude X rows that project to rows that contain NAs
  prj_rows <- (which(is.na(rowSums(X))) - p)
  na_rows <- prj_rows[which(prj_rows > 0)]
  # xt_dist[na_rows, ] <- NA
  xt_dist[, na_rows] <- NA

  # Exclude X rows that contain a projection of the focal value
  rep_rows <- rep(seq_rows, each = length(lags[[1]]))
  prj_rows <- rep_rows + p + lags[[1]]
  na_mat <- matrix(c(rep_rows, prj_rows), ncol = 2)[which(prj_rows <= nrow(X)),]
  xt_dist[na_mat] <- NA # Specify [row, col] pairs

  #----------------- Create neighbour index matrix ----------------------------#

  nbr_dist <- t(apply(xt_dist, 1, sort, na.last = TRUE))
  nbr_inds <- t(apply(xt_dist, 1, order))
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

  #----------------- Make forecasts -------------------------------------------#

  X_forecast <- sapply(X = seq_rows,
                        FUN = function(X, l, m) sum(m[, X] * l[X, ]),
                        m = c_matrix,
                        l = prj_lags)

  X_forecast[is.nan(X_forecast)] <- NA_real_

  #----------------- Prepare return values ------------------------------------#

  X_observed <- X[, 1]
  xt_rho <- cor(X_observed, X_forecast, use = "pairwise.complete.obs")
  xt_rmse <- sqrt(mean((X_observed - X_forecast)^2, na.rm = TRUE))
  xt_dim <- length(lags_size)
  xt_theta <- theta
  E <- xt_dim
  xt_results <- data.frame(xt_dim = xt_dim,
                           E = E,
                           xt_theta = xt_theta,
                           xt_rho = xt_rho,
                           xt_rmse = xt_rmse,
                           stringsAsFactors = FALSE)

  #----------------- Return a list --------------------------------------------#

  structure(
    list(
      N = N,
      N_observed = N[,1], # Check that response var. is in first column
      N_forecast = NULL, # TODO: recover raw forecast from centred & scaled
      X = X,
      X_observed = X_observed,
      X_forecast = X_forecast,
      X_distance = xt_dist,
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
      results = xt_results
    ),
    class = "pbsEDM"
  )
}
