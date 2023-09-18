#' Empirical Dynamic Modeling Forecasts
#'
#' @param X [matrix()] a state space reconstruction in which the rows
#'   are points in the state space
#' @param distance [matrix()] of allowed neighbour distances
#' @param beyond [logical()]
#'
#' @author Luke A. Rogers
#'
#' @return [numeric()] [vector()] of forecast values
#' @export
#'
state_space_forecasts <- function (X, distance, beyond = FALSE) {

  # Check arguments ------------------------------------------------------------

  # Create neighbour index matrix ----------------------------------------------

  num_nbrs <- ncol(X) + 1
  seq_nbrs <- seq_len(num_nbrs)
  nbr_inds <- t(apply(distance, 1, order))[, seq_nbrs, drop = FALSE]
  nbr_inds[which(rowSums(!is.na(distance)) < num_nbrs), ] <- NA
  nbr_inds <- rbind(nbr_inds, array(NA, dim = c(1, num_nbrs)))

  # Create neighbour matrices --------------------------------------------------

  nbr_vals <- t(apply(nbr_inds, 1, function (x, y) y[x, 1], y = X))
  nbr_dist <- t(apply(distance, 1, sort, na.last = TRUE))[, seq_nbrs]
  nbr_wts <- t(apply(nbr_dist, 1, function (x) exp(-x / x[1])))
  nbr_wts <- rbind(nbr_wts, array(NA, dim = c(1, num_nbrs)))

  # Project neighbour matrices -------------------------------------------------

  proj_inds <- create_lags(nbr_inds, 1L) + 1L
  proj_vals <- t(apply(proj_inds, 1, function (x, y) y[x, 1], y = X))
  proj_wts <- create_lags(nbr_wts, 1L)

  # Compute X_forecast ---------------------------------------------------------

  X_forecast <- as.vector(rowSums(proj_vals * proj_wts) / rowSums(proj_wts))

  # Forecast beyond ssr? -------------------------------------------------------

  if (!beyond) {
    X_forecast <- X_forecast[seq_len(nrow(X))]
  }

  # Return X_forecast ----------------------------------------------------------

  return(X_forecast)
}
