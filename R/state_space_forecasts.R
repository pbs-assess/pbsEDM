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
  seq_nbrs <- 1:num_nbrs
  # For each t* (row), the nearest neighbour's index goes in first column,
  # second nearest in the second column, up to num_nbrs (drop neighbours we
  # don't care about). (TODO probably set to 1
  # for multiview embedding, we only care about the nearest in each view). TODO
  # That may be wrong (don't think so though), as further down it's taking different components of each
  # point in ssr.
  nbr_inds <- t(apply(distance, 1, order))[, seq_nbrs, drop = FALSE]

  # Count number of NA's in each row, if not enough neighbours then no valid
  # neighbours for that focal point
  nbr_inds[which(rowSums(!is.na(distance)) < num_nbrs), ] <- NA

  # Add a row of NA's, presumably to get filled in
  nbr_inds <- rbind(nbr_inds, array(NA, dim = c(1, num_nbrs)))

  # Create neighbour matrices --------------------------------------------------
  # TODO Confused, this next bit seems wrong. For simple example I have row 10 having in
  #  nbr_inds of     8    6    4    5    9    7. So t=8 is the nearest
  # neighbour, t = 6 the next one, etc.
  # This seems to then take R_t (first dimension of ssr, X here) for t=8, but
  # then next element (S_t_0 in my example, second column of X) for t = 6, etc.
  #  So mixing up definitions, as you shouldn't link R_8 with S_6_0 etc.
  #  And if just want the actual nearest neighbour and value of R_t we can use
  # that. Even TODO select it based on the best solution of R_t{*+1} not the
  # full state space.
  nbr_vals <- t(apply(nbr_inds, 1, function (x, y) y[x, 1], y = X))

  # TODO Not sure what next bit is, but then does weights which is simples, not
  # done according to Hao and Sugihara (2016). Andy to write simple function
  # based on this, doing what H&S 2016 describe.
  # AHA - this is taken from pbsEDM() and only been slightly changed.
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
