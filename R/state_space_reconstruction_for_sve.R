#' State Space Reconstruction for Single View Embedding (and then Multiview)
#'
#' Given data of variables, a list of lags, and a defined response variable
#' (TODO: never quite sure why one has to be specified, have asked Luke, think
#' it's slightly inefficient but keeping it for now), create matrix (S.10)
#' from our first manuscript. Absolute numbers are first differenced (for each
#' column, i.e. each original variable), then standardised (centred and scaled),
#' and then the lagged matrix is created.
#'
#' Rewriting from `state_space_reconstruction()` to include all the steps we
#' gave in Appendix S1 of first manuscript.
'
#'
#' @param data [matrix()] or [data.frame()] with variables as named columns
#' @param response [character()] column name of the response variable
#' @param lags [list()] of a named vector of lags for each explanatory variable
#'
#' @author Andrew M. Edward and Luke A. Rogers
#'
#' @return [state_space_reconstruction()] [matrix()] with unlagged response
#'   and lagged explanatory variables that are first differenced and then centred on their means and divided by
#'   their respective standard deviations, with automatically generated column
#'   names in the following style: for a variable called `S_t` and lags of 0, 1,
#'   2, the scaled variable and lagged names are `S_t_s_0, `S_t_s_1, and
#'   `S_t_s_2, where the last number is the lag and `s` stands for scaled. First
#'   column is the scaled response variable.
#'
#' @export
#'
#' @examples
#' d <- data.frame(x = 1:10, y = 11:20)
#' state_space_reconstruction(d, response = "x", lags = list(y = c(0, 1, 2, 3)))
#'
state_space_reconstruction_for_sve <- function(data,
                                               response,
                                               lags){

  # Define values --------------------------------------------------------------
  col_names <- c(response, names(lags))
  lag_sizes <- unlist(lags, use.names = FALSE)
  lag_names <- rep(names(lags), lengths(lags))

  # Create Z -------------------------------------------------------------------
  # Also doing first differencing now, as per our write up.
  data_first_differenced <- rbind(diff(as.matrix(data)),
                                  rep(NA, ncol(data)))

  Z <- as.matrix(data_first_differenced[, col_names, drop = FALSE])
                                        # TODO creates two R_t
                                        # columns, don't think we need those
  Z_means <- apply(Z, 2, mean, na.rm = TRUE)
  Z_sds <- apply(Z, 2, stats::sd, na.rm = TRUE)

  # Create Y -------------------------------------------------------------------

  Y <- t((t(Z) - Z_means) / Z_sds)

  # Create X -------------------------------------------------------------------

  X <- cbind(
    Y[, response, drop = FALSE],
    create_lags(
      Y[, lag_names, drop = FALSE],
      lag_sizes
    )
  )

  # Rename with lags explicitly given, and "_s" to denote scaled (and
  # first-differenced), to distinguish from raw values. So like going from N_{t,k}
  # to Y_{t-lag,k} in our first manuscript Appendix.
  colnames(X) <- c(paste0(response, "_s"), paste0(lag_names, "_s_", lag_sizes))

  return(structure(X, class = "state_space_reconstruction"))
}
