#' State Space Reconstruction for Single View Embedding (and then Multiview)
#'
#' Given data of variables, a list of lags, and a defined response variable
#' (TODO: never quite sure why one has to be specified, have asked Luke), create matrix (S.10)
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
#'   and lagged explanatory variables that are first differenced and then centred on their means and scaled by
#'   their respective standard deviations, with automatically generated column
#'   names.
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
browser()
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
  colnames(X) <- c(response, paste0(lag_names, "_", lag_sizes))
# TODO here we can add s in for scaled, as they're not the same as the original data
  # Return ssr -----------------------------------------------------------------

  return(structure(X, class = "state_space_reconstruction"))
}
