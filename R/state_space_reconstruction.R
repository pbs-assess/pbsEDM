#' State Space Reconstruction (SSR)
#'
#' @description Rows are centred and scaled points in the state-space
#'   reconstruction.
#'
#' @param data [matrix()] with variables as named columns
#' @param response [character()] column name of the response variable
#' @param lags [list()] of a named vector of lags for each explanatory variable
#'
#' @author Luke A. Rogers
#'
#' @return [state_space_reconstruction()] [matrix()] with unlagged response
#'   and lagged explanatory variables centred on their means and scaled by
#'   their respective standard deviations, with automatically generated column
#'   names.
#'
#' @export
#'
#' @examples
#' d <- data.frame(x = 1:10, y = 11:20)
#' state_space_reconstruction(d, response = "x", lags = list(y = c(0, 1, 2, 3)))
#'
state_space_reconstruction <- function (data, response, lags) {

  # Check arguments ------------------------------------------------------------


  # Define values --------------------------------------------------------------

  col_names <- c(response, names(lags))
  lag_sizes <- unlist(lags, use.names = FALSE)
  lag_names <- rep(names(lags), lengths(lags))

  # Create Z -------------------------------------------------------------------

  Z <- as.matrix(data[, col_names, drop = FALSE])
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

  # Return ssr -----------------------------------------------------------------

  return(structure(X, class = "state_space_reconstruction"))
}
