#' State Space Reconstruction for Single View Embedding (and then Multiview)
#'
#' Given data of variables and a list of lags,
#' create matrix with rows being time and columns being each dimension in the
#' state space; i.e. (S.10) from our first manuscript. Absolute numbers of the
#' variables are first differenced (for each column, i.e. each original
#' variable), then standardised (centred and scaled), and then the lagged matrix
#' is created.
#'
#' Rewriting from `state_space_reconstruction()` to include all the steps we
#' gave in Appendix S1 of first manuscript, and to not have the response as an
#' output (as it was getting used in the distance calculations later).
#'
#' If this function changes then need to update `untransform_predictions()` also.
#'
#' @param data [matrix()] or [data.frame()] with variables as named columns
#' @param lags [list()] of a named vector of lags for each explanatory variable
#' @param response_only if TRUE then `data` can only have one column and
#'   `lags` a list of the form `lags = list(R_t = 0)` (only the response
#'   variable and must be lag of 0). The variable specified is considered the
#'   response variable and is being rescaled and renamed here, and renamed, for
#'   example, `R_t_s`.
#'
#' @author Andrew M. Edward and Luke A. Rogers
#'
#' @return [state_space_reconstruction()] [matrix()] with lagged variables that
#'   are first differenced and then centred on their means and divided by  their
#'   respective standard deviations, with automatically generated column
#'   names in the following style: for a variable called `S_t` and lags of 0, 1,
#'   2, the scaled variable and lagged names are `S_t_s_0, `S_t_s_1, and
#'   `S_t_s_2, where the last number is the lag and `s` stands for scaled.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' d <- data.frame(x = 1:10, y = 11:20)
#' state_space_reconstruction(d, lags = list(y = c(0, 1, 2, 3)))
#'
#' # If thinking of x as being the response variable, but not having it in the
#' # embedding, then do
#' response_scaled <- state_space_reconstruction_for_sve(d, lags = list(x = 0))
#' names(response_scaled) <- substr(names(response_scaled), 1,
#'   nchar(names(response_scaled)) - 2)   TODO check these
#'
#' }
state_space_reconstruction_for_sve <- function(data,
                                               lags,
                                               response_only = FALSE){

  # Define values --------------------------------------------------------------
  col_names <- names(lags)
  lag_sizes <- unlist(lags, use.names = FALSE)
  lag_names <- rep(names(lags), lengths(lags))

  if(response_only){
    stopifnot(length(lag_sizes) == 1)
    stopifnot(lag_sizes == 0)
  }


  # Create Z -------------------------------------------------------------------
  # Also doing first differencing now, as per our write up.
  data_first_differenced <- rbind(diff(as.matrix(data)),
                                  rep(NA, ncol(data)))

  Z <- as.matrix(data_first_differenced[, col_names, drop = FALSE])
  Z_means <- apply(Z, 2, mean, na.rm = TRUE)
  Z_sds <- apply(Z, 2, stats::sd, na.rm = TRUE)

  # Create Y -------------------------------------------------------------------

  Y <- t((t(Z) - Z_means) / Z_sds)

  # Create X -------------------------------------------------------------------

  X <- create_lags(Y[,
                     lag_names,
                     drop = FALSE],
                   lag_sizes)

  # Rename with lags explicitly given, and "_s" to denote scaled (and
  # first-differenced), to distinguish from raw values. So like going from N_{t,k}
  # to Y_{t-lag,k} in our first manuscript Appendix.
  if(response_only){
    colnames(X) <- paste0(lag_names, "_s")   # e.g just R_t_s
  } else {
    colnames(X) <- paste0(lag_names, "_s_", lag_sizes)
  }

  return(structure(X, class = "state_space_reconstruction"))
}
