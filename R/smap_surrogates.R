#' Perform Surrogate Test For Nonlinearity
#'
#' See pbsSmap vignette for example use.
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
#'
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
