##' Do the `pbsEDM()` calculation for vector of `E` values and then return the
##' optimal `E`.
##'
##' Only works for one variable, non-differenced $N_t$, for now, as I haven't thought how to incorporate
##' the lags variable properly yet. Calls `pbsEDM_Evec` and then just returns
##' the results for the optimal `E`, based on highest `rho`.
##'
##' @param N_t [vector] Vector of non-differenced values $N_t$, with time assumed
##'   to be `1:length(N_t)`.
##' @param lags_to_test vec Vector of lags to try (gets passed onto `E_vec` in
##'   `pbsEDM_Evec`; see Issue 28).
##' @param ... Further options to pass to `pbsEDM_Evec()` and then onto `pbsEDM()`.
##' @return A `pbsEDM` list corresponding to the optimal `E`, based on `N_rho`.
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'   aa <- pbsEDM_optimal_E(NY_lags_example$N_t)
##' }
pbsEDM_optimal_E <- function(N_t,
                             lags_to_test = 1:10,
                             ...){

  res_all_E <- pbsEDM_Evec(N_t,
                           E_vec = lags_to_test,
                           ...)

  N_rho_vals <- vector(length = length(lags_to_test))

  for(i in 1:length(lags_to_test)){
    N_rho_vals[i]  <- res_all_E[[i]]$results$N_rho
  }

  index_of_max_N_rho <- which.max(N_rho_vals)

  return(res_all_E[[index_of_max_N_rho]])
}
