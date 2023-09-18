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
