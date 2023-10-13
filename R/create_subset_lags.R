#' Create Subset Lags
#'
#' Creates all combinations, though not in quite the logical order one might
#'   expect. Could resolve by reording the code maybe, but it works (Andy
#'   manually checked a detailed example, essentially the one given).
#'
#' @param lags [list()] of a named vector of lags for each explanatory
#'   variable.
#' @param response character string; if not NULL then any subset that does not
#'   include `response` with a lag of 0 is excluded; this if for running
#'   `pbsEDM` (which requires response variable in the state space) in
#'   `multivariate_simplex()`.
#'
#' @return [list()] with one lag [list()] for each subset state space
#'   reconstruction
#' @export
#'
#' @examples
#' s <- create_subset_lags(list(a = c(0, 1, 2), b = c(0, 1)))
#' s_with_response <- create_subset_lags(list(a = c(0, 1, 2), b = c(0, 1)), response = "b")
create_subset_lags <- function (lags,
                                response_name = NULL) {
  # Check arguments

  # Create subset lags
  len <- length(unlist(lags))
  num <- 2^len - 1
  subsets_list <- list()
  includes_response <- logical(num)  # default is FALSE

  for (i in seq_len(num)) {
    inds <- as.logical(binary(i, digits = len))
    subs <- unlist(utils::as.relistable(lags))
    subs[!rev(inds)] <- NA_real_
    subs <- utils::relist(subs)
    subs <- lapply(subs, function(x) x[!is.na(x)])
    for (j in rev(seq_along(subs))) {
      if (length(subs[[j]]) == 0) {
        subs[[j]] <- NULL
      }
    }
    subsets_list[[i]] <- subs

    # Which subsets include the response with lag of 0
    if(!is.null(response_name)){
      if(response_name %in% names(subs)){
        if(0 %in% subs[[response_name]]){
          includes_response[i] <- TRUE
        }
      }
    }
  }

  if(!is.null(response_name)){
    out <- subsets_list[includes_response]     # only include ones with response lag 0;
                                               # out is list but this works
  } else {
    out <- subsets_list
  }

  # Return a list of subset lags lists
  return(out)
}
