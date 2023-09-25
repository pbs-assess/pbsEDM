#' Create Subset Lags
#'
#' Creates all combinations, though not in quite the logical order one might
#'   expect. Could resolve by reording the code maybe, but it works (Andy
#'   manually checked a detailed example, essentially the one given).
#'
#' @param lags [list()] of a named vector of lags for each explanatory
#'   variable.
#'
#' @return [list()] with one lag [list()] for each subset state space
#'   reconstruction
#' @export
#'
#' @examples
#' create_subset_lags(list(a = c(0, 1, 2), b = c(0, 1)))
#'
create_subset_lags <- function (lags) {
  # Check arguments

  # Create subset lags
  len <- length(unlist(lags))
  num <- 2^len - 1
  out <- list()
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
    out[[i]] <- subs
  }
  # Return a list of subset lags lists
  return(out)
}
