#' Make a tibble from MVE results to then use in a plot, indicating the axes used in the top multiview embedding reconstructions
#'
#' Takes the results list from [multiview_embedding()] and indicates the
#' variables and lags used in the top reconstruction subsets, with the first
#' column being the respective rho value.
#'
#' @param res_mve list object resulting from a call to [multiview_embedding()]
#' @param ... passed onto [create_named_lags()], namely `order` to order the
#' variables consistently.
#' @examples
#' \dontrun{
#'   make_mve_spaces_tibble(age11_res)
#' }
#'
make_mve_spaces_tibble <- function(res_mve,
                                   ...){
  all_lags <- create_named_lags(res_mve$subset_lags,
                                ...)

  num_top <- length(res_mve$lags_of_top_subsets)   # How many 'top' subsets.

  # Create a data frame with rho values in first column and then columns for each lag variable
  df <- data.frame(rho = res_mve$rho_each_top_subset)  # They are ordered by rho

  # Initialize columns for each lag variable with FALSE
  for (lag_name in all_lags) {
    df[[lag_name]] <- FALSE
  }

  # Fill in TRUE/FALSE based on which lags appear in each top subset
  for (i in seq_len(num_top)) {
    subset_lags <- res_mve$lags_of_top_subsets[[i]]
    # Get the named lags for this subset
    subset_named_lags <- create_named_lags(list(subset_lags), ...)
    # Mark TRUE for each lag in this subset
    df[i, subset_named_lags] <- TRUE
  }

  # Convert to tibble
  result <- tibble::as_tibble(df)

  return(result)
}
