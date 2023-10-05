#' Multiview Embedding
#'
#' Take a multivariate data set containing variables through time, a response
#' variable (the variable we care about, such as population number), and a set
#' of lags that we wish to consider. Based on Ye and Sugihara (2016)'s
#' description.
#'
#' For all allowed lags, this function builds every possible state space
#' reconstruction (Figure 1C in Ye and Sugihara), of all possible dimensions. So
#' variable 1 with 0 lag, variable 1 with 0 lag and variable 2 with 0 lag,
#' variable 1 with 0 lag and variable 2 with 1 lag, etc. up to variable 1 with
#' all lag allowed from `lags` and variable 2 with all lags allowed from `lags`.
#' TODO I feel that some combinations should be duplicated, in the sense that
#' variable 1 with lag 1 and variable 2 with lag 2 is the same as variable 1
#' with lag 0 and variable 2 with lag 1 (i.e. if you have nothing with lag of 0
#' then you can shift everything), but this will reduce the data a little --
#' look into as may end up keeping state space reconstructions that are
#' essentially the same, just shifted in time. If there are $N$ total
#' variable-lag combinations, then there should be $2^N - 1$ possible
#' reconstructions (each one is either in or out, minus them all being out), but
#' this may get reduced with what is described above.
#'
#' This function calls `single_view_embedding_for_sve()` (the `for_sve` was just
#' to distinguish from earlier functions), which creates the state space
#' reconstruction, calculates the distances between points, makes predictions
#' for all allowable focal times $t^*$ taking into account which neighbours
#' should be candidates for nearest neighbours - predictions are just from the
#' single nearest neighbour (as per Y&S, rather than full Simplex), calculate
#' predicted values of the response variable both scaled and unscaled. Here we
#' will calculate metrics of the fit, based on just the response variable (as
#' that is what we are interested in), and pick the best performing ones - the
#' square root of the total number of reconstructions, as per Y&S. Then use the
#' average of those to make the actual forecast for the time step after the data.
#'
#' Returns **
#'
#' @param data [tibble()] or [data.frame()] with named [numeric()] columns
#' @param response [character()] column name of the response variable in
#'   \code{data}
#' @param lags [list()] of a named vector of lags for each explanatory
#'   variable.
#'
#' @author Andrew M. Edwards and Luke A. Rogers
#'
#' @return [list()] **TODO
#' @export
#'
multiview_embedding <- function(data,
                                response,
                                lags){

  # Create subset lags ---------------------------------------------------------
  subset_lags <- create_subset_lags(lags)

  num_subsets <- length(subset_lags)

  response_each_subset <- list()
  rho_each_subset <- rep(NA, num_subsets)      # rho based on unscaled (original)
  rho_s_each_subset <- rep(NA, num_subsets)    # rho based on scaled

  # Do single view embedding for each subset, calc rho both ways
  for(i in 1:num_subsets){
    # print(i)

    response_calc <- single_view_embedding_for_sve(data = data,
                                                   response = response,#"R_t",
                                                   lags = subset_lags[[i]])

    rho_each_subset[i] <- cor(response_calc[, response],
                              response_calc[, paste0(response,
                                                     "_predicted")],
                              use = "pairwise.complete.obs")

    rho_s_each_subset[i] <- cor(response_calc[, paste0(response, "_s")],
                                response_calc[, paste0(response,
                                                       "_s_predicted")],
                                use = "pairwise.complete.obs")

    response_each_subset[[i]] <- response_calc
  }

  sqrt_num_subsets <- round(sqrt(num_subsets))  # approx 2^(N/2) I think

  order_of_subsets_for_rho_index <- order(rho_each_subset, decreasing = TRUE)
  # Gives, in order, the index of the 1st, 2nd, 3rd, ... highest rho's.
  # > which.max(rho_each_subset)
  #  4229
  # > order_of_subsets_for_rho_index[1]
  #  4229

# browser()
  # Index of the subsets with the highest sqrt_num_subsets rho values
  # Think ties (unlikely given accuracy) will just get ignored and the first one
  # taken. But still use length of this not sqrt_num_subsets in
     # loop below
  # Indices of the top rho's, e.g. 4229 13073 12825 13104 ...:
  top_subsets_for_rho_index <-
    order_of_subsets_for_rho_index[1:sqrt_num_subsets]

 # The actual top rho's (in order of size)
  rho_each_top_subset <- rho_each_subset[top_subsets_for_rho_index]

  R_t_predicted_from_each_top_subset <- matrix(nrow = nrow(response_calc),
                                               ncol =
                                                 length(top_subsets_for_rho_index))
                                        # rows are time, columns are
                                        # each top_subsets_for_rho_index in order

  lags_of_top_subsets <- list()        # List of list of each top lag

  for(i in 1:length(top_subsets_for_rho_index)){

    actual_subset_index <- top_subsets_for_rho_index[i]

    R_t_predicted_from_each_top_subset[, i] <-
      dplyr::pull(response_each_subset[[actual_subset_index]], paste0(response,
                                                           "_predicted"))
    lags_of_top_subsets[[i]] <- subset_lags[[actual_subset_index]]
  }

  # Take the mean for each t* across all the top subsets
  R_t_predicted_from_mve <- rowMeans(R_t_predicted_from_each_top_subset,
                                     na.rm = FALSE)

  # Take response from the last response_calc (all the same), not data
  #  as want forecast value (to be an NA)
  rho_prediction_from_mve <- cor(dplyr::pull(response_calc,
                                             response),
                                 R_t_predicted_from_mve,
                                 use = "pairwise.complete.obs")

  # want to return:
  return(list(
    subset_lags = subset_lags,          # all the combinations of lags
    rho_each_subset = rho_each_subset,  # rho for each subset
    rho_s_each_subset = rho_s_each_subset, # rho_s for each subset
    response_each_subset = response_each_subset,
    top_subsets_for_rho_index = top_subsets_for_rho_index, # indices of the
                                        # subsets with the highest
                                        # sqrt(total number of subsets) rho
                                        # values
    rho_each_top_subset = rho_each_top_subset, # rho for
                                        # the top subsets (original order is retained)
    R_t_predicted_from_each_top_subset = R_t_predicted_from_each_top_subset,
    lags_of_top_subsets = lags_of_top_subsets,
    R_t_predicted_from_mve = R_t_predicted_from_mve,  # R_t predicted by taking
                                        # mean; will contain NaN for ones we
                                        # can't predict, and has one for
                                        # forecasting next time step.
    rho_prediction_from_mve = rho_prediction_from_mve)) # rho from comparing
                                        # R_t_predicted_from_mve with original data
}
