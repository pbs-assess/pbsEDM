#' Multiview Simplex -- Do Simplex for Permutations of Lags
#'
#' Take a multivariate data set containing variables through time, a response
#' variable (the variable we care about, such as population number), and a set
#' of lags that we wish to consider. Apply Simplex to each subset of lags. Based
#' on `multiview_embedding()` which was based on Ye and Sugihara (2016)'s
#' description of that; not sure if anyone has done this approach before (though
#' seems fairly logical to try it).
#'
#' For all allowed lags, this function builds every possible state space
#' reconstruction (Figure 1C in Ye and Sugihara), of all possible dimensions. So
#' variable 1 with 0 lag, variable 1 with 0 lag and variable 2 with 0 lag,
#' variable 1 with 0 lag and variable 2 with 1 lag, etc. up to variable 1 with
#' all lag allowed from `lags` and variable 2 with all lags allowed from `lags`.
##'  If there are $N$ total
#' variable-lag combinations, then there should be $2^N - 1$ possible
#' reconstructions (each one is either in or out, minus them all being out), but
#' this may get reduced with what is described above.
#'
#' This function then calls `pbsEDM_optimal()`, for each possible state space
#' reconstruction, which creates the state space
#' reconstruction and makes Simplex predictions
#' for all allowable focal times $t^*$ taking into account which neighbours
#' should be candidates for nearest neighbours, and calculates
#' predicted values of the response variable both scaled and unscaled. Here we
#' will calculate metrics of the fit, based on just the response variable (as
#' that is what we are interested in), and pick the best performing ones - the
#' square root of the total number of reconstructions, as per Y&S for multiview embedding. Then use the
#' average of those to make the actual forecast for the time step after the
#' data.
#'
#' TODO Also an add option to just use the best Simplex projection. Might be
#' best to do that as default.
#'
#' Returns **
#'
#' @param data [tibble()] or [data.frame()] with named [numeric()] columns
#' @param lags [list()] of a named vector of lags for each explanatory
#'   variable.
#' @param response [character()] column name of the response variable in
#'   \code{data}
#'
#' @author Andrew M. Edwards
#'
#' @return [list()] **TODO
#' @examples
#' \dontrun{
#' # Stupid example as Y_t is already first differences of N_t; was using to
#'   debug. Best comes out as rho=0.97 for lags of 0 and 0, but then can't do
#'   forecast (as Y_t[100] not defined).
#'  res <- multivariate_simplex(data = NY_lags_example, response = "N_t", lags = list("N_t" = 0:2, "Y_t" = 0:1))
#' }
#' @export
#'
multivariate_simplex <- function(data,
                                 lags,
                                 response = NULL){

  # Move response to first column of data, as required by pbsEDM().
  if(!is.null(response)){
    data <- dplyr::relocate(data,
                            response)
  }

  # Create subset lags ---------------------------------------------------------
  subset_lags <- create_subset_lags(lags,
                                    response_name = response)   # Ignore ones
                                        # that don't have 0 lag of response variable

  num_subsets <- length(subset_lags)

  response_each_subset <- list()
  rho_each_subset <- rep(NA, num_subsets)      # rho based on unscaled (original)
  rho_s_each_subset <- rep(NA, num_subsets)    # rho based on scaled


  # Do simplex for each subset, calc rho both ways
  for(i in 1:num_subsets){
    # print(i)

    # This was response_calc in mve
    simplex_res <- pbsEDM(N = data,
                          lags = subset_lags[[i]],
                          first_difference = TRUE,
                          centre_and_scale = TRUE)


    ## response_calc <- single_view_embedding_for_sve(data = data,
    ##                                                response = response,#"R_t",
    ##                                                lags = subset_lags[[i]])

    # TODO need a check in pbsEDM for correlations, probably
    ## if(all(is.na(response_calc))){        # Because lags are highly correlated in
    ##                                  #  state_space_reconstruction_for_sve();
    ##                                  #  actually will be a single NA but need to
    ##                                  #  check like this
    ##   rho_each_subset[i] <- NA
    ##   rho_s_each_subset[i] <- NA
    ##   response_each_subset[[i]] <- NA    # Likely need a switch later to do with
    ##                                     # this single NA that is not same size tibble
    ##                                     # as non-NA ones
#    } else {

    rho_each_subset[i] <- simplex_res$results$N_rho

    # For mve did this, should be same
    #      cor(response_calc[, response],
    #                          response_calc[, paste0(response,
    #                                                 "_predicted")],
    #                          use = "pairwise.complete.obs")

    rho_s_each_subset[i] <- simplex_res$results$X_rho

    response_each_subset[[i]] <- simplex_res$N_forecast
 #   }
  }

  # Now just pick the top rho from all subsets, and use that to make the
  # forecast. Could later add an option to average the top ones, so leaving
  # those calcs in from multiview_embedding()

  # sqrt_num_subsets <- round(sqrt(num_subsets))  # approx 2^(N/2) I think

  # order_of_subsets_for_rho_index <- order(rho_each_subset, decreasing = TRUE)

  # Indices of the top rho's, e.g. 4229 13073 12825 13104 ...:
  # top_subsets_for_rho_index <-
  #   order_of_subsets_for_rho_index[1:sqrt_num_subsets]

  # The actual top rho's (in order of size)
  # rho_each_top_subset <- rho_each_subset[top_subsets_for_rho_index]

  # This shouldn't be needed; if do then change data
  # response_predicted_from_each_top_subset <- matrix(nrow = nrow(data) + 1,
  #                                              ncol =
  #                                                length(top_subsets_for_rho_index))
                                        # rows are time, columns are
                                        # each top_subsets_for_rho_index in order

  # lags_of_top_subsets <- list()        # List of list of each top lag

  # for(i in 1:length(top_subsets_for_rho_index)){

  #  actual_subset_index <- top_subsets_for_rho_index[i]

  #  response_predicted_from_each_top_subset[, i] <-
  #    dplyr::pull(response_each_subset[[actual_subset_index]], paste0(response,
  #                                                         "_predicted"))
  #  lags_of_top_subsets[[i]] <- subset_lags[[actual_subset_index]]
  #}

  # Take the mean for each t* across all the top subsets
  # response_predicted_from_mve <- rowMeans(response_predicted_from_each_top_subset,
  #                                   na.rm = FALSE)

  # So here we currently just want the single best performing subset and use
  #  that alone, based on rho - TODO check, that may be for all the axes, not
  #  just for response variable. May want to calculate it here based on just the
  #  response. TODO TODO


  top_subset_for_rho_index <- which.max(rho_each_subset)

  rho_top_subset <- max(rho_each_subset)
  lags_of_top_subset <- subset_lags[[top_subset_for_rho_index]]
# Some of response_each_subset have NA in final one, but some don't? That was
# for the example given in help, which wasn't a great one (see above); so code
# was working properly.

  response_predicted_from_top_subset <- response_each_subset[[top_subset_for_rho_index]]

  return(list(
    subset_lags = subset_lags,          # all the combinations of lags
    # rho_top_subset = rho_top_subset,    # rho for top subset, based on
                                        # untransformed response
    # rho_s_top_subset = rho_s_each_subset, # rho_s for each subset
    # response_top_subset = response_predicted_from_top_subset,
    top_subset_for_rho_index = top_subset_for_rho_index, # index of the
                                        # subset with the highest rho
    rho_each_subset = rho_each_subset, # rho for
                                       # the top subsets (original order is
                                       # retained)

    # response_predicted_from_each_top_subset = response_predicted_from_each_top_subset,
    lags_of_top_subset = lags_of_top_subset,
    response_predicted_from_multivariate_simplex = response_predicted_from_top_subset,
    rho_prediction_from_multivariate_simplex = rho_top_subset))
}
