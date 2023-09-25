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

browser()
  # Do single view embedding for each subset, calc rho both ways
  for(i in 1:num_subsets){
print(i)
# fails at St = 1 being the only lag
    response_calc <- single_view_embedding_for_sve(data = simulated_small_3,
                                                   response = "R_t",
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


  # HERE
  # Take subsets with top sqrt() rho's

  #  use those to get forecast,
  #   taking the mean of all the forecasts
  # Compare with known value outside of this function.

  # HERE
  # Weight redm forecasts ------------------------------------------------------

  # Return results object ------------------------------------------------------

  return(
    structure(
      list(
        data = data,
        observed = c(dplyr::pull(data, response), NA),
        forecast = c(rep(NA_real_, index - 1), weighted$results$forecast),
        response = response,
        lags = lags,
        index = index,
        buffer = buffer,
        window = window,
        metric = metric,
        beyond = beyond,
        n_weight = n_weight,
        raw_forecasts = forecasts,
        ranks = weighted$ranks,
        summary = weighted$summary,
        hindsight = weighted$hindsight,
        results = weighted$results
      ),
      class = "multiview_embedding"
    )
  )
}
