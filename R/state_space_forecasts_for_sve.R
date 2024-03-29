#' Forecasts for Single View Embedding Based on Only the Nearest Neighbour, Not Simplex
#'
#' Adapting, simplifying and correcting `state_space_forecasts()` to be called
#'   from `single_view_embedding()`. Based on what's described in Hao and
#'   Sugihara (2016; Science 353:922). This function will do some of the
#'   following steps in one go (basically looping through the t* values).
#'
#' Steps are (basically a modified version of our Table 3 in manuscript 1) TODO
#' add to manuscript, move elsewhere
#' - take a single view embedding (input into here), which is matrix X in
#' manuscript 1, with each row representing time and each column a component
#' (dimension) of the state space
#' - pick focal time `t&*` to make projection from, we need to know `x_{t^*+1}` (see
#' allowable focal values `t^*` in manuscript 1)
#' - define the set of candidate nearest neighbours to `x_t^*`
#' - find the single nearest neighbour of the focal point `x_t^*`
#' - make a projection of `x_{t^*+1}` based on that neighbour (no weighting
#' needed, though we could add that in if necessary)
#' - Repeat for all appropriate focal times `t^*`
#' - Calculate the correlation coefficient between predictions and known
#' observations - TODO make this just based on our response variable or on all
#' dimensions??? First seems sensible, but second seems to be capturing more of
#' the ethos of Takens' theorem.
#' - In Simplex the above would be repeated for different `E`, in multiview we
#' will repeat this for different views (state space reconstructions) in
#' `multiview_embedding()`, which will call the current function.
#'
#' Return the projected values of the ssr, to then be used to calculate the
#' correlation coefficient (or metric of interest) and make the forecast if this
#' particular sve is deemed good.
#'
#' @param ssr_input [matrix()] a state space reconstruction in which the rows
#'   are points in the state space
#' @param distance [matrix()] of distances of allowed neighbours
#'
#' @param lags_of_response_variable NULL or vector of lags being considered for
#'   the response variable, to help eliminate invalid candidate nearest neighbours
#' @param response_s vector of values of scaled response variable, which may or
#'   may not be in the `ssr_input`
#' @author Andrew M. Edwards and Luke A. Rogers
#'
#' @return [vector()] where index of each element corresponds to a `t^*`, and
#'   the value of the element gives the index of the projection (i.e. projection
#'   for `t^*+1`, based on where the nearest neighbour went). Because just
#'   using the single nearest neighbour, can just return indices and deal with
#'   of indices of projections for each `t^* + 1`projected values of t,
#'   including for `T+1`.
#' @export
#'
state_space_forecasts_for_sve <- function(ssr_input,
                                          distance,
                                          lags_of_response_variable,
                                          response_s){
#                                          min_number_neighbours = 1, # minimum
                                            # number of neighbours in order to
                                            # pick the closest, TODO implement


  # Andy adapted code based on
  # the above description. Remove un-needed arguments, just doing the algorithm.
  # TODO may want a minimum number of candidate nearest neighbours to choose
  # from? See distances_small - maybe just require a simple rule. Think about
  # how this relates to Figure 4 of manuscript 1. Though may not be needed.

  # Create neighbour index matrix ----------------------------------------------

  # We just want the t value of the nearest neighbour for each t*, using valid
  # t* times and valid candidate nearest neighbours.

  nearest_neighbour_index <- rep(NA, nrow(distance))

  for(t_star in 1:nrow(distance)){
    # If want to have a minimum number neighbours then do it here: but not
    # needed now it's symmetric and we can have future neighbours, which we do
    # for real data.

    # Work out nearest neighbours for each t_star, but ignore (leave
    #  nearest_neighbour_index[t_star] as NA):
    #  - t_star's for which distance is a row of NA's (distance calculation has
    #  taken care of ssr rows with any NA's; this includes first ones and
    #  row T = nrow(distance)).
    # - Was going to ignore t_star = T-1 because we don't know where it goes (so can't
    #  test how well this ssr performs), but we do need it to make the actual
    #  forecast into the future - calculation of rho will ignore this value
    #  anyway because the known value at T will be NA. Also, may be okay to have
    #  that value when we have no zero lags.
    if(!all(is.na(distance[t_star, ]))){
       # & t_star != nrow(distance) - 1){

      # Invalid candidate nearest neighbours are rows (see first manuscript, page 16):
      # - t* itself (already have NA as distance to itself, so that's taken care
      # of)
      # - x_t that are not fully defined; again, taken care of above as row of
      # NA's in distance
      # - T - 1 as we need to know where it goes to use for prediction; also taken
      # care of in distances (T - 1 column is NA's)  HERE TODO - not true if we
      # have no non-zero lags.
      # - we think we should not use any X_t that includes the quantity we are trying
      # to predict. For univariate this was anything including Y_{t^*+1}; for
      # multivariate it includes any row with Y_{t^*+1, 1}, where 1 is response
      # variable, and then any row
      # that goes to that (so any row with a Y_{t*, 1}. Had thought this was
      # rows t*+1 to t* + max_lag + 1, but that's wrong (considers all
      # variables, but we're only excluding response one, plus assumes
      # continuously incrementing lags (and I think assumes a 0 lag). Now have
      # considered just lags of the response variable that end up at
      # response_{t^*+1} and response_{t^*} which goes to the former. See notes.

      # - ALSO, neighbours t for which response_{t+1} is undefined, which needs
      #   to be made explicit if the response variable (think with 0 lag) is not in the
      #   state space

      # So for t_star we look at the row distance[t_star, ], but need to set
      #  values corresponding to invalid candidate neighbours to NA.
      # which.min will ignores NA's later

      candidate_neighbours_distances <- as.vector(distance[t_star, ])

      # For each lag, l, of the response variable, need to discard t* + l and
      # t* + l + 1 (see written notes, will type them up). lags_of_response_variable is a vector, so can use to index
      if(!is.null(lags_of_response_variable)){
        candidate_neighbours_distances[t_star + lags_of_response_variable + 1] <- NA
        candidate_neighbours_distances[t_star + lags_of_response_variable] <- NA
      }

      # neighbours t for which response_{t+1} is undefined, including T
      response_t_plus_1_undefined <- union(which(is.na(response_s)),
                                           length(response_s))

      # Take off one that would lead to index of 0, this works if empty, though
      # shouldn't be as have now included the final one above
      response_t_plus_1_undefined <-
        response_t_plus_1_undefined[(response_t_plus_1_undefined - 1) > 0.5]


      if(length(response_t_plus_1_undefined) != 0){
        candidate_neighbours_distances[response_t_plus_1_undefined - 1] <- NA
      }

      nearest_neighbour_index[t_star] <- which.min(candidate_neighbours_distances)
    }
  }

  # Then just want to add one to the indices to produce the projections for
  # each one.

  projection_index <- nearest_neighbour_index + 1     # And this will just
                                                      # include the final time step, so a
                                                      # forecast for T+1

  return(projection_index)
}
