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
                                          distance){
#                                          min_number_neighbours = 1, # minimum
                                            # number of neighbours in order to
                                            # pick the closest, TODO implement


  # TODO not changed any code yet, just writing the help. Adapt code based on
  # the above description. Remove un-needed arguments, just do the algorithm.
  # TODO may want a minimum number of candidate nearest neighbours to choose
  # from? See distances_small - maybe just require a simple rule. Think about
  # how this relates to Figure 4 of manuscript 1.

  # Create neighbour index matrix ----------------------------------------------
browser()


  #num_nbrs <- ncol(ssr_input) + 1
  #seq_nbrs <- 1:num_nbrs
  # For each t* (row), the nearest neighbour's index goes in first column,
  # second nearest in the second column, up to num_nbrs (drop neighbours we
  # don't care about). (TODO probably set to 1
  # for multiview embedding, we only care about the nearest in each view). TODO
  # That may be wrong (don't think so though), as further down it's taking different components of each
  # point in ssr.
  # Andy confused what is happening here. Think best just to start again, then
  # maybe delete all this.
  #nbr_inds <- t(apply(distance, 1, order))[, seq_nbrs, drop = FALSE]


  # We just want the t value of the nearest neighbour for each t*, using valid
  # t* times.

  nearest_neighbour_index <- rep(NA, nrow(distance))
  for(i in 1:nrow(distance)){
    # If want to have a minimum number neighbours then do it here: but not
    # needed now it's symmetric and we can have future neighbours, which we do
    # for real data.
    if(!all(is.na(distance[i, ]))){  #ignore rows of NA's
      nearest_neighbour_index[i] <- which.min(distance[i, ])
    }
  }

  # Jump to end

  # Add a row of NA's, presumably to get filled in
  # nbr_inds <- rbind(nbr_inds, array(NA, dim = c(1, num_nbrs)))

  # Create neighbour matrices --------------------------------------------------
  # TODO Confused, this next bit seems wrong. For simple example I have row 10 having in
  #  nbr_inds of     8    6    4    5    9    7. So t=8 is the nearest
  # neighbour, t = 6 the next one, etc.
  # This seems to then take R_t (first dimension of ssr, ssr_input here) for t=8, but
  # then next element (S_t_0 in my example, second column of ssr_input) for t = 6, etc.
  #  So mixing up definitions, as you shouldn't link R_8 with S_6_0 etc.
  #  And if just want the actual nearest neighbour and value of R_t we can use
  # that. Even TODO select it based on the best solution of R_t{*+1} not the
  # full state space.
#  nbr_vals <- t(apply(nbr_inds, 1, function (x, y) y[x, 1], y = ssr_input))

  # TODO Not sure what next bit is, but then does weights which is simplex, not
  # done according to Hao and Sugihara (2016). Andy to write simple function
  # based on this, doing what H&S 2016 describe.
  # AHA - this is taken from pbsEDM() and only been slightly changed.
#  nbr_dist <- t(apply(distance, 1, sort, na.last = TRUE))[, seq_nbrs]
#  nbr_wts <- t(apply(nbr_dist, 1, function (x) exp(-x / x[1])))
#  nbr_wts <- rbind(nbr_wts, array(NA, dim = c(1, num_nbrs)))

  # Project neighbour matrices -------------------------------------------------

#  proj_inds <- create_lags(nbr_inds, 1L) + 1L
#  proj_vals <- t(apply(proj_inds, 1, function (x, y) y[x, 1], y = ssr_input))
#  proj_wts <- create_lags(nbr_wts, 1L)

  # Compute ssr_input_forecast ---------------------------------------------------------

#  ssr_input_forecast <- as.vector(rowSums(proj_vals * proj_wts) / rowSums(proj_wts))

  # Forecast beyond ssr? -------------------------------------------------------

#  if (!beyond) {
#    ssr_input_forecast <- ssr_input_forecast[seq_len(nrow(ssr_input))]
#  }

  # Return ssr_input_forecast ----------------------------------------------------------

  #  return(ssr_input_forecast)



  # TODO exclude invalid t*'s

  # TODO exclude invalid candidate nearest neighbours - likely needs to be done
  #  in earlier loop, will depend on t* given our Figure 4.

  # Then just want to add one to the indices to produce the projections for
  # each one.

  projection_index <- nearest_neighbour_index + 1     # And this will just
                                                      # include the final time step, so a
                                                      # forecast for T+1

  return(projection_index)
}
