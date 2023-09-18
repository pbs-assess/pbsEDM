#' State Space Distance Matrix
#'
#' @details Row index corresponds to focal point time. Column index
#'   corresponds to neighbour point time. The value represents the distance
#'   from the focal point to the neighbour point. Disallowed focal point
#'   and neighbour combinations have value NA.
#'
#' @param ssr [matrix()] a state space reconstruction in which the rows
#'   are points in the state space
#' @param index [integer()] time index of the first value to forecast
#' @param buffer [integer()] number of values to forecast before \code{index}
#'
#' @author Luke A. Rogers
#'
#' @return [matrix()] of allowed neighbour distances, rows corresponding to
#'   focal point time, columns to neighbour point time.
#' @export
#'
#'
state_space_distances <- function (ssr, index = 50L, buffer = 10L) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_integerish(index, lower = 1, upper = nrow(ssr) + 1, len = 1)
  checkmate::assert_integerish(buffer, lower = 1, upper = 10, len = 1)
  checkmate::assert_number(index - buffer, lower = 1)

  # Compute distances ----------------------------------------------------------

  # Avoid partial component distances
  ssr_na <- ssr
  ssr_na[is.na(rowSums(ssr)), ] <- NA_real_

  # Compute the distance matrix (careful if debugging as have 'distances' in vignette)
  distances <- as.matrix(stats::dist(ssr_na))    # Usual Euclidean L2-norm
                                        # (Pythagoras); as.matrix is more intiuitve
  # distances_all <- as.matrix(stats::dist(ssr_na,
  #                                      diag = TRUE,
  #                                      upper = TRUE))  # These options fill in
  #                                      the diagonal as 0's (even those with
  #                                      NA's though) and fills in upper
  #                                      triangle. Slightly more intuitive
  #                                      for understanding, though not needed
  #                                      for calculations. Just has NA's in
  #                                      first three rows and columns.
  # Though think this is what is happening here, except for the 0's maybe.

  # Exclude focal point and future neighbours ----------------------------------
  # This seems to really just change diag to NA not 0.
  distances[upper.tri(distances, diag = TRUE)] <- NA_real_

  # Exclude points by using NA ----------------------------------------



  # - TODO does this deal with our Aspect from first manuscript???

  na_rows <- which(is.na(rowSums(ssr)))                # Rows with NA's
  na_proj <- subset(na_rows - 1L, na_rows - 1L > 0)    # Valid rows that project to na_rows
  # Exclude focal points that contain missing values
  distances[na_rows, ] <- NA_real_      # Think already covered automatically,
                                        # but keep; not sure how distance would
                                        # be calculated.

  # Neighbours that contain missing values
  distances[, na_rows] <- NA_real_      # Again, not sure if they'd be defined.

  # Neighbours that project to points that contain missing values. This seems
  # like it wouldn't always fall out automatically.
  distances[, na_proj] <- NA_real_

  # Exclude focal points in the training set -----------------------------------
  #  Not quite sure about this, but with index = 2 and buffer = 1 (and five dimensions)this just
  #  refers to rows that have NA's. May be more important in other situations.
  # distances[seq_len(index - buffer - 1L), ] <- NA_real_
  distances[1:(index - buffer - 1), ] <- NA_real_

  # Return the distance matrix -------------------------------------------------

  return(distances)
}
