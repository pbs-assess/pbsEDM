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
#' @return [matrix()] of allowed neighbour distances
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

  # Compute the distance matrix
  distances <- as.matrix(stats::dist(ssr_na))

  # Exclude focal point and future neighbours ----------------------------------
# ***TODO*** we don't want to do this either

  distances[upper.tri(distances, diag = TRUE)] <- NA_real_

  # Exclude points with a missing value ----------------------------------------
  # - Neighbours of focal points that themselves contain missing values
  # - Neighbours that contain missing values
  # - Neighbours that project to points that contain missing values
  # - TODO does this deal with our Aspect from first manuscript???
  na_rows <- which(is.na(rowSums(ssr)))
  na_proj <- subset(na_rows - 1L, na_rows - 1L > 0)
  distances[na_rows, ] <- NA_real_
  distances[, na_rows] <- NA_real_
  distances[, na_proj] <- NA_real_

  # Exclude focal points in the training set -----------------------------------

  distances[seq_len(index - buffer - 1L), ] <- NA_real_

  # Return the distance matrix -------------------------------------------------

  return(distances)
}
