#' State Space Distance Matrix for Single View Embedding
#'
#' Given a state-space reconstruction, return the distances between focal points
#'   (rows) and each candidate nearest neighbour, having ruled out non-candidate
#'   potential neighbours. Adapted from `state_space_distances()`.
#'
#'
#' @param ssr [matrix()] a state space reconstruction in which the rows
#'   are points in the state space, each row is a time
#'
#' @author Luke A. Rogers
#'
#' @return [matrix()] of allowed neighbour distances, rows corresponding to
#'   focal point time, columns to neighbour point time. The value represents the
#'   distance from the focal point to the neighbour point. Disallowed focal
#'   point and neighbour combinations have value NA. Original
#'   `state_space_distances()` had an upper triangular matrix to exclude future
#'   points, but not doing that here -- want to utilise the full state-space
#'   reconstruction.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From vignette:
#' simulated <- EDMsimulate::salmon_sim()
#' ssr_5 <- state_space_reconstruction_for_sve(data = simulated[1:15, ], lags = list(S_t = 1))
#' # TODO - correc this after correcting function, need to better deal with
#'   neighbours that don't hav a projection - hadn't thought of possibility of
#'   being the final row.
#'
#' # So no lag of 0, so final row is complete, so distances come out with only
#'   NA's in first column and first row, not final row like when we have a lag
#'   of 0 --- may be wrong once correct it).
#' ssr_5
#' So just the S_t_s_1 column, which is correct.
#' 'distances_5 <- state_space_distances_for_sve(ssr_5)
#' }
#'
state_space_distances_for_sve <- function(ssr){

  # Check arguments ------------------------------------------------------------

  # checkmate::assert_integerish(index, lower = 1, upper = nrow(ssr) + 1, len = 1)
  # checkmate::assert_integerish(buffer, lower = 1, upper = 10, len = 1)
  # checkmate::assert_number(index - buffer, lower = 1)

  # Compute distances ----------------------------------------------------------

  # Avoid partial component distances (any row with an NA just becomes all NA's,
  # since the point does not exist)
  ssr_na <- ssr
  ssr_na[is.na(rowSums(ssr)), ] <- NA

  # Compute the distance matrix (careful if debugging as have 'distances' in
  # vignette), which has dimension nrow(ssr_na) * nrow(ssr_na). Columns are
  # time points also. This automatically makes the first C columns all NA's if the
  # first C rows are all NA's (except diagonal is 0, this gets fixed next).

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

  # Exclude focal point (i.e. distance from itself)
  diag(distances) <- NA

  # This does not deal with our Aspect 2 from first manuscript because
  #  are just calculating the distances here, not prescribing a specific focal point.

  na_rows <- which(is.na(rowSums(ssr)))         # Incomplete rows of ssr

  # Exclude focal points that contain missing values
  distances[na_rows, ] <- NA      # Think already covered automatically,
                                  # but keep just in case; not sure how distance would
                                  # be calculated.

  # Neighbours that contain missing values
  distances[, na_rows] <- NA      # Again, not sure if they'd be defined.

  # Neighbours that project to points that contain missing values. This seems
  # like it wouldn't always fall out automatically.
  na_proj <- subset(na_rows - 1,
                    na_rows - 1 > 0)            # Valid rows that project to na_rows

  distances[, na_proj] <- NA

  # Neigbhours that project to points that are undefined, as is the case for the
  # final time points since we do not know where it goes, therefore it should
  # not be a neighbour.
  distances[, nrow(distances)] <- NA

  return(distances)
}
