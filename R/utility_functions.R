#' Create a Tibble of Lagged Copies of a Named Tibble Column
#'
#' @param tbl A tibble with a named column to be lagged
#' @param col_name The name of a column in tbl
#' @param embed_dim An integer embedding dimension corresponding to the number
#'     of lagged columns to create
#' @param lag_size The number of time steps separating successive lags
#'
#' @return A tibble of successively lagged columns
#' 
#' @importFrom magrittr %>%
#' @export
#'
#' @examples pbs_make_lags(tibble::tibble(x = 1:15), "x", 5, 1)
#' 
pbs_make_lags <- function(tbl,
                          col_name,
                          embed_dim,
                          lag_size) {
  # Check arguments
  assertthat::assert_that(
    tibble::is_tibble(tbl),
    col_name %in% names(tbl),
    assertthat::is.count(embed_dim),
    assertthat::is.count(lag_size)
  )
  
  # Create time_series vector and ts_index tibble
  time_series <- dplyr::pull(tbl, col_name)
  ts_index <- tibble::tibble(index = seq_along(time_series))
  
  # Return lagged tibble
  lapply(
    X = seq_len(embed_dim) - 1,
    FUN = function(X, ts, n) dplyr::lag(ts, X * n),
    ts = time_series,
    n = lag_size
  ) %>%
    dplyr::bind_cols(ts_index) %>%
    tidyr::drop_na() %>%
    dplyr::right_join(ts_index,
                      by = c("index")) %>%
    dplyr::select(-index) %>%
    magrittr::set_colnames(paste0("lag_", as.character(seq_len(embed_dim) - 1)))
}

#' Calculate Distances Among Tibble Row Vectors
#'
#' @param row_tbl A tibble of row vectors
#'
#' @return A long tibble of distances between pairs of row vectors
#' 
#' @importFrom magrittr %>%
#' @export
#'
#' @examples pbs_calc_dist(pbs_make_lags(tibble::tibble(x = 1:100), "x", 3, 1))
#' 
pbs_calc_dist <- function(row_tbl) {
  
  # Check argument
  assertthat::assert_that(tibble::is_tibble(row_tbl))
  
  # Calculate distances among row vectors
  row_tbl %>%
    as.matrix() %>%
    stats::dist(diag = FALSE, upper = TRUE) %>%
    broom::tidy() %>%
    dplyr::rename(focal_ind = item2, nbr_ind = item1) %>%
    dplyr::select(focal_ind, nbr_ind, distance) %>%
    tidyr::drop_na()
}


#' Attempt to Return a Tibble with a Named Column
#'
#' @param data An object to be returned as a tibble (a tibble, vector with 
#'     null dimensions, or object coercable by tibble::as_tibble)
#' @param col_name A name to identify or apply to a column (character)
#'
#' @return A tibble with a named column
#' 
#' @importFrom magrittr %>%
#' @export
#'
#' @examples pbs_make_tibble(1:10, "x")
#' 
pbs_make_tibble <- function(data, col_name) {
  
  # Is data a tibble with column col_name?
  if (tibble::is_tibble(data)) {
    assertthat::assert_that(col_name %in% names(data))
    data
  } else if (is.vector(data)) {
    assertthat::assert_that(is.character(col_name))
    tibble::tibble(data) %>% magrittr::set_colnames(col_name)
  } else {
    assertthat::assert_that(col_name %in% names(data))
    tibble::as_tibble(data)
  }
}

#' Make a Tibble of Nearest Neighbours
#'
#' @param lag_dist Distances between neighbours (tibble with columns:
#'     'focal_ind', 'nbr_ind', and 'distance')
#' @param time_ind The time index of the focal vector (numeric scalar)
#' @param rel_lib_ind The allowable neighbour indices (numeric vector)
#' @param embed_dim The embedding dimension (numeric scalar)
#' @param pred_dist The prediction distance (numeric scalar)
#' @param max_nbrs The maximum number of neighbours (numeric scalar). It
#'     should be the embed_dim + 1 for pbs_simplex() and length(rel_lib_ind)
#'     for pbs_s_map().
#'
#' @return An tibble of ordered neighbour indices
#' 
#' @importFrom magrittr %>%
#' @export
#'
#' @examples 1:100 %>% pbs_make_tibble("x") %>% 
#'     pbs_make_lags("x", 5, 1) %>% pbs_calc_dist() %>%
#'     pbs_make_nbrs(50, c(6:49, 56:99), 5, 1, 6)
#' 
pbs_make_nbrs <- function(lag_dist, 
                          time_ind, 
                          rel_lib_ind, 
                          embed_dim, 
                          pred_dist,
                          max_nbrs = embed_dim + 1) {
  
  # Check arguments
  assertthat::assert_that(
    tibble::is_tibble(lag_dist),
    all.equal(names(lag_dist), c("focal_ind", "nbr_ind", "distance")),
    assertthat::is.count(time_ind),
    assertthat::is.count(embed_dim),
    assertthat::is.count(pred_dist),
    assertthat::is.count(max_nbrs),
    is.numeric(rel_lib_ind),
    is.vector(rel_lib_ind)
  )
  
  # Return tibble of ordered neighbours
  lag_dist %>%
    dplyr::filter(focal_ind %in% time_ind,
                  nbr_ind %in% rel_lib_ind) %>%
    dplyr::arrange(distance) %>%
    dplyr::mutate(dist_rank = dplyr::row_number()) %>%
    dplyr::filter(dist_rank <= max_nbrs) %>%
    dplyr::mutate(weight = exp(-distance / distance[1]),
                  proj_focal_ind = focal_ind + pred_dist,
                  proj_nbr_ind = nbr_ind + pred_dist)
}