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
#' @param rel_ind The allowable neighbour indices (numeric vector)
#' @param embed_dim The embedding dimension (numeric scalar)
#' @param pred_dist The prediction distance (numeric scalar)
#' @param max_nbrs The maximum number of neighbours (numeric scalar). It
#'     should be the embed_dim + 1 for pbs_simplex() and length(rel_ind)
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
                          rel_ind, 
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
    is.numeric(rel_ind),
    is.vector(rel_ind)
  )
  
  # Return tibble of ordered neighbours
  lag_dist %>%
    dplyr::filter(focal_ind %in% time_ind,
                  nbr_ind %in% rel_ind) %>%
    dplyr::arrange(distance) %>%
    dplyr::mutate(dist_rank = dplyr::row_number()) %>%
    dplyr::filter(dist_rank <= max_nbrs) %>%
    dplyr::mutate(focal_proj = focal_ind + pred_dist,
                  nbr_proj = nbr_ind + pred_dist)
}


#' Exclude Indicies Associated with a Specified Time Index
#'
#' @param time_ind The time index of the focal vector (numeric scalar)
#' @param pred_dist The prediction distance (numeric scalar)
#' @param lag_size The number of time steps separating successive lags
#' @param embed_dim The embedding dimension (numeric scalar)
#' @param one_sided Should only the focal value be excluded? (logical)
#'
#' @return A numeric vector
#'
#' @examples util_exclude_indices(10, 1, 1, 5, T)
#' 
util_exclude_indices <- function(time_ind, 
                                 pred_dist, 
                                 lag_size, 
                                 embed_dim,
                                 one_sided = T) {
  # Check arguments
  assertthat::assert_that(
    assertthat::is.count(time_ind),
    assertthat::is.count(pred_dist),
    assertthat::is.count(lag_size),
    assertthat::is.count(embed_dim),
    is.logical(one_sided)
  )
  if (one_sided) {
    seq(from = time_ind + pred_dist, by = lag_size, length.out = embed_dim)
  } else {
    sort(
      union(
        seq(from = time_ind + pred_dist, by = lag_size, length.out = embed_dim),  
        seq(from = time_ind + pred_dist, by = -lag_size, length.out = embed_dim)
      )
    )
  }
}

#' Make a Tibble of Lagged Columns of a Data Frame
#' 
#'
#' @param data A data frame with a named numeric columns to be lagged
#' @param lags A named list giving the lags to use for each corresponding
#'     column. List names must match column names. (list of numeric vectors)
#'
#' @return A tibble of lagged columns
#'
#' @importFrom magrittr %>%
#'
#' @examples 
#' dat <- data.frame(x = 1:10, y = 11:20, z = 21:30)
#' make_lag_tibble(dat, list(x = 0:2, y = 0:1, z = 0))
#' 
make_lag_tibble <- function(data, lags) {
  
  # Check arguments
  stopifnot(
    is.data.frame(data),
    is.list(lags),
    is.numeric(unlist(lags)),
    all(is.element(names(lags), names(data)))
  )
  
  # Initialize vectors
  lags_vector <- unlist(lags)
  col_vector <- rep(names(lags), lengths(lags))
  names_vector <- paste0(col_vector, "_lag", lags_vector)
  
  # Initialize data list
  data_list <- mapply(FUN = dplyr::pull,
                      var = col_vector,
                      MoreArgs = list(.data = data),
                      SIMPLIFY = FALSE)
  
  # Return lag tibble
  mapply(FUN = dplyr::lag,
         x = data_list,
         n = lags_vector,
         SIMPLIFY = FALSE) %>%
    dplyr::bind_cols() %>%
    magrittr::set_colnames(names_vector)
}

#' Make a Tibble of Euclidean Distances Between Rows of a Matrix
#'
#' @param mat A matrix of row vectors (matrix)
#'
#' @return A tibble of Euclidean distances
#' 
#' @importFrom magrittr %>%
#'
#' @examples make_dist_tibble(matrix(c(1:6, NA, 8:24), nrow = 6))
#' 
make_dist_tibble <- function(mat) {
  
  # Check argument
  stopifnot(
    is.numeric(mat),
    is.matrix(mat)
  )
  
  # Replace rows with NAs by NA rows
  mat_na <- mat %>% as.data.frame() %>% tibble::as_tibble() %>%
    dplyr::mutate(ind = seq_len(nrow(.))) %>%
    tidyr::drop_na() %>%
    dplyr::right_join(tibble::tibble(ind = seq_len(nrow(mat))), by = "ind") %>%
    dplyr::select(-ind) %>%
    as.matrix()

  # Calculate distances among row vectors
  stats::dist(mat_na, diag = FALSE, upper = TRUE) %>%
    broom::tidy() %>%
    dplyr::rename(focal = item2, nbr = item1) %>%
    dplyr::select(focal, nbr, distance) %>%
    tidyr::drop_na()
}

#' Make Global Allowable Indices from a Matrix of Lagged Column Row Vectors
#'
#' @param mat A matrix of lagged columns defining row vectors (matrix)
#' @param from A vector of indices to consider forecasting from (numeric)
#' @param into A vector of indices to consider forecasting into (numeric)
#' @param dist The forecast distance (numeric scalar)
#'
#' @return A tibble with two columns of indices
#' 
#' @importFrom magrittr %>%
#'
#' @examples
#' dat <- data.frame(x = c(1:7, NA, 9:15), y = 11:25)
#' mat <- as.matrix(combine_lag_tibbles(dat, c("x", "y"), c(3, 2), 1))
#' make_global_indices(mat)
#' 
make_global_indices <- function(mat, 
                                from = seq_len(nrow(mat)), 
                                dist = 1L) {
  
  # Check arguments
  stopifnot(
    is.matrix(mat),
    is.numeric(from),
    is.vector(from),
    is.numeric(dist)
  )
  
  # Make 'from' indices
  global_from <- tibble::as_tibble(mat) %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    dplyr::mutate(na_col = purrr::pmap_dbl(., sum, na.rm = FALSE)) %>%
    dplyr::mutate(buffer = dplyr::lead(na_col, dist)) %>%
    tidyr::drop_na() %>%
    dplyr::pull(index)
  
  # Make 'into' indices
  global_into <- tibble::as_tibble(mat) %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    dplyr::mutate(na_col = purrr::pmap_dbl(., sum, na.rm = FALSE)) %>%
    dplyr::mutate(buffer = dplyr::lag(na_col, dist)) %>%
    tidyr::drop_na() %>%
    dplyr::pull(index)
  
  # Return indices
  dplyr::bind_cols(from = global_from, into = global_into)
}

#' Make Local Allowable Indices for Nearest Neighbour Vectors
#'
#' @param index Index for the focal vector (numeric scalar)
#' @param from A vector of indices to consider forecasting from (numeric)
#' @param lags Named list of numeric vectors giving lags to use for each
#'     column. List names must match column names. (list of numeric vectors)
#' @param dist The forecast distance (numeric scalar)
#' @param symm Symmetric exclusion radius? (logical scalar)
#'
#' @return A vector of allowable indices for nearest neighbours
#' 
#' @importFrom magrittr %>%
#'
#' @examples make_local_indices(8, 1:15, list(x = 0:2))
#' 
make_local_indices <- function(index,
                               from,
                               lags,
                               dist = 1,
                               symm = FALSE) {
  
  # Check arguments
  stopifnot(
    is.numeric(index),
    is.numeric(from),
    is.vector(from),
    is.list(lags),
    is.numeric(unlist(lags)),
    is.numeric(dist),
    is.logical(symm)
  )
  
  # Specify unique lags
  lags_vector <- sort(unique(unlist(lags)))
  
  # Make local indices
  if (symm) {
    from %>% setdiff(index) %>%
      setdiff(index + dist + lags_vector) %>%
      setdiff(index + dist - lags_vector)
  } else {
    from %>% setdiff(index) %>% 
      setdiff(index + dist + lags_vector)
  }
}

make_neighbours <- function() {
  
}

make_forecast <- function(index,
                          lag_matrix,
                          dist_tibble,
                          global_indices) {
  
}

