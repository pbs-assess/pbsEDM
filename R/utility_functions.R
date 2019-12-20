#' Make a Tibble of Lagged Columns of a Data Frame
#' 
#'
#' @param data_frame A data frame with a named numeric columns to be lagged
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
make_lag_tibble <- function(data_frame, lags) {
  
  # Check arguments
  stopifnot(
    is.data.frame(data_frame),
    is.list(lags),
    is.numeric(unlist(lags)),
    all(is.element(names(lags), names(data_frame)))
  )
  
  # Initialize vectors
  lags_vector <- unlist(lags)
  col_vector <- rep(names(lags), lengths(lags))
  names_vector <- paste0(col_vector, "_lag", lags_vector)
  
  # Initialize data list
  data_list <- mapply(FUN = dplyr::pull,
                      var = col_vector,
                      MoreArgs = list(.data = data_frame),
                      SIMPLIFY = FALSE)
  
  # Return lag tibble
  mapply(FUN = dplyr::lag,
         x = data_list,
         n = lags_vector,
         SIMPLIFY = FALSE) %>%
    dplyr::bind_cols() %>%
    magrittr::set_colnames(names_vector)
}

#' Make a Tibble of Euclidean Distances Between Rows of a Tibble
#'
#' @param tbl A tibble of row vectors (tibble)
#'
#' @return A tibble of Euclidean distances
#' 
#' @importFrom magrittr %>%
#'
#' @examples 
#'   tbl <- tibble::tibble(x = 1:6, y = c(NA, 8:12), z = 13:18)
#'   dist_tibble <- make_dist_tibble(tbl)
#' 
make_dist_tibble <- function(tbl) {
  
  # Check argument
  stopifnot(
    tibble::is_tibble(tbl)
  )
  
  # Initialize matrix; replace rows with NAs by NA rows
  mat_na <- tbl %>% 
    dplyr::mutate(ind = seq_len(nrow(.))) %>%
    tidyr::drop_na() %>%
    dplyr::right_join(tibble::tibble(ind = seq_len(nrow(tbl))), by = "ind") %>%
    dplyr::select(-ind) %>%
    as.matrix()

  # Calculate distances among row vectors
  stats::dist(mat_na, diag = FALSE, upper = TRUE) %>%
    broom::tidy() %>%
    dplyr::rename(from_focal = item2, from_nbr = item1) %>%
    dplyr::select(from_focal, from_nbr, distance) %>%
    tidyr::drop_na()
}

#' Make Global Allowable Indices from a Data Frame of Lagged Column Row Vectors
#'
#' @param lag_tibble A data frame of lagged columns defining row vectors
#' @param from A vector of indices to consider forecasting from (numeric)
#' @param into A vector of indices to consider forecasting into (numeric)
#' @param dist The forecast distance (numeric scalar)
#'
#' @return A tibble with two columns of indices
#' 
#' @importFrom magrittr %>%
#'
#' @examples
#' dat <- data.frame(x = c(1:7, NA, 9:20), y = 21:40)
#' lag_tibble <- make_lag_tibble(dat, list(x = 0:2, y = 0:1))
#' make_global_indices(lag_tibble)
#' 
make_global_indices <- function(lag_tibble, 
                                from = seq_len(nrow(lag_tibble)), 
                                dist = 1L) {
  
  # Check arguments
  stopifnot(
    tibble::is_tibble(lag_tibble),
    is.numeric(from),
    is.vector(from),
    is.numeric(dist)
  )
  
  # Make 'from' indices
  from_global <- lag_tibble %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    dplyr::mutate(na_col = purrr::pmap_dbl(., sum, na.rm = FALSE)) %>%
    dplyr::mutate(buffer = dplyr::lead(na_col, dist)) %>%
    tidyr::drop_na() %>%
    dplyr::pull(index)
  
  # Make 'into' indices
  into_global <- lag_tibble %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    dplyr::mutate(na_col = purrr::pmap_dbl(., sum, na.rm = FALSE)) %>%
    dplyr::mutate(buffer = dplyr::lag(na_col, dist)) %>%
    tidyr::drop_na() %>%
    dplyr::pull(index)
  
  # Return indices
  dplyr::bind_cols(from = from_global, into = into_global)
}

#' Make Local Allowable Indices for Nearest Neighbour Vectors
#'
#' @param from_index Index for the focal vector (numeric scalar)
#' @param from_global A vector of indices to consider forecasting from (numeric)
#' @param lags Named list of numeric vectors giving lags for each column. 
#'     List names must match column names. (list of numeric vectors)
#' @param forecast_distance The forecast distance (numeric scalar)
#' @param symmetric_exclusion Symmetric exclusion radius? (logical scalar)
#'
#' @return A vector of allowable indices for nearest neighbours
#' 
#' @importFrom magrittr %>%
#'
#' @examples make_local_indices(8, 1:15, list(x = 0:2))
#' 
make_local_indices <- function(from_index,
                               from_global,
                               lags,
                               forecast_distance = 1L,
                               symmetric_exclusion = FALSE) {
  
  # Check arguments
  stopifnot(
    is.numeric(from_index),
    is.numeric(from_global),
    is.vector(from_global),
    is.list(lags),
    is.numeric(unlist(lags)),
    is.numeric(forecast_distance),
    is.logical(symmetric_exclusion)
  )
  
  # Specify unique lags
  lags_vector <- sort(unique(unlist(lags)))
  
  # Make local indices
  if (symmetric_exclusion) {
    from_global %>% setdiff(from_index) %>%
      setdiff(from_index + forecast_distance + lags_vector) %>%
      setdiff(from_index + forecast_distance - lags_vector)
  } else {
    from_global %>% setdiff(from_index) %>% 
      setdiff(from_index + forecast_distance + lags_vector)
  }
}

#' Make a Tibble of Nearest Neighbour Indices
#'
#' @param from_index The time index of the focal vector
#' @param from_local The allowable neighbour indices (numeric vector)
#' @param dist_tibble A tibble with columns 'focal', 'nbr' and 'distance'
#' @param max_neighbours The maximum number of neighbours
#' @param forecast_distance The forecast distance (numeric scalar)
#'
#' @return A tibble
#'
#' @examples
#'     dist_tibble <- tibble::tibble(x = 1:100) %>% 
#'     make_lag_tibble(list(x = 0:4)) %>% 
#'     make_dist_tibble()
#'     make_neighbours(50, c(6:49, 56:99), dist_tibble, 6)
#' 

make_neighbours <- function(from_index,
                            from_local,
                            distance_tibble,
                            max_neighbours,
                            forecast_distance = 1L) {
  # Check arguments
  stopifnot(
    is.numeric(from_index),
    is.numeric(from_local),
    is.vector(from_local),
    is.numeric(forecast_distance),
    tibble::is_tibble(distance_tibble),
    is.numeric(max_neighbours)
  )
  
  # Return a tibble of neighbours ordered by distance
  distance_tibble %>%
    dplyr::filter(from_focal %in% from_index, from_nbr %in% from_local) %>%
    dplyr::arrange(distance) %>%
    dplyr::mutate(dist_rank = dplyr::row_number()) %>%
    dplyr::filter(dist_rank <= max_neighbours) %>%
    dplyr::mutate(into_focal = from_focal + forecast_distance, 
                  into_nbr = from_nbr + forecast_distance)
}


#' Make a Simplex Projection Forecast from the Current Time Index
#'
#' @param from_index The index from which to forecast (numeric scalar)
#' @param from_global A superset of neighbour indices (numeric vector)
#' @param lags Named list of numeric vectors giving lags for each column. 
#'     List names must match column names. (list of numeric vectors)
#' @param lag_tibble A tibble of lagged numeric columns (tibble)
#' @param distance_tibble A long-format tibble of Euclidean distances between
#'     row vectors in lag_tibble (tibble)
#' @param forecast_distance The forecast distance (numeric scalar)
#' @param symmetric_exclusion Remove local indices symmectric around the
#'     forecast value? (Logical scalar)
#'
#' @return A tibble with numeric columns: index, obsertation, and forecast; and
#'     a list column with the neighbours tibble
#'
#' @examples
#' tbl <- tibble::tibble(x = simple_ts)
#' from_index <- 15
#' lags <- list(x = 0:3)
#' lag_tibble <- make_lag_tibble(tbl, lags)
#' from_global <- dplyr::pull(make_global_indices(lag_tibble), from)
#' distance_tibble <- make_dist_tibble(lag_tibble)
#' forecast_distance <- 1
#' forecast <- make_simplex_forecast(from_index,
#'                                   from_global,
#'                                   lags,
#'                                   distance_tibble,
#'                                   forecast_distance,
#'                                   symmetric_exclusion = FALSE)
#' 
make_simplex_forecast <- function(from_index,
                                  from_global,
                                  lags,
                                  lag_tibble,
                                  distance_tibble,
                                  forecast_distance,
                                  symmetric_exclusion = FALSE) {
  
  # Check arguments
  stopifnot(
    is.numeric(from_index),
    is.numeric(from_global),
    is.vector(from_global),
    is.list(lags),
    is.numeric(unlist(lags)),
    tibble::is_tibble(lag_tibble),
    tibble::is_tibble(distance_tibble),
    is.logical(symmetric_exclusion)
  )
  
  # Specify local indices
  from_local <- make_local_indices(from_index,
                                   from_global,
                                   lags,
                                   forecast_distance,
                                   symmetric_exclusion)
  
  # Identify nearest neighbours
  nbr <- make_neighbours(from_index,
                         from_local,
                         distance_tibble,
                         max_neighbours = length(unlist(lags)) + 1,
                         forecast_distance)
  
  # Calculate weights
  nbr_wts <- nbr %>% dplyr::mutate(weight = exp(-distance / distance[1]))
  
  # Calculate into index
  into_index <- from_index + forecast_distance
  
  # Forecast value
  if (nrow(nbr_wts) > length(unlist(lags))) {
    
    # Pull vectors
    into_neighbour_indices <- nbr_wts %>% dplyr::pull(into_nbr)
    into_neighbour_values <- lag_tibble %>% 
      dplyr::filter(dplyr::row_number() %in% into_neighbour_indices) %>%
      dplyr::mutate(row_order = order(into_neighbour_indices)) %>%
      dplyr::arrange(row_order) %>%
      dplyr::pull(1)
    into_neighbour_weights <- nbr_wts %>% dplyr::pull(weight)
    
    # Make forecast
    into_forecast <- sum(into_neighbour_values * into_neighbour_weights) /
      sum(into_neighbour_weights)

  } else {
    into_forecast <- NA_real_
  }
  
  # Make forecast
  tibble::tibble(
    index = into_index,
    observation = lag_tibble[[from_index + forecast_distance, 1]],
    forecast = into_forecast,
    neighbours = list(nbr_wts)
  )
}

#' Make a S-Map Forecast from the Current Time Index
#'
#' @param from_index The index from which to forecast (numeric scalar)
#' @param from_global A superset of neighbour indices (numeric vector)
#' @param lags Named list of numeric vectors giving lags for each column. 
#'     List names must match column names. (list of numeric vectors)
#' @param local_weight State-dependence parameter (numeric scalar)
#' @param lag_tibble A tibble of lagged numeric columns (tibble)
#' @param distance_tibble A long-format tibble of Euclidean distances between
#'     row vectors in lag_tibble (tibble)
#' @param forecast_distance The forecast distance (numeric scalar)
#' @param symmetric_exclusion Remove local indices symmectric around the
#'     forecast value? (Logical scalar)
#'
#' @return A tibble with numeric columns: index, obsertation, and forecast; and
#'     a list column with the neighbours tibble
#'
#' @examples
#' tbl <- tibble::tibble(x = simple_ts)
#' from_index <- 15
#' lags <- list(x = 0:3)
#' local_weight <- 0
#' lag_tibble <- make_lag_tibble(tbl, lags)
#' from_global <- dplyr::pull(make_global_indices(lag_tibble), from)
#' distance_tibble <- make_dist_tibble(lag_tibble)
#' forecast_distance <- 1
#' forecast <- make_smap_forecast(from_index,
#'                                from_global,
#'                                lags,
#'                                local_weight,
#'                                lag_tibble,
#'                                distance_tibble,
#'                                forecast_distance,
#'                                symmetric_exclusion = FALSE)
#' 
make_smap_forecast <- function(from_index,
                               from_global,
                               lags,
                               local_weight,
                               lag_tibble,
                               distance_tibble,
                               forecast_distance,
                               symmetric_exclusion = FALSE) {
  
  # Check arguments
  stopifnot(
    is.numeric(from_index),
    is.numeric(from_global),
    is.vector(from_global),
    is.list(lags),
    is.numeric(unlist(lags)),
    is.numeric(local_weight),
    tibble::is_tibble(lag_tibble),
    tibble::is_tibble(distance_tibble),
    is.logical(symmetric_exclusion)
  )
  
  # Specify local indices
  from_local <- make_local_indices(from_index,
                                   from_global,
                                   lags,
                                   forecast_distance,
                                   symmetric_exclusion)
  
  # Identify nearest neighbours
  nbr <- make_neighbours(from_index,
                         from_local,
                         distance_tibble,
                         max_neighbours = nrow(lag_tibble),
                         forecast_distance)
  
  # Calculate weights
  nbr_wts <- nbr %>% dplyr::mutate(
    weight = exp(-local_weight * distance / mean(distance, na.rm = TRUE))
  )

  # Calculate into index
  into_index <- from_index + forecast_distance
  
  # Compute number of neighbours and embedding dimension
  num_nbrs <- nrow(nbr_wts)
  embedding_dimension <- length(unlist(lags))
  
  # Forecast value
  if (num_nbrs > embedding_dimension) {
    
    # Pull vectors
    into_neighbour_indices <- nbr_wts %>% dplyr::pull(into_nbr)
    into_neighbour_values <- lag_tibble %>% 
      dplyr::filter(dplyr::row_number() %in% into_neighbour_indices) %>%
      dplyr::mutate(row_order = order(into_neighbour_indices)) %>%
      dplyr::arrange(row_order) %>%
      dplyr::pull(1)
    into_neighbour_weights <- nbr_wts %>% dplyr::pull(weight)
    
    # Calculate the B vector and A matrix
    a_mat <- matrix(NA_real_, nrow = num_nbrs, ncol = embedding_dimension)
    b_vec <- rep(NA_real_, num_nbrs)
    for (i in seq_len(num_nbrs)) {
      b_vec[i] <- nbr_wts$weight[i] * into_neighbour_values[i]
      for (j in seq_len(embedding_dimension)) {
        a_mat[i, j] <- nbr_wts$weight[i] * lag_tibble[[nbr_wts$from_nbr[i], j]]
      }
    }
    
    # Solve for C using Singluar Value Decomposition
    a_svd <- svd(a_mat) # m x n
    d_inv <- diag(1 / a_svd$d) # diag n x n
    u_mat <- a_svd$u # m x n
    v_mat <- a_svd$v # n x n
    c_vec <- v_mat %*% d_inv %*% t(u_mat) %*% b_vec

    # Make forecast
    into_forecast <- sum(c_vec * lag_tibble[from_index, ])
      
  } else {
    into_forecast <- NA_real_
  }
  
  # Make forecast
  tibble::tibble(
    index = into_index,
    observation = lag_tibble[[from_index + forecast_distance, 1]],
    forecast = into_forecast,
    neighbours = list(nbr_wts)
  )
}



#------------------ Functions to Remove Below Here ----------------------------#



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




