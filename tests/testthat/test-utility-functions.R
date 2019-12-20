context("test-utility-functions.R")

test_that("pbs_make_lags() returns a tibble with the correct dimensions", {
	embed_dim <- 10
	tbl <- tibble::tibble(x = 1:100)
	lag_tbl <- pbs_make_lags(tbl, "x", embed_dim, 1)
	expect_true(tibble::is_tibble(lag_tbl))
  expect_equal(nrow(lag_tbl), nrow(tbl))
  expect_equal(ncol(lag_tbl), embed_dim)
})

test_that("pbs_calc_dist() returns a tibble with three columns", {
	lag_tbl <- pbs_make_lags(tibble::tibble(x = 1:100), "x", 5, 1)
	lag_dist <- pbs_calc_dist(lag_tbl)
	expect_true(tibble::is_tibble(lag_dist))
	expect_equal(ncol(lag_dist), 3)
	expect_equal(names(lag_dist), c("focal_ind", "nbr_ind", "distance"))
})

test_that("pbs_make_tibble returns a tibble with named column", {
	col_name <- "unlikely_name"
	num_vec <- 1:10
	old_tbl <- tibble::tibble(num_vec) %>% 
		magrittr::set_colnames(col_name) %>%
		pbs_make_tibble(col_name)
	new_tbl <- pbs_make_tibble(num_vec, col_name)
	expect_true(tibble::is_tibble(old_tbl))
	expect_true(col_name %in% names(old_tbl))
	expect_true(tibble::is_tibble(new_tbl))
	expect_true(col_name %in% names(new_tbl))
})

test_that("util_exclude_indices() returns a numeric vector", {
  call_true <- util_exclude_indices(10, 1, 1, 5, T)
  call_false <- util_exclude_indices(10, 1, 1, 5, F)
  expect_true(is.numeric(call_true))
  expect_true(is.vector(call_true))
  expect_true(is.numeric(call_false))
  expect_true(is.vector(call_false))
})

test_that("make_lag_tibble() returns a tibble of correct dimensions", {
  dat <- data.frame(x = 1:10, y = 11:20, z = 21:30)
  ll <- list(x = 0:2, y = 0:1, z = 0)
  num_rows <- nrow(dat)
  num_cols <- length(unlist(ll))
  lag_tibble <- make_lag_tibble(dat, ll)
  expect_true(tibble::is_tibble(lag_tibble))
  expect_true(nrow(lag_tibble) == num_rows)
  expect_true(ncol(lag_tibble) == num_cols)
})

test_that("make_dist_tibble() returns a tibble of correct dimensions", {
  tbl <- tibble::tibble(x = 1:6, y = c(NA, 8:12), z = 13:18)
  num_rows <- nrow(tidyr::drop_na(tbl)) * (nrow(tidyr::drop_na(tbl)) - 1)
  num_cols <- ncol(tbl)
  dist_tibble <- make_dist_tibble(tbl)
  expect_true(tibble::is_tibble(dist_tibble))
  expect_true(nrow(dist_tibble) == num_rows)
  expect_true(ncol(dist_tibble) == num_cols)
})

test_that("make_global_indices() returns a tibble of correct dimensions", {
  dat <- tibble::tibble(x = c(1:7, NA, 9:20), y = 21:40)
  lags <- list(x = 0:2, y = 0:1)
  lag_tibble <- make_lag_tibble(dat, lags)
  global_indices <- make_global_indices(lag_tibble)
  expect_true(tibble::is_tibble(global_indices))
  expect_true(nrow(global_indices) == 13)
  expect_true(ncol(global_indices) == 2)
})

test_that("make_local_indices() returns a vector of correct length", {
  from_index <- 8
  from_global <- 1:15
  lags <- list(x = 0:2)
  local_indices_f <- make_local_indices(from_index, from_global, lags)
  local_indices_t <- make_local_indices(from_index, from_global, lags, 1, TRUE)
  expect_true(is.vector(local_indices_f))
  expect_true(length(local_indices_f) == 11)
  expect_true(is.vector(local_indices_t))
  expect_true(length(local_indices_t) == 10)
})

test_that("make_neighbours() returns a tibble of correct dimensions", {
  tbl <- tibble::tibble(x = 1:20)
  lags <- list(x = 0:3)
  lag_tibble <- make_lag_tibble(tbl, lags)
  dist_tibble <- make_dist_tibble(lag_tibble)
  from_global <- make_global_indices(lag_tibble) %>% dplyr::pull(from)
  from_local <- make_local_indices(11, from_global, lags)
  nbs <- make_neighbours(15, from_local, dist_tibble, length(unlist(lags)) + 1)
  expect_true(tibble::is_tibble(nbs))
  expect_true(nrow(nbs) == length(unlist(lags)) + 1)
  expect_true(ncol(nbs) == 6)
})

test_that("make_simplex_forecast() returns a tibble of correct dimensions", {
  tbl <- tibble::tibble(x = simple_ts)
  from_index <- 15
  lags <- list(x = 0:3)
  lag_tibble <- make_lag_tibble(tbl, lags)
  from_global <- dplyr::pull(make_global_indices(lag_tibble), from)
  distance_tibble <- make_dist_tibble(lag_tibble)
  forecast_distance <- 1
  forecast <- make_simplex_forecast(from_index,
                                    from_global,
                                    lags,
                                    lag_tibble,
                                    distance_tibble,
                                    forecast_distance,
                                    symmetric_exclusion = FALSE)
  expect_true(tibble::is_tibble(forecast))
  expect_true(nrow(forecast) == 1)
  expect_true(ncol(forecast) == 4)
}) 

test_that("make_smap_forecast() returns a tibble of correct dimensions", {
  tbl <- tibble::tibble(x = simple_ts)
  from_index <- 15
  lags <- list(x = 0:3)
  local_weight <- 0
  lag_tibble <- make_lag_tibble(tbl, lags)
  from_global <- dplyr::pull(make_global_indices(lag_tibble), from)
  distance_tibble <- make_dist_tibble(lag_tibble)
  forecast_distance <- 1
  forecast <- make_smap_forecast(from_index,
                                 from_global,
                                 lags,
                                 local_weight,
                                 lag_tibble,
                                 distance_tibble,
                                 forecast_distance,
                                 symmetric_exclusion = FALSE)
  expect_true(tibble::is_tibble(forecast))
  expect_true(nrow(forecast) == 1)
  expect_true(ncol(forecast) == 4)
}) 
