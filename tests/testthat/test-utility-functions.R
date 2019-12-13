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
  num_rows <- 16
  num_cols <- 3
  lag_tibble <- make_lag_tibble(data.frame(x = 1:num_rows), "x", num_cols, 2)
  expect_true(tibble::is_tibble(lag_tibble))
  expect_true(nrow(lag_tibble) == num_rows)
  expect_true(ncol(lag_tibble) == num_cols)
})

test_that("combine_lag_tibbles() returns a tibble of correct dimensions", {
  dat <- data.frame(x = 1:15, y = 11:25, z = 21:35)
  dims <- c(3, 2, 1)
  num_rows <- nrow(dat)
  num_cols <- sum(dims)
  lag_tibble <- combine_lag_tibbles(dat, c("x", "y", "z"), dims, 1)
  expect_true(tibble::is_tibble(lag_tibble))
  expect_true(nrow(lag_tibble) == num_rows)
  expect_true(ncol(lag_tibble) == num_cols)
})

test_that("make_dist_tibble() returns a tibble of correct dimensions", {
  vals <- c(1:6, NA, 8:24)
  rows <- 6
  mat <- matrix(vals, rows)
  num_rows <- (5 - 1) * (6 - 1)
  num_cols <- 3
  dist_tibble <- make_dist_tibble(mat)
  expect_true(tibble::is_tibble(dist_tibble))
  expect_true(nrow(dist_tibble) == num_rows)
  expect_true(ncol(dist_tibble) == num_cols)
})

test_that("make_global_indices() returns a tibble of correct dimensions", {
  x_vals <- c(1:7, NA, 9:15)
  y_vals <- 11:25
  dat <- data.frame(x = x_vals, y = y_vals)
  dims <- c(3, 2)
  mat <- as.matrix(combine_lag_tibbles(dat, c("x", "y"), dims, 1))
  global_indices <- make_global_indices(mat)
  expect_true(tibble::is_tibble(global_indices))
  expect_true(nrow(global_indices) == 8)
  expect_true(ncol(global_indices) == 2)
})

test_that("make_local_indices() returns a vector of correct length", {
  index <- 8
  dim <- 3
  from <- 1:15
  local_indices_f <- make_local_indices(index, dim, from)
  local_indices_t <- make_local_indices(index, dim, from, symm = TRUE)
  expect_true(is.vector(local_indices_f))
  expect_true(length(local_indices_f) == 11)
  expect_true(is.vector(local_indices_t))
  expect_true(length(local_indices_t) == 10)
})




