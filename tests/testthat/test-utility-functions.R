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


