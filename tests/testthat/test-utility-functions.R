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

