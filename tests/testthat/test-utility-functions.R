context("test-utility-functions.R")

test_that("pbs_make_lags() returns a tibble with the correct number of rows", {
	embed_dim <- 10
	tbl <- tibble::tibble(x = 1:100)
	lag_tbl <- pbs_make_lags(tbl, "x", embed_dim, 1)
	expect_true(tibble::is_tibble(lag_tbl))
  expect_equal(nrow(lag_tbl), nrow(tbl))
})

	