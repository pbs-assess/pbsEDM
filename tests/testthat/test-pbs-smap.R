context("test-pbs-smap.R")

test_that("pbs_smap() returns a tibble with the correct dimensions", {
	data_frame <- data.frame(x = simple_ts)
	lags <- list(x = 0:1)
	smap_tbl <- pbs_smap(data_frame, lags)
	expect_true(tibble::is_tibble(smap_tbl))
	expect_equal(nrow(smap_tbl), 1)
	expect_equal(ncol(smap_tbl), 8)
})