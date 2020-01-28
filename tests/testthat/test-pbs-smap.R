context("test-pbs-smap.R")

test_that("pbs_smap() returns a tibble with the correct dimensions", {
	data_frame <- data.frame(x = simple_ts)
	lags <- list(x = 0:1)
	smap_list <- pbs_smap(data_frame, lags, show_calculations = TRUE)
	expect_true(is.list(smap_list))
	expect_equal(length(smap_list), 12)
})
