context("test-methods.R")

test_that("pbsEDM() returns a list", {
	nt <- data.frame(x = simple_ts)
	lags <- list(x = 0:1)
	out_list <- pbsEDM(nt, lags)
	expect_true(is.list(out_list))
})
