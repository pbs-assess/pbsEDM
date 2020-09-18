context("test-methods.R")

test_that("pbsEDM() returns a list", {
	N <- data.frame(x = simple_ts)
	lags <- list(x = 0:1)
	out_list <- pbsEDM(N, lags)
	expect_true(is.list(out_list))

        out_list2 <- pbsEDM(N,
                            lags,
                            centre_and_scale = TRUE)
	expect_true(is.list(out_list2))
})

test_that("pbsSmap() returns a list", {
	N <- data.frame(x = simple_ts)
	lags <- list(x = 0:1)
	out_list <- pbsSmap(N, lags)
	expect_true(is.list(out_list))

      	out_list2 <- pbsSmap(N,
                             lags,
                             centre_and_scale = TRUE)
	expect_true(is.list(out_list2))

      	out_list3 <- pbsSmap(N,
                             lags,
                             first_difference = TRUE)
	expect_true(is.list(out_list3))
})
