context("test-methods.R")

test_that("pbsEDM() returns a list", {
	nt <- data.frame(x = simple_ts)
	lags <- list(x = 0:1)
	out_list <- pbsEDM(nt, lags)
	expect_true(is.list(out_list))

        out_list2 <- pbsEDM(nt,
                            lags,
                            centre_and_scale = TRUE)
	expect_true(is.list(out_list2))
})

test_that("pbsSMAP() returns a list", {
	nt <- data.frame(x = simple_ts)
	lags <- list(x = 0:1)
	out_list <- pbsSMAP(nt, lags)
	expect_true(is.list(out_list))

      	out_list2 <- pbsSMAP(nt,
                             lags,
                             centre_and_scale = TRUE)
	expect_true(is.list(out_list2))

      	out_list3 <- pbsSMAP(nt,
                             lags,
                             first_difference = TRUE)
	expect_true(is.list(out_list3))
})
