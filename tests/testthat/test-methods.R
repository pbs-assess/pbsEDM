context("test-methods.R")

# pbsEDM
test_that("pbsEDM() returns a list of class 'pbsEDM'", {
	N <- data.frame(N1 = 1:10)
	lags <- list(N1 = 0:1)
	m1 <- pbsEDM(N, lags)
	expect_true(is.list(m1))
	expect_s3_class(m1, "pbsEDM")
})

test_that("pbsEDM() throws an error if lags[[1]][1] != 0L", {
	N <- data.frame(N1 = 1:10)
	lags <- list(N1 = 1:2)
	expect_error(pbsEDM(N, lags))
})

# pbsSmap
test_that("pbsSmap() returns a list of class 'pbsEDM' and 'pbsSmap'", {
	N <- data.frame(N1 = 1:10)
	lags <- list(N1 = 0:1)
	m1 <- pbsSmap(N, lags)
	expect_true(is.list(m1))
	expect_s3_class(m1, "pbsEDM")
	expect_s3_class(m1, "pbsSmap")
})

test_that("pbsSmap() throws an error if lags[[1]][1] != 0L", {
	N <- data.frame(N1 = 1:10)
	lags <- list(N1 = 1:2)
	expect_error(pbsSmap(N, lags))
})
