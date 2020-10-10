context("test-utilities.R")
# pbsN()
test_that("pbsN() errors when N is a <0 x 0 matrix>", {
	
})
test_that("pbsN() errors when N is a data frame with 0 columns and 0 rows", {
	
})
test_that("pbsN() errors when the columns of N are not numeric", {
	
})
test_that("pbsN() errors when the elements of lags are not numeric", {
	
})
test_that("pbsN() errors when the elements of lags are not vectors", {
	
})
test_that("pbsN() errors when names(lags) are not in colnames(N)", {
	
})
test_that("pbsN() errors when there are duplicates in names(lags)", {
	
})
test_that("pbsN() errors when there are duplicates in colnames(N)", {
	
})
test_that("pbsN() errors when lags[[1]][1] is not equal to 0L", {
	
})
test_that("pbsN() errors when p is not a scalar integer", {
	
})
# pbsZ()


test_that("pbsLag() returns a vector or matrix depending on the input", {
	x1 <- c(1:10)
	x2 <- matrix(rep(1:10, 2), ncol = 2)
	m1 <- pbsLag(x1, 1)
	m2 <- pbsLag(x2, c(1, 2))
	expect_true(is.numeric(m1))
	expect_true(is.matrix(m2))
})	

test_that("pbsLag() returns the correct valued vector or matrix", {
	x1 <- c(1:10)
	x2 <- matrix(rep(1:10, 2), ncol = 2)
	t1 <- c(NA, 1:9)
	t2 <- matrix(c(NA, 1:9, NA, NA, 1:8), ncol = 2)
	m1 <- pbsLag(x1, 1)
	m2 <- pbsLag(x2, c(1, 2))
	expect_equal(m1, t1)
	expect_equal(m2, t2)
})	

test_that("pbsLag() recycles input lags of length one only", {
	x3 <- matrix(rep(1:10, 3), ncol = 3)
	t31 <- matrix(rep(c(NA, NA, NA, 1:7), 3), ncol = 3)
	t32 <- matrix(c(NA, NA, NA, 1:7, NA, NA, 1:8, 1:10), ncol = 3)
	t33 <- matrix(c(NA, NA, NA, 1:7, NA, NA, 1:8, NA, 1:9), ncol = 3)
	m31 <- pbsLag(x3, 3)
	m32 <- pbsLag(x3, 3:2)
	m33 <- pbsLag(x3, 3:1)
	expect_equal(m31, t31)
	expect_equal(m32, t32)
	expect_equal(m33, t33)
})

