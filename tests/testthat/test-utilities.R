context("test-utilities.R")

test_that("pbsLAG() returns a vector or matrix depending on the input", {
	x1 <- c(1:10)
	x2 <- matrix(rep(1:10, 2), ncol = 2)
	m1 <- pbsLAG(x1, 1)
	m2 <- pbsLAG(x2, c(1, 2))
	expect_true(is.numeric(m1))
	expect_true(is.matrix(m2))
})	

test_that("pbsLAG() returns the correct valued vector or matrix", {
	x1 <- c(1:10)
	x2 <- matrix(rep(1:10, 2), ncol = 2)
	t1 <- c(NA, 1:9)
	t2 <- matrix(c(NA, 1:9, NA, NA, 1:8), ncol = 2)
	m1 <- pbsLAG(x1, 1)
	m2 <- pbsLAG(x2, c(1, 2))
	expect_equal(m1, t1)
	expect_equal(m2, t2)
})	

test_that("pbsLAG() recycles input lags of length one only", {
	x3 <- matrix(rep(1:10, 3), ncol = 3)
	t31 <- matrix(rep(c(NA, NA, NA, 1:7), 3), ncol = 3)
	t32 <- matrix(c(NA, NA, NA, 1:7, NA, NA, 1:8, 1:10), ncol = 3)
	t33 <- matrix(c(NA, NA, NA, 1:7, NA, NA, 1:8, NA, 1:9), ncol = 3)
	m31 <- pbsLAG(x3, 3)
	m32 <- pbsLAG(x3, 3:2)
	m33 <- pbsLAG(x3, 3:1)
	expect_equal(m31, t31)
	expect_equal(m32, t32)
	expect_equal(m33, t33)
})

