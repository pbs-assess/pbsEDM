context("test-utilities.R")

# Test pbsN --------------------------------------------------------------------

test_that("pbsN() returns a matrix", {
	# Vector
	N <- 1:10
	lags <- list(x = c(0, 1, 2))
	expect_true(class(pbsN(N, lags))[1] == "matrix")	
	# Matrix
	N <- as.matrix(1:10)
	colnames(N) <- "x"
	lags <- list(x = c(0, 1, 2))
	expect_true(class(pbsN(N, lags))[1] == "matrix")	
	# Data frame
	N <- data.frame(x = 1:10)
	lags <- list(x = c(0, 1, 2))
	expect_true(class(pbsN(N, lags))[1] == "matrix")	
})

# Test pbsZ --------------------------------------------------------------------

test_that("pbsZ() returns a matrix", {
	lags <- list(x = c(0, 1, 2))
	N <- pbsN(1:10, lags)
	expect_true(class(pbsZ(N, FALSE))[1] == "matrix")	
	expect_true(class(pbsZ(N, TRUE))[1] == "matrix")	
})

# Test pbsY --------------------------------------------------------------------

test_that("pbsY() returns a matrix", {
	lags <- list(x = c(0, 1, 2))
	N <- pbsN(sample(1:10), lags)
	Z <- pbsZ(N, TRUE)
	expect_true(class(pbsY(Z, FALSE))[1] == "matrix")	
	expect_true(class(pbsY(Z, TRUE))[1] == "matrix")	
})

# Test pbsX --------------------------------------------------------------------

test_that("pbsX() returns a matrix", {
	lags <- list(x = c(0, 1, 2))
	N <- pbsN(sample(1:10), lags)
	Z <- pbsZ(N, TRUE)
	Y <- pbsY(Z, TRUE)
	expect_true(class(pbsX(Y, lags))[1] == "matrix")	
})

# Test pbsSSR ------------------------------------------------------------------

# test_that("pbsSSR() returns a matrix of class 'pbsSSR'", {
# 	lags <- list(x = c(0, 1, 2))
# 	N <- pbsN(sample(1:10), lags)
# 	expect_true(class(pbsSSR(N, lags))[1] == "matrix")	
# 	expect_true(class(pbsSSR(N, lags))[2] == "pbsSSR")	
# })

# Test pbsLag ------------------------------------------------------------------

test_that("pbsLag() returns a vector for vector input", {
	x <- 0
	y <- pbsLag(x)
	expect_vector(y)
})

test_that("pbsLag() returns a matrix for matrix input", {
	x <- matrix(0, nrow = 1, ncol = 1)
	y <- pbsLag(x)
	expect_true(is.matrix(y))
})
# TODO: Update from here
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





