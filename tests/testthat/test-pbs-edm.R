context("test-pbs-edm.R")

test_that("pbs_edm() returns a tibble with the correct dimensions", {
  data_frame <- data.frame(x = simple_ts)
  lags <- list(x = 0:1)
  edm_tbl <- pbs_edm(data_frame, lags)
  expect_true(tibble::is_tibble(edm_tbl))
  expect_equal(nrow(edm_tbl), 1)
  expect_equal(ncol(edm_tbl), 9)
})
