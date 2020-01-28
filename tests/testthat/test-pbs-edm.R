context("test-pbs-edm.R")

test_that("pbs_edm() returns a list with the correct length", {
  data_frame <- data.frame(x = simple_ts)
  lags <- list(x = 0:1)
  edm_list <- pbs_edm(data_frame, lags, show_calculations = TRUE)
  expect_true(is.list(edm_list))
  expect_equal(length(edm_list), 12)
})
