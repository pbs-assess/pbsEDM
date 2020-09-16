context("test-plots.R")

test_that("plot_time_series() runs on example",{
  expect_null(plot_time_series(Nx_lags_orig$Nt,
                               X.or.N = "N"))
})
