context("test-plots.R")

test_that("plotting functions run on simple examples",{
  expect_null(plot_time_series(Nx_lags_orig$Nt,
                               X.or.N = "N"))

  s3d <- scatterplot3d::scatterplot3d(rnorm(5), rnorm(5), rnorm(5))
  expect_equal(length(gets3dusr(s3d)), 6)

  aa <- pbsEDM(Nx_lags_orig,
               lags = list(Nt = 0:1),
               first_difference = TRUE)
  expect_null(plot.pbsEDM(aa))


  aa_Evec <- pbsEDM_Evec(Nx_lags_orig$Nt)
  expect_null(plot_pbsEDM_Evec(aa_Evec))


})
