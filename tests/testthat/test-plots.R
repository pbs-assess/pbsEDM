context("test-plots.R")

test_that("plotting functions run on simple examples",{
  expect_invisible(plot_time_series(Nx_lags_orig$Nt,
                               X.or.N = "N"))

  s3d <- scatterplot3d::scatterplot3d(rnorm(5), rnorm(5), rnorm(5))
  expect_equal(length(gets3dusr(s3d)), 6)

  aa <- pbsEDM(Nx_lags_orig,
               lags = list(Nt = 0:1),
               first_difference = TRUE)
  expect_invisible(plot.pbsEDM(aa,
                          late.col = "pink")) # to use "...", but they don't
                                        # need to get tested (don't show up in
                                        # red on the covr report). But "pink" works.

  # Change some defaults for functions that get called above with plot.pbsEDM()
  expect_invisible(plot_phase_2d(Nx_lags_orig$Nt,
                            X.or.N = "N",
                            cobwebbing = FALSE))
  expect_invisible(plot_phase_3d(aa,
                            last.time.to.plot = 1))

  aa_Evec <- pbsEDM_Evec(Nx_lags_orig$Nt)
  expect_invisible(plot_pbsEDM_Evec(aa_Evec))
  expect_invisible(plot_rho_Evec(aa_Evec))

  expect_is(plot_explain_edm_movie(aa,
                                   tstar = 15),
            "list")
})
