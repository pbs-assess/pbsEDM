context("test-earlier-functions.R")

test_that("original function runs",{
  expect_is(EDM_pred_E_2(dplyr::select(Nx_lags_orig,
                                       -c("t", "my.pred"))),
            "list")
})
