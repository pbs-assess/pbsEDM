context("test-pbs-edm.R")

test_that("pbs_edm() returns a tibble with the correct dimensions", {
  data_frame <- data.frame(x = simple_ts, y = simple_ts, z = simple_ts)
  lags <- list(x = 0:3, y = 0:1, z = 0)
  edm_00 <- pbs_edm(data_frame, lags)
  edm_01 <- pbs_edm(data_frame, lags, include_neighbours = FALSE)
  edm_02 <- pbs_edm(data_frame, 
                    lags, 
                    include_forecasts = FALSE,
                    include_neighbours = FALSE)
  expect_true(tibble::is_tibble(edm_00))
  expect_true(tibble::is_tibble(edm_01))
  expect_true(tibble::is_tibble(edm_02))
  expect_true(nrow(edm_00) == 1)
  expect_true(nrow(edm_01) == 1)
  expect_true(nrow(edm_02) == 1)
  expect_true(ncol(edm_00) == 5)
  expect_true(ncol(edm_01) == 4)
  expect_true(ncol(edm_02) == 3)
})
