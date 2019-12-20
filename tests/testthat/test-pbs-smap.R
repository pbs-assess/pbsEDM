context("test-pbs-smap.R")

test_that("pbs_smap() returns a tibble of correct dimensions", {
  data_frame <- data.frame(x = simple_ts)
  lags <- list(x = 0:3)
  local_weight <- 0
  smap_00 <- pbs_smap(data_frame, lags, local_weight)
  expect_true(tibble::is_tibble(smap_00))
  expect_true(nrow(smap_00) == 1)
  expect_true(ncol(smap_00) == 5)
})
