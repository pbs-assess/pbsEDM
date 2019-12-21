context("test-pbs-previous.R")

test_that("pbs_simplex_1() returns a list of 2 tibbles and a list", {
  value <- pbs_simplex_1(tibble::tibble(x = 1:100), "x")
  expect_true(is.list(value))
  expect_equal(length(value), 3)
  expect_equal(names(value), c("stats_tbl", "pred_tbl", "nbr_list"))
  expect_true(tibble::is_tibble(value[[1]]))
  expect_true(tibble::is_tibble(value[[2]]))
  expect_true(is.list(value[[3]]))
})


test_that("pbs_s_map_1() returns a list of 2 tibbles and a list", {
  value <- pbs_s_map_1(tibble::tibble(x = 1:100), "x")
  expect_true(is.list(value))
  expect_equal(length(value), 3)
  expect_equal(names(value), c("stats_tbl", "pred_tbl", "nbr_list"))
  expect_true(tibble::is_tibble(value[[1]]))
  expect_true(tibble::is_tibble(value[[2]]))
  expect_true(is.list(value[[3]]))
})

test_that("pbs_edm_1() returns a tibble with the correct dimensions", {
  data_frame <- data.frame(x = simple_ts, y = simple_ts, z = simple_ts)
  lags <- list(x = 0:3, y = 0:1, z = 0)
  edm_00 <- pbs_edm_1(data_frame, lags)
  edm_01 <- pbs_edm_1(data_frame, lags, include_neighbours = FALSE)
  edm_02 <- pbs_edm_1(data_frame, 
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


test_that("pbs_smap_1() returns a tibble of correct dimensions", {
  data_frame <- data.frame(x = simple_ts)
  lags <- list(x = 0:3)
  local_weight <- 0
  smap_00 <- pbs_smap_1(data_frame, lags, local_weight)
  expect_true(tibble::is_tibble(smap_00))
  expect_true(nrow(smap_00) == 1)
  expect_true(ncol(smap_00) == 5)
})