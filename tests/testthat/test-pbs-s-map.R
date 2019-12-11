context("test-pbs-s-map.R")

test_that("pbs_s_map() returns a list of 2 tibbles and a list", {
  value <- pbs_s_map(tibble::tibble(x = 1:100), "x")
  expect_true(is.list(value))
  expect_equal(length(value), 3)
  expect_equal(names(value), c("stats_tbl", "pred_tbl", "nbr_list"))
  expect_true(tibble::is_tibble(value[[1]]))
  expect_true(tibble::is_tibble(value[[2]]))
  expect_true(is.list(value[[3]]))
})
