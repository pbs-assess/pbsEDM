context("test-create-named-lags.R")

# Test create_named_lags ---------------------------------------------------

test_that("create_named_lags() generates correct lag names", {
  lags_list <- list(
    list(age11 = c(1)),
    list(age11 = c(1, 2)),
    list(age11 = c(1))  # duplicate
  )
  result <- create_named_lags(lags_list)

  expect_equal(
    result,
    c("age11_lag_1", "age11_lag_2")
  )
})

test_that("create_named_lags() handles multiple variables", {
  lags_list <- list(
    list(age11 = c(1)),
    list(temp = c(0, 1)),
    list(age11 = c(1, 2))
  )
  result <- create_named_lags(lags_list)

  expect_equal(
    result,
    c("age11_lag_1", "age11_lag_2", "temp_lag_0", "temp_lag_1")
  )
})

test_that("create_named_lags() removes duplicates correctly", {
  lags_list <- list(
    list(x = c(1)),
    list(x = c(1)),
    list(x = c(1)),
    list(y = c(2))
  )
  result <- create_named_lags(lags_list)

  expect_equal(
    result,
    c("x_lag_1", "y_lag_2")
  )
  expect_equal(length(result), 2)
})

test_that("create_named_lags() respects specified order parameter", {
  lags_list <- list(
    list(zebra = c(1)),
    list(apple = c(0, 1)),
    list(zebra = c(1, 2))
  )
  result <- create_named_lags(lags_list, order = c("zebra", "apple"))

  expect_equal(
    result,
    c("zebra_lag_1", "zebra_lag_2", "apple_lag_0", "apple_lag_1")
  )
})

test_that("create_named_lags() handles partial order specification", {
  lags_list <- list(
    list(zebra = c(1)),
    list(apple = c(1)),
    list(monkey = c(1))
  )
  result <- create_named_lags(lags_list, order = c("zebra"))

  expect_equal(
    result,
    c("zebra_lag_1", "apple_lag_1", "monkey_lag_1")
  )
})

test_that("create_named_lags() works properly on our example results", {
  result <- create_named_lags(age11_res$subset_lags)

  expect_equal(result,
               c("age11_lag_1",
                 "age11_lag_2",
                 "age11_lag_3",
                 "pdo_winter_mean_lag_1",
                 "pink_lag_1",
                 "pink_lag_2",
                 "sst_spring_lag_1"))
})
