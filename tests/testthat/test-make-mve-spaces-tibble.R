context("test-make-mve-spaces-tibble.R")

test_that("make_mve_spaces_tibble creates tibble correctly", {
  # Test with age11_res example
  result <- make_mve_spaces_tibble(age11_res)

  # Check that it has the correct number of rows
  expect_equal(nrow(result), length(age11_res$rho_each_top_subset))

  # Check that rho values are correct
  expect_equal(result$rho, age11_res$rho_each_top_subset)

  # Check that all other columns are logical
  for (col in colnames(result)[-1]) {
    expect_true(is.logical(result[[col]]))
  }
})

test_that("make_mve_spaces_tibble respects order parameter", {
  # Test with custom order
  result_ordered <- make_mve_spaces_tibble(age11_res, order = c("pink"))

  # Check column order (pink columns should come before others)
  col_names <- colnames(result_ordered)
  pink_cols <- grep("^pink_", col_names)
  other_cols <- grep("^pink_", col_names, invert = TRUE)[-1]  # Exclude rho

  if (length(pink_cols) > 0 && length(other_cols) > 0) {
    expect_true(max(pink_cols) < min(other_cols))
  }
})
