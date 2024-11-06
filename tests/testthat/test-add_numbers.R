test_that("Basic check", {
  expect_equal(add_numbers(2, 3), 5)
  expect_equal(add_numbers(2, -3), -1)
  expect_equal(add_numbers(-2, 3), 1)
  expect_equal(add_numbers(2.5, 3), 5.5)
})

test_that("Check for non-numeric", {
  expect_error(add_numbers("a", 3), "`x` must be a number")
  expect_error(add_numbers(2, "b"), "`y` must be a number")
  expect_error(add_numbers(2, list(a = 20,
                                   b = 10)), "`y` must be a number")
})
