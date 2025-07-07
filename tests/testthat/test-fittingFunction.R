# --- Unit Tests for fittingFunction ---

test_that("fittingFunction calculates the squared difference correctly", {
  # Test with known values to check the basic calculation
  expect_equal(fittingFunction(lambda = 2, q = 0.5,
                               nu = 2, sigmaSquaredHat = 10),
               50.6176743)
  
  # Another test case
  expect_equal(fittingFunction(lambda = 5, q = 0.1,
                               nu = 4, sigmaSquaredHat = 5),
               5.9006287)
})

test_that("fittingFunction output is always non-negative", {
  # The result is squared so it should never be less than zero
  
  # Case where sigmaSquaredHat is larger
  result1 <- fittingFunction(lambda = 2, q = 0.5,
                             nu = 2, sigmaSquaredHat = 10)
  expect_gte(result1, 0) # gte means "greater than or equal to"
  
  # Case where sigmaSquaredHat is smaller
  result2 <- fittingFunction(lambda = 10, q = 0.1,
                             nu = 2, sigmaSquaredHat = 1)
  expect_gte(result2, 0)
})