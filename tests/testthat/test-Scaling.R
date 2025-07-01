# --- Tests for scaleData_internal ---

test_that("scaleData_internal scales a numeric vector correctly", {
  vec <- c(10, 20, 30, 40)
  result <- scaleData_internal(vec)

  # Expected: centre = 25, range = 30
  # (10-25)/30 = -0.5
  # (40-25)/30 = 0.5
  expect_equal(result$scaledData, c(-0.5, -1/6, 1/6, 0.5))
  expect_equal(result$centres, 25)
  expect_equal(result$ranges, 30)
  expect_type(result, "list")
})

test_that("scaleData_internal scales a numeric matrix column-wise", {
  mat <- matrix(c(10, 20, 100, 200), nrow = 2) # Two columns
  result <- scaleData_internal(mat)

  # Col 1: centre=15, range=10 -> scaled: -0.5, 0.5
  # Col 2: centre=150, range=100 -> scaled: -0.5, 0.5
  expectedScaled <- matrix(c(-0.5, 0.5, -0.5, 0.5), nrow = 2)
  expect_equal(result$centres, c(15, 150))
  expect_equal(result$ranges, c(10, 100))
})

test_that("scaleData_internal throws an error for zero-range data", {
  vec <- rep(5, 4)
  expect_error(scaleData_internal(vec), "Range is zero. Cannot scale data.")
})


# --- Tests for applyScaling_internal ---

test_that("applyScaling_internal applies scaling correctly", {
  # Get scaling parameters from a "training" set
  trainMat <- matrix(c(0, 10, 100, 200), nrow = 2)
  scalingParams <- scaleData_internal(trainMat)

  # Apply to a new "test" set
  testMat <- matrix(c(5, 15, 50, 250), nrow = 2)
  scaledTestMat <- applyScaling_internal(testMat,
                                         scalingParams$centres,
                                         scalingParams$ranges)

  # Expected for testMat using trainMat params:
  # Col 1 (centre=5, range=10): (5-5)/10=0, (15-5)/10=1
  # Col 2 (centre=150, range=100): (50-150)/100=-1, (250-150)/100=1
  expectedScaled <- matrix(c(0, 1, -1, 1), nrow = 2)
  expect_true(scaledTestMat[1,1] == expectedScaled[1,1])
})

test_that("applyScaling_internal warns for NA parameters", {
  mat <- matrix(1:4, nrow = 2)
  centres <- c(1, NA)
  ranges <- c(1, NA)

  expect_warning(
    applyScaling_internal(mat, centres, ranges),
    "Scaling parameters for column 2"
  )
})
