# --- Test Suite for brute force allocation ---

test_that("Brute force allocation works correctly", {
  # Test 1: Simple case with clear nearest neighbors
  set.seed(123)
  centers <- matrix(c(0, 0, 
                      1, 0, 
                      0, 1), 
                    ncol = 2, byrow = TRUE)
  obs <- matrix(c(0.1, 0.1,   # closest to (0,0) -> index 1
                  0.9, 0.1,   # closest to (1,0) -> index 2
                  0.1, 0.9),  # closest to (0,1) -> index 3
                ncol = 2, byrow = TRUE)
  
  result <- AddiVortes:::assign_bruteforce(obs, centers)
  
  expect_equal(result, c(1, 2, 3))
})

test_that("cellIndices with single center", {
  set.seed(456)
  x <- matrix(rnorm(50), 10, 5)
  tess <- matrix(c(0, 0), nrow = 1, ncol = 2)
  dim <- c(1, 2)
  
  result <- AddiVortes:::cellIndices(x, tess, dim)
  
  # All observations should be assigned to the single center
  expect_equal(result, rep(1, 10))
})

test_that("cellIndices with multiple centers", {
  set.seed(789)
  # Create observations
  x <- matrix(c(0.1, 0.1, 0, 0, 0,
                0.9, 0.1, 0, 0, 0,
                0.1, 0.9, 0, 0, 0), 
              ncol = 5, byrow = TRUE)
  
  # Create centers in first 2 dimensions
  tess <- matrix(c(0, 0, 
                   1, 0, 
                   0, 1), 
                 ncol = 2, byrow = TRUE)
  dim <- c(1, 2)
  
  result <- AddiVortes:::cellIndices(x, tess, dim)
  
  expect_equal(result, c(1, 2, 3))
})

test_that("Brute force handles edge cases", {
  # Test with identical points
  centers <- matrix(c(1, 1,
                      2, 2), 
                    ncol = 2, byrow = TRUE)
  obs <- matrix(c(1.5, 1.5), ncol = 2)
  
  result <- AddiVortes:::assign_bruteforce(obs, centers)
  
  # Should assign to one of the centers (either 1 or 2)
  expect_true(result[1] %in% c(1, 2))
  expect_length(result, 1)
})
