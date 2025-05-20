# --- Test Suite for sample_mu_values ---

set.seed(123211)

test_that("Simple Tesselation", {
  Tess <- list(matrix(0,1,1))

  results <- Sample_mu_values(1, Tess,
                              1, 3,
                              1, 1)

  expect_equal(results, 0.146,
               tolerance = 0.01)
})
