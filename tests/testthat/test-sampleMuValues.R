# --- Test Suite for sample_mu_values ---

test_that("Simple Tesselation", {
  withr::local_seed(123211)
  tess <- list(matrix(0, 1, 1))

  results <- sampleMuValues(
    1, tess,
    1, 3,
    1, 1
  )

  expect_equal(results, 0.146,
               tolerance = 0.01
  )
})
