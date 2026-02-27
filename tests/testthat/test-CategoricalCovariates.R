# --- Tests for categorical covariate support ---

test_that("AddiVortes fits with one binary categorical covariate (d = 2)", {
  set.seed(42)
  n <- 30
  # Two numeric covariates plus one binary categorical (d = 2: "A", "B")
  x <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n),
    cat2 = sample(c("A", "B"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  # Model should fit without error and produce finite RMSE
  expect_true(is.finite(fit$inSampleRmse))

  # catEncoding should be stored and reference the correct column
  expect_false(is.null(fit$catEncoding))
  expect_equal(fit$catEncoding$catColIndices, 3L)

  # A categorical variable with d = 2 produces d - 1 = 1 binary column
  # Total encoded columns: 2 numeric + 1 binary = 3
  expect_equal(length(fit$xCentres), 3L)
})

test_that("AddiVortes fits with one 3-level categorical covariate (d = 3)", {
  set.seed(7)
  n <- 40
  # Two numeric covariates plus one categorical with 3 levels (d = 3)
  x <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n),
    cat3 = sample(c("low", "mid", "high"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  expect_true(is.finite(fit$inSampleRmse))
  expect_false(is.null(fit$catEncoding))

  # A categorical variable with d = 3 produces d - 1 = 2 binary columns
  # Total encoded columns: 2 numeric + 2 binary = 4
  expect_equal(length(fit$xCentres), 4L)
})

test_that("AddiVortes fits with both a d=2 and a d=3 categorical covariate", {
  set.seed(123)
  n <- 50
  x <- data.frame(
    x1    = rnorm(n),
    cat2  = sample(c("yes", "no"), n, replace = TRUE),
    x2    = runif(n),
    cat3  = sample(c("red", "green", "blue"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  expect_true(is.finite(fit$inSampleRmse))
  expect_false(is.null(fit$catEncoding))

  # cat2 -> 1 binary column; cat3 -> 2 binary columns; 2 numeric -> 2 columns
  # Total encoded columns: 2 + 1 + 2 = 5
  expect_equal(length(fit$xCentres), 5L)

  # catColIndices: columns 2 and 4 in the original data frame
  expect_equal(sort(fit$catEncoding$catColIndices), c(2L, 4L))
})

test_that("predict works after fitting with categorical covariates", {
  set.seed(99)
  n <- 40
  x <- data.frame(
    x1   = rnorm(n),
    cat2 = sample(c("A", "B"), n, replace = TRUE),
    cat3 = sample(c("p", "q", "r"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  # Predict on new data with the same column structure
  x_new <- data.frame(
    x1   = rnorm(10),
    cat2 = sample(c("A", "B"), 10, replace = TRUE),
    cat3 = sample(c("p", "q", "r"), 10, replace = TRUE),
    stringsAsFactors = FALSE
  )

  preds <- predict(fit, x_new, showProgress = FALSE)
  expect_equal(length(preds), 10L)
  expect_true(all(is.finite(preds)))
})

test_that("catScaling controls the binary indicator value", {
  set.seed(55)
  n <- 30
  x <- data.frame(
    x1   = rnorm(n),
    cat2 = sample(c("A", "B"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit_default <- AddiVortes(y, x, m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
                             catScaling = 1, showProgress = FALSE)
  fit_half    <- AddiVortes(y, x, m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
                             catScaling = 0.5, showProgress = FALSE)

  # Both models should store their respective catScaling
  expect_equal(fit_default$catEncoding$catScaling, 1)
  expect_equal(fit_half$catEncoding$catScaling, 0.5)
})

test_that("factor columns are handled like character columns", {
  set.seed(77)
  n <- 30
  x <- data.frame(
    x1   = rnorm(n),
    cat2 = factor(sample(c("A", "B"), n, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
                    showProgress = FALSE)

  expect_true(is.finite(fit$inSampleRmse))
  # factor with d = 2 levels produces 1 binary column; total = 1 numeric + 1 binary = 2
  expect_equal(length(fit$xCentres), 2L)
})

test_that("unseen categories at predict time map to the reference level (all zeros)", {
  set.seed(11)
  n <- 40
  x <- data.frame(
    x1   = rnorm(n),
    cat3 = sample(c("low", "mid", "high"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  # "unknown" is not a level seen during training; it should map to all zeros
  # (i.e., treated as the reference level)
  x_new <- data.frame(
    x1   = c(0.5, -0.5),
    cat3 = c("unknown", "low"),   # "unknown" is unseen; "low" is the reference
    stringsAsFactors = FALSE
  )
  # Both rows should yield finite predictions (not errors or NAs)
  preds <- predict(fit, x_new, showProgress = FALSE)
  expect_equal(length(preds), 2L)
  expect_true(all(is.finite(preds)))

  # Encoding the "unknown" row should produce the same binary columns as the
  # reference level "low" (all indicator columns = 0)
  enc <- encodeCategories_internal(x_new, encoding = fit$catEncoding)$encoded
  # Binary columns (for "mid" and "high") should both be 0 for the "unknown" row
  binary_cols <- fit$catEncoding$encodedBinaryCols
  expect_true(all(enc[1, binary_cols] == 0))
  expect_true(all(enc[2, binary_cols] == 0))
})

test_that("categorical-only covariates (no numeric columns) work", {
  set.seed(22)
  n <- 40
  x <- data.frame(
    grp2 = sample(c("A", "B"), n, replace = TRUE),
    grp3 = sample(c("x", "y", "z"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  expect_true(is.finite(fit$inSampleRmse))
  # grp2 -> 1 binary col, grp3 -> 2 binary cols; total = 3 encoded cols
  expect_equal(length(fit$xCentres), 3L)

  x_new <- data.frame(
    grp2 = c("A", "B"),
    grp3 = c("x", "z"),
    stringsAsFactors = FALSE
  )
  preds <- predict(fit, x_new, showProgress = FALSE)
  expect_equal(length(preds), 2L)
  expect_true(all(is.finite(preds)))
})

test_that("4-level categorical produces 3 binary columns", {
  set.seed(33)
  n <- 50
  x <- data.frame(
    x1   = rnorm(n),
    cat4 = sample(c("Q1", "Q2", "Q3", "Q4"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  expect_true(is.finite(fit$inSampleRmse))
  # cat4 (d=4) -> 3 binary cols; 1 numeric -> 1 col; total = 4
  expect_equal(length(fit$xCentres), 4L)
  expect_equal(fit$catEncoding$catColIndices, 2L)
})

test_that("encodeCategories_internal produces correct binary values", {
  x <- data.frame(
    num  = c(1, 2, 3, 4),
    cat2 = c("A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  result <- encodeCategories_internal(x, catScaling = 2)
  enc <- result$encoded

  # Column layout: numeric 'num', then binary 'cat2_B' (reference = "A")
  expect_equal(ncol(enc), 2L)
  expect_equal(colnames(enc), c("num", "cat2_B"))

  # "A" rows should have cat2_B = 0; "B" rows should have cat2_B = 2
  expect_equal(enc[, "cat2_B"], c(0, 2, 0, 2))
  # Numeric column should be unchanged
  expect_equal(enc[, "num"], c(1, 2, 3, 4))
})

test_that("encodeCategories_internal applies stored encoding to new data", {
  x_train <- data.frame(
    cat3 = c("low", "mid", "high", "low", "mid"),
    stringsAsFactors = FALSE
  )
  result_train <- encodeCategories_internal(x_train, catScaling = 1)
  encoding <- result_train$encoding

  x_new <- data.frame(cat3 = c("mid", "high", "low"), stringsAsFactors = FALSE)
  result_new <- encodeCategories_internal(x_new, encoding = encoding)
  enc <- result_new$encoded

  # Levels sorted alphabetically: "high" (ref, dropped), "low", "mid"
  # Binary cols: cat3_low, cat3_mid
  expect_equal(colnames(enc), c("cat3_low", "cat3_mid"))
  expect_equal(enc[1, "cat3_low"], 0)   # "mid" row
  expect_equal(enc[1, "cat3_mid"], 1)   # "mid" row
  expect_equal(enc[2, "cat3_low"], 0)   # "high" row (reference)
  expect_equal(enc[2, "cat3_mid"], 0)   # "high" row (reference)
  expect_equal(enc[3, "cat3_low"], 1)   # "low" row
  expect_equal(enc[3, "cat3_mid"], 0)   # "low" row
})

test_that("catScaling must be a single positive number", {
  set.seed(1)
  x <- data.frame(x1 = rnorm(20), cat2 = sample(c("A", "B"), 20, replace = TRUE))
  y <- rnorm(20)

  expect_error(
    AddiVortes(y, x, m = 3, totalMCMCIter = 10, mcmcBurnIn = 3,
               catScaling = 0, showProgress = FALSE),
    "must be a single positive number"
  )
  expect_error(
    AddiVortes(y, x, m = 3, totalMCMCIter = 10, mcmcBurnIn = 3,
               catScaling = -1, showProgress = FALSE),
    "must be a single positive number"
  )
  expect_error(
    AddiVortes(y, x, m = 3, totalMCMCIter = 10, mcmcBurnIn = 3,
               catScaling = c(1, 2), showProgress = FALSE),
    "must be a single positive number"
  )
})

test_that("predict.AddiVortes accepts data.frame for model with categorical covariates", {
  set.seed(44)
  n <- 40
  x <- data.frame(
    x1   = rnorm(n),
    cat2 = sample(c("A", "B"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  x_new_df <- data.frame(x1 = rnorm(5), cat2 = c("A", "B", "A", "B", "A"),
                          stringsAsFactors = FALSE)

  # Should not error — data.frame is accepted when model has catEncoding
  preds <- predict(fit, x_new_df, showProgress = FALSE)
  expect_equal(length(preds), 5L)
  expect_true(all(is.finite(preds)))
})

test_that("catEncoding stores correct origNCols and origColNames", {
  set.seed(66)
  n <- 30
  x <- data.frame(
    x1   = rnorm(n),
    grp  = sample(c("a", "b", "c"), n, replace = TRUE),
    x2   = runif(n),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
                    showProgress = FALSE)

  expect_equal(fit$catEncoding$origNCols, 3L)
  expect_equal(fit$catEncoding$origColNames, c("x1", "grp", "x2"))
  # grp (d=3) -> 2 binary cols; total encoded = 1 + 2 + 1 = 4
  expect_equal(length(fit$xCentres), 4L)
})
