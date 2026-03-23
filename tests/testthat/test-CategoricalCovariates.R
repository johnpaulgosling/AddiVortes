# --- Tests for categorical covariate support ---

test_that("AddiVortes fits with one binary categorical covariate (d = 2)", {
  skip_on_cran()
  withr::local_seed(42)
  n <- 30
  # Two numeric covariates plus one binary categorical (d = 2: "A", "B")
  x <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n),
    cat2 = sample(c("A", "B"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
    showProgress = FALSE
  )

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
  skip_on_cran()
  withr::local_seed(7)
  n <- 40
  # Two numeric covariates plus one categorical with 3 levels (d = 3)
  x <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n),
    cat3 = sample(c("low", "mid", "high"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
    showProgress = FALSE
  )

  expect_true(is.finite(fit$inSampleRmse))
  expect_false(is.null(fit$catEncoding))

  # A categorical variable with d = 3 produces d - 1 = 2 binary columns
  # Total encoded columns: 2 numeric + 2 binary = 4
  expect_equal(length(fit$xCentres), 4L)
})

test_that("AddiVortes fits with both a d=2 and a d=3 categorical covariate", {
  skip_on_cran()
  withr::local_seed(123)
  n <- 50
  x <- data.frame(
    x1 = rnorm(n),
    cat2 = sample(c("yes", "no"), n, replace = TRUE),
    x2 = runif(n),
    cat3 = sample(c("red", "green", "blue"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
    showProgress = FALSE
  )

  expect_true(is.finite(fit$inSampleRmse))
  expect_false(is.null(fit$catEncoding))

  # cat2 -> 1 binary column; cat3 -> 2 binary columns; 2 numeric -> 2 columns
  # Total encoded columns: 2 + 1 + 2 = 5
  expect_equal(length(fit$xCentres), 5L)

  # catColIndices: columns 2 and 4 in the original data frame
  expect_equal(sort(fit$catEncoding$catColIndices), c(2L, 4L))
})

test_that("predict works after fitting with categorical covariates", {
  skip_on_cran()
  withr::local_seed(99)
  n <- 40
  x <- data.frame(
    x1 = rnorm(n),
    cat2 = sample(c("A", "B"), n, replace = TRUE),
    cat3 = sample(c("p", "q", "r"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
    showProgress = FALSE
  )

  # Predict on new data with the same column structure
  x_new <- data.frame(
    x1 = rnorm(10),
    cat2 = sample(c("A", "B"), 10, replace = TRUE),
    cat3 = sample(c("p", "q", "r"), 10, replace = TRUE),
    stringsAsFactors = FALSE
  )

  preds <- predict(fit, x_new, showProgress = FALSE)
  expect_equal(length(preds), 10L)
  expect_true(all(is.finite(preds)))
})

test_that("catScaling controls the binary indicator value", {
  skip_on_cran()
  withr::local_seed(55)
  n <- 30
  x <- data.frame(
    x1 = rnorm(n),
    cat2 = sample(c("A", "B"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit_default <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
    catScaling = 1, showProgress = FALSE
  )
  fit_half <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
    catScaling = 0.5, showProgress = FALSE
  )

  # Both models should store their respective catScaling
  expect_equal(fit_default$catEncoding$catScaling, 1)
  expect_equal(fit_half$catEncoding$catScaling, 0.5)
})

test_that("factor columns are handled like character columns", {
  skip_on_cran()
  withr::local_seed(77)
  n <- 30
  x <- data.frame(
    x1 = rnorm(n),
    cat2 = factor(sample(c("A", "B"), n, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
    showProgress = FALSE
  )

  expect_true(is.finite(fit$inSampleRmse))
  # factor with d = 2 levels produces 1 binary column; total = 1 numeric + 1 binary = 2
  expect_equal(length(fit$xCentres), 2L)
})

test_that("unseen categories at predict time map to the reference level (all zeros)", {
  skip_on_cran()
  withr::local_seed(11)
  n <- 40
  x <- data.frame(
    x1 = rnorm(n),
    cat3 = sample(c("low", "mid", "high"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
    showProgress = FALSE
  )

  # "unknown" is not a level seen during training; it should map to all zeros
  # (i.e., treated as the reference level)
  x_new <- data.frame(
    x1 = c(0.5, -0.5),
    cat3 = c("unknown", "low"), # "unknown" is unseen; "high" is the reference
    stringsAsFactors = FALSE
  )
  # Both rows should yield finite predictions (not errors or NAs)
  preds <- predict(fit, x_new, showProgress = FALSE)
  expect_equal(length(preds), 2L)
  expect_true(all(is.finite(preds)))

  # Encoding the "unknown" row should produce the same binary columns as the
  # reference level "high" (alphabetically first, so dropped; all indicator = 0)
  enc <- encodeCategories_internal(x_new, encoding = fit$catEncoding)$encoded
  # Levels sorted alphabetically: "high" (ref, dropped), "low", "mid"
  # Binary columns: cat3_low, cat3_mid
  # Row 1: "unknown" (unseen) -> treated as reference "high" -> all zeros
  binary_cols <- fit$catEncoding$encodedBinaryCols
  expect_true(all(enc[1, binary_cols] == 0))
  # Row 2: "low" -> cat3_low = 1 (non-reference), cat3_mid = 0
  expect_equal(unname(enc[2, binary_cols]), c(1, 0))
})

test_that("categorical-only covariates (no numeric columns) work", {
  skip_on_cran()
  withr::local_seed(22)
  n <- 40
  x <- data.frame(
    grp2 = sample(c("A", "B"), n, replace = TRUE),
    grp3 = sample(c("x", "y", "z"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
    showProgress = FALSE
  )

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


test_that("encodeCategories_internal produces correct binary values", {
  x <- data.frame(
    num = c(1, 2, 3, 4),
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
  # enc[row, col_name] returns a named scalar in R; ignore_attr = TRUE
  # compares values only, ignoring the column name attached to the result.
  expect_equal(enc[1, "cat3_low"], 0, ignore_attr = TRUE) # "mid" row
  expect_equal(enc[1, "cat3_mid"], 1, ignore_attr = TRUE) # "mid" row
  expect_equal(enc[2, "cat3_low"], 0, ignore_attr = TRUE) # "high" row (reference)
  expect_equal(enc[2, "cat3_mid"], 0, ignore_attr = TRUE) # "high" row (reference)
  expect_equal(enc[3, "cat3_low"], 1, ignore_attr = TRUE) # "low" row
  expect_equal(enc[3, "cat3_mid"], 0, ignore_attr = TRUE) # "low" row
})

test_that("catScaling must be a single positive number", {
  skip_on_cran()
  withr::local_seed(1)
  x <- data.frame(x1 = rnorm(20), cat2 = sample(c("A", "B"), 20, replace = TRUE))
  y <- rnorm(20)

  expect_error(
    AddiVortes(y, x,
      m = 3, totalMCMCIter = 10, mcmcBurnIn = 3,
      catScaling = 0, showProgress = FALSE
    ),
    "must be a single positive number"
  )
  expect_error(
    AddiVortes(y, x,
      m = 3, totalMCMCIter = 10, mcmcBurnIn = 3,
      catScaling = -1, showProgress = FALSE
    ),
    "must be a single positive number"
  )
  expect_error(
    AddiVortes(y, x,
      m = 3, totalMCMCIter = 10, mcmcBurnIn = 3,
      catScaling = c(1, 2), showProgress = FALSE
    ),
    "must be a single positive number"
  )
})


test_that("catEncoding stores correct origNCols and origColNames", {
  skip_on_cran()
  withr::local_seed(66)
  n <- 30
  x <- data.frame(
    x1 = rnorm(n),
    grp = sample(c("a", "b", "c"), n, replace = TRUE),
    x2 = runif(n),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
    showProgress = FALSE
  )

  expect_equal(fit$catEncoding$origNCols, 3L)
  expect_equal(fit$catEncoding$origColNames, c("x1", "grp", "x2"))
  # grp (d=3) -> 2 binary cols; total encoded = 1 + 2 + 1 = 4
  expect_equal(length(fit$xCentres), 4L)
})

test_that("all stored tess centres for binary dimensions lie in [0, catScaling]", {
  skip_on_cran()
  withr::local_seed(88)
  n <- 40
  cs <- 0.75
  x <- data.frame(
    x1 = rnorm(n),
    cat3 = sample(c("a", "b", "c"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x,
    m = 8, totalMCMCIter = 60, mcmcBurnIn = 20,
    catScaling = cs, showProgress = FALSE
  )

  binaryCols <- fit$catEncoding$encodedBinaryCols

  # Collect all binary-dimension centre values across every posterior sample
  all_binary_vals <- unlist(lapply(seq_along(fit$posteriorTess), function(s) {
    mapply(function(tess_mat, dim_vec) {
      local_bin <- which(dim_vec %in% binaryCols)
      if (length(local_bin) > 0) as.vector(tess_mat[, local_bin]) else numeric(0)
    }, fit$posteriorTess[[s]], fit$posteriorDim[[s]], SIMPLIFY = FALSE)
  }))

  expect_true(all(all_binary_vals >= 0 & all_binary_vals <= cs))
})
