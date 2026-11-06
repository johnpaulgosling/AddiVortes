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

  # Categorical columns are moved to the end for membership grouping, so
  # catColIndices refer to the reordered frame (x1, x2, cat2, cat3).
  expect_equal(sort(fit$catEncoding$catColIndices), c(3L, 4L))
  expect_equal(
    fit$catEncoding$origColNames,
    c("x1", "x2", "cat2", "cat3")
  )
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

test_that("covariateStructure assigns distinct membership for auto-detected categories", {
  x <- data.frame(
    age = rnorm(10),
    income = runif(10),
    region = factor(sample(c("East", "North", "South", "West"), 10, replace = TRUE)),
    product = factor(sample(c("Basic", "Deluxe", "Premium"), 10, replace = TRUE))
  )

  san <- AddiVortes:::covariateStructure_internal(
    x, structure = rep("E", 4), membership = NULL, one.hot = FALSE
  )

  expect_equal(san$structure, c("E", "E", "C", "C"))
  # Euclidean columns share a group; each categorical has its own group
  expect_equal(san$membership[1], san$membership[2])
  expect_true(san$membership[3] != san$membership[1])
  expect_true(san$membership[4] != san$membership[1])
  expect_true(san$membership[3] != san$membership[4])
  expect_equal(colnames(san$data), c("age", "income", "region", "product"))
})

test_that("AddiVortes fits with Eskin distance when cat.onehot is FALSE", {
  skip_on_cran()
  withr::local_seed(87)
  n <- 40
  x <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n),
    cat2 = factor(sample(c("A", "B"), n, replace = TRUE)),
    cat3 = factor(sample(c("low", "mid", "high"), n, replace = TRUE),
                  levels = c("high", "low", "mid"))
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x,
    m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
    cat.onehot = FALSE, showProgress = FALSE
  )

  expect_true(is.finite(fit$inSampleRmse))
  # No one-hot expansion: still 4 covariate columns
  expect_equal(length(fit$xCentres), 4L)
  expect_null(fit$catEncoding)

  x_new <- x[1:5, ]
  preds <- predict(fit, x_new, showProgress = FALSE)
  expect_equal(length(preds), 5L)
  expect_true(all(is.finite(preds)))
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
  # Categorical columns are moved to the end before encoding
  expect_equal(fit$catEncoding$origColNames, c("x1", "x2", "grp"))
  # grp (d=3) -> 2 binary cols; total encoded = 1 + 1 + 2 = 4
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

# Nearest centre under Eskin distance. Column layout is given by the reduced
# metric: one group per entry of `metric`, with sizes in `members`.
eskin_nearest <- function(query, centres, metric, members, ncat = NULL, weight = NULL) {
  storage.mode(query) <- "double"
  storage.mode(centres) <- "double"
  AddiVortes:::cellIndices(
    query,
    centres,
    dim = seq_len(ncol(query)),
    metric = as.integer(metric),
    members = as.integer(members),
    ncat = ncat,
    eskinWeight = weight
  )
}

# Bracket one categorical mismatch against a Euclidean rival.
# Query is (0, category 1). Centre 1 mismatches the category at the same
# Euclidean coordinate. Centre 2 matches the category.
expect_eskin_cost <- function(ncat, weight, expected) {
  query <- matrix(c(0, 1), nrow = 1)
  below <- matrix(c(0, 2, sqrt(expected - 1e-8), 1), nrow = 2, byrow = TRUE)
  above <- matrix(c(0, 2, sqrt(expected + 1e-8), 1), nrow = 2, byrow = TRUE)
  expect_equal(
    eskin_nearest(query, below, c(0L, 2L), c(1L, 1L), ncat, weight),
    2L
  )
  expect_equal(
    eskin_nearest(query, above, c(0L, 2L), c(1L, 1L), ncat, weight),
    1L
  )
}

expect_eskin_match_is_free <- function(ncat, weight) {
  query <- matrix(c(0, 1), nrow = 1)
  centres <- matrix(c(sqrt(1e-6), 1, 0, 2), nrow = 2, byrow = TRUE)
  expect_equal(
    eskin_nearest(query, centres, c(0L, 2L), c(1L, 1L), ncat, weight),
    1L
  )
}

minimal_eskin_fit <- function(metric_aug, member_aug, eskin, tess, mu) {
  p <- length(metric_aug)
  new_AddiVortes(
    posteriorTess = list(list(tess)),
    posteriorDim = list(list(seq_len(ncol(tess)))),
    posteriorSigma = 1,
    posteriorPred = list(list(as.numeric(mu))),
    xCentres = rep(0, p),
    xRanges = rep(1, p),
    yCentre = 0,
    yRange = 1,
    inSampleRmse = NA_real_,
    metric = rep(0L, p),
    members = NULL,
    metric_aug = as.integer(metric_aug),
    member_aug = as.integer(member_aug),
    eskin = eskin
  )
}

test_that("default Eskin weight charges 2/d^2 and a match charges 0", {
  expect_eskin_cost(2L, 1, 0.5)
  expect_eskin_cost(4L, 1, 0.125)
  expect_eskin_cost(24L, 1, 2 / 576)
  expect_eskin_match_is_free(2L, 1)
  expect_eskin_match_is_free(4L, 1)
  expect_eskin_match_is_free(24L, 1)
  # Omitting the weight is the same as weight 1.
  expect_eskin_cost(2L, NULL, 0.5)
})

test_that("Eskin weight d^2/2 makes every mismatch contribute 1", {
  expect_eskin_cost(2L, 2^2 / 2, 1)
  expect_eskin_cost(4L, 4^2 / 2, 1)
  expect_eskin_cost(24L, 24^2 / 2, 1)
  expect_eskin_match_is_free(2L, 2)
  expect_eskin_match_is_free(4L, 8)
  expect_eskin_match_is_free(24L, 288)
})

test_that("two factors use the weight that belongs to each column", {
  query <- matrix(c(1, 1), nrow = 1)
  centres <- matrix(c(2, 1, 1, 2), nrow = 2, byrow = TRUE)
  metric <- c(2L, 2L)
  members <- c(1L, 1L)
  # Column 1 (d = 2) costs w1 * 0.5. Column 2 (d = 4) costs w2 * 0.125.
  expect_equal(
    eskin_nearest(query, centres, metric, members, c(2L, 4L), c(1, 1)),
    2L
  )
  expect_equal(
    eskin_nearest(query, centres, metric, members, c(2L, 4L), c(1, 32)),
    1L
  )
  # A scalar is recycled.
  expect_equal(
    eskin_nearest(query, centres, metric, members, c(2L, 4L), 8),
    2L
  )
})

test_that("predict uses the training level count and the stored weight", {
  levs <- paste0("L", 1:4)
  x_new <- data.frame(
    cont = c(0, 0),
    grp = factor(c("L1", "L2"), levels = levs)
  )
  tess <- matrix(c(0, 2, sqrt(0.3), 1), nrow = 2, byrow = TRUE)
  fit <- minimal_eskin_fit(
    c(0L, 2L),
    c(1L, 2L),
    list(columns = "grp", ncat = 4L, weight = 1),
    tess,
    c(0, 1)
  )

  # Training d = 4 costs 0.125, so the mismatch centre wins for L1.
  # L2 matches that centre. Recomputing d from this hold-out would use d = 2
  # and cost 0.5, which would send L1 to the other centre.
  expect_equal(predict(fit, x_new, showProgress = FALSE), c(0, 0))

  fit_recomputed <- fit
  fit_recomputed$eskin <- NULL
  expect_equal(predict(fit_recomputed, x_new, showProgress = FALSE), c(1, 0))

  fit$eskin$weight <- 4^2 / 2
  expect_equal(predict(fit, x_new, showProgress = FALSE), c(1, 0))
})

test_that("predict applies a different Eskin weight to each column", {
  x_new <- data.frame(
    a = factor("A", levels = c("A", "B")),
    b = factor("A", levels = c("A", "B"))
  )
  tess <- matrix(c(2, 1, 1, 2), nrow = 2, byrow = TRUE)
  fit <- minimal_eskin_fit(
    c(2L, 2L),
    c(1L, 2L),
    list(columns = c("a", "b"), ncat = c(2L, 2L), weight = c(1, 100)),
    tess,
    c(10, 20)
  )
  expect_equal(predict(fit, x_new, showProgress = FALSE), 10)

  fit$eskin$weight <- c(100, 1)
  expect_equal(predict(fit, x_new, showProgress = FALSE), 20)
})

test_that("AddiVortes stores training Eskin level counts and weights", {
  skip_on_cran()
  n <- 12
  x <- data.frame(
    z = rnorm(n),
    b = factor(rep(letters[1:4], length.out = n), levels = letters[1:4]),
    a = factor(rep(c("A", "B"), length.out = n), levels = c("A", "B"))
  )
  y <- rnorm(n)

  fit <- AddiVortes(
    y, x,
    m = 2, totalMCMCIter = 4, mcmcBurnIn = 1,
    cat.onehot = FALSE, eskinWeight = c(0.5, 3),
    showProgress = FALSE
  )
  expect_equal(fit$eskin$columns, c("b", "a"))
  expect_equal(fit$eskin$ncat, c(4L, 2L))
  expect_equal(fit$eskin$weight, c(0.5, 3))
  expect_null(fit$catEncoding)

  fit_default <- AddiVortes(
    y, x,
    m = 2, totalMCMCIter = 4, mcmcBurnIn = 1,
    cat.onehot = FALSE, showProgress = FALSE
  )
  expect_equal(fit_default$eskin$weight, c(1, 1))

  # The stored count is the largest training code, not an unused tail level.
  x_partial <- data.frame(
    grp = factor(rep(c("A", "B"), length.out = n), levels = c("A", "B", "C", "D"))
  )
  fit_partial <- AddiVortes(
    rnorm(n), x_partial,
    m = 2, totalMCMCIter = 4, mcmcBurnIn = 1,
    cat.onehot = FALSE, showProgress = FALSE
  )
  expect_equal(fit_partial$eskin$ncat, 2L)
  expect_equal(fit_partial$eskin$columns, "grp")
})

test_that("eskinWeight rejects non-positive values and the wrong length", {
  x <- data.frame(z = rnorm(8), g = factor(rep(c("A", "B"), 4)))
  y <- rnorm(8)
  expect_error(
    AddiVortes(
      y, x,
      m = 2, totalMCMCIter = 2, mcmcBurnIn = 0,
      cat.onehot = FALSE, eskinWeight = 0, showProgress = FALSE
    ),
    "positive finite"
  )
  expect_error(
    AddiVortes(
      y, x,
      m = 2, totalMCMCIter = 2, mcmcBurnIn = 0,
      cat.onehot = FALSE, eskinWeight = -1, showProgress = FALSE
    ),
    "positive finite"
  )
  expect_error(
    AddiVortes(
      y, x,
      m = 2, totalMCMCIter = 2, mcmcBurnIn = 0,
      cat.onehot = FALSE, eskinWeight = NA_real_, showProgress = FALSE
    ),
    "positive finite"
  )
  expect_error(
    AddiVortes(
      y, x,
      m = 2, totalMCMCIter = 2, mcmcBurnIn = 0,
      cat.onehot = FALSE, eskinWeight = c(1, 2), showProgress = FALSE
    ),
    "one value per categorical column"
  )
})

test_that("one-hot and Euclidean fits ignore eskinWeight", {
  skip_on_cran()
  x <- data.frame(z = rnorm(20), g = factor(rep(c("A", "B"), 10)))
  y <- rnorm(20)

  fit_default <- withr::with_seed(7, AddiVortes(
    y, x,
    m = 3, totalMCMCIter = 12, mcmcBurnIn = 4, showProgress = FALSE
  ))
  fit_ignored <- withr::with_seed(7, AddiVortes(
    y, x,
    m = 3, totalMCMCIter = 12, mcmcBurnIn = 4,
    eskinWeight = 9, showProgress = FALSE
  ))
  expect_null(fit_ignored$eskin)
  expect_identical(fit_default$posteriorTess, fit_ignored$posteriorTess)
  expect_identical(fit_default$posteriorPred, fit_ignored$posteriorPred)
  expect_equal(fit_default$inSampleRmse, fit_ignored$inSampleRmse)
  expect_no_error(withr::with_seed(7, AddiVortes(
    y, x,
    m = 3, totalMCMCIter = 12, mcmcBurnIn = 4,
    eskinWeight = -1, showProgress = FALSE
  )))

  xm <- matrix(rnorm(40), 20, 2)
  ym <- rnorm(20)
  fit_euc <- withr::with_seed(3, AddiVortes(
    ym, xm,
    m = 3, totalMCMCIter = 10, mcmcBurnIn = 3, showProgress = FALSE
  ))
  fit_euc_weight <- withr::with_seed(3, AddiVortes(
    ym, xm,
    m = 3, totalMCMCIter = 10, mcmcBurnIn = 3,
    eskinWeight = 5, showProgress = FALSE
  ))
  expect_null(fit_euc_weight$eskin)
  expect_identical(fit_euc$posteriorTess, fit_euc_weight$posteriorTess)
  expect_identical(fit_euc$posteriorPred, fit_euc_weight$posteriorPred)
})

test_that("an explicit Eskin weight of 1 matches the default", {
  skip_on_cran()
  x <- data.frame(
    z = rnorm(20),
    g = factor(rep(c("A", "B", "C", "D"), 5))
  )
  y <- rnorm(20)
  fit_default <- withr::with_seed(11, AddiVortes(
    y, x,
    m = 3, totalMCMCIter = 12, mcmcBurnIn = 4,
    cat.onehot = FALSE, showProgress = FALSE
  ))
  fit_one <- withr::with_seed(11, AddiVortes(
    y, x,
    m = 3, totalMCMCIter = 12, mcmcBurnIn = 4,
    cat.onehot = FALSE, eskinWeight = 1, showProgress = FALSE
  ))
  expect_identical(fit_default$posteriorTess, fit_one$posteriorTess)
  expect_identical(fit_default$posteriorDim, fit_one$posteriorDim)
  expect_identical(fit_default$posteriorPred, fit_one$posteriorPred)
  expect_equal(fit_default$inSampleRmse, fit_one$inSampleRmse)
  expect_equal(fit_one$eskin$ncat, 4L)
  expect_equal(fit_one$eskin$weight, 1)
})

test_that("classification predict reuses the stored Eskin specification", {
  skip_on_cran()
  x <- data.frame(
    z = rnorm(30),
    g = factor(rep(c("A", "B"), 15))
  )
  y <- factor(rep(c("no", "yes"), 15))
  fit <- AddiVortes(
    y, x,
    m = 2, totalMCMCIter = 8, mcmcBurnIn = 2,
    cat.onehot = FALSE, showProgress = FALSE
  )
  expect_equal(fit$eskin$ncat, 2L)
  expect_equal(fit$eskin$weight, 1)
  preds <- predict(fit, x[1:4, ], type = "class", showProgress = FALSE)
  expect_length(preds, 4L)
  expect_true(all(preds %in% c("no", "yes")))
})
