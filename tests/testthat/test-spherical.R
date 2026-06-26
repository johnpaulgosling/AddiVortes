test_that("covariateStructure places longitude last within a spherical group", {
  lat <- c(-pi / 4, pi / 4, 0)
  lon <- c(-pi, pi, 0)
  df <- data.frame(longitude = lon, colatitude = lat)

  san <- covariateStructure_internal(df, c("S", "S"))
  extents <- apply(san$data, 2, function(col) diff(range(col)))

  expect_lte(extents[1], pi)
  expect_gt(extents[2], pi)
  expect_equal(san$data[, 1], lat)
  expect_equal(san$data[, 2], lon)
  expect_equal(names(san$data), c("colatitude", "longitude"))
})

test_that("covariateStructure keeps names when spherical order is already correct", {
  lat <- c(-pi / 4, pi / 4, 0)
  lon <- c(-pi, pi, 0)
  df <- data.frame(colatitude = lat, longitude = lon)

  san <- covariateStructure_internal(df, c("S", "S"))

  expect_equal(names(san$data), c("colatitude", "longitude"))
  expect_equal(san$data$colatitude, lat)
  expect_equal(san$data$longitude, lon)
})

test_that("fit preprocessing leaves spherical columns in original radians", {
  lat <- c(-0.3, 0.4, 0.1)
  lon <- c(-2.5, 2.8, 0.5)
  x <- cbind(lat, lon)
  metric <- c(1L, 1L)

  xScalingResult <- scaleData_internal(x)
  xScaled <- xScalingResult$scaledData
  xScaled[, metric != 0] <- x[, metric != 0]

  expect_equal(xScaled[, 1], lat)
  expect_equal(xScaled[, 2], lon)
  expect_false(isTRUE(all.equal(xScalingResult$scaledData[, 1], lat)))
})

test_that("predict preprocessing matches fit for spherical columns", {
  lat <- c(-0.3, 0.4)
  lon <- c(-2.5, 2.8)
  newdata <- cbind(lat, lon)
  metric <- c(1L, 1L)
  centres <- c(mean(lat), mean(lon))
  ranges <- c(diff(range(lat)), diff(range(lon)))

  xNewScaled <- applyScaling_internal(newdata, centres, ranges)
  xNewScaled[, metric != 0] <- newdata[, metric != 0]

  expect_equal(xNewScaled[, 1], lat)
  expect_equal(xNewScaled[, 2], lon)
})

test_that("spherical fit and predict preserve coordinate values", {
  skip_on_cran()
  withr::local_seed(42)

  n <- 10
  lat <- runif(n, -pi / 4, pi / 4)
  lon <- runif(n, -pi, pi)
  x <- data.frame(lat = lat, lon = lon)
  y <- rnorm(n)

  fit <- AddiVortes(
    y, x,
    m = 2,
    totalMCMCIter = 4,
    mcmcBurnIn = 0,
    metric = "S",
    showProgress = FALSE
  )

  xnew <- data.frame(lat = lat[1:2], lon = lon[1:2])
  pred <- predict(fit, as.matrix(xnew))

  expect_length(pred, 2)
  expect_true(all(is.finite(pred)))
})
