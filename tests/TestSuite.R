require(AddiVortes)

# Test 1 -----------------------------------------------------------------------

set.seed(111333)
X <- matrix(rnorm(100), 10, 10)
Y <- rnorm(10)
X_test <- matrix(rnorm(100), 10, 10)
Y_test <- rnorm(10)

T1_result <- AddiVortes(Y, X, 10,
  90, 10,
  6, 0.85, 3, 0.8, 3, 25,
  IntialSigma = "Linear"
)

T1_result[[8]]
Predicts <- predict(
  T1_result,
  X_test
)
PredictsQ <- predict(
  T1_result,
  X_test,
  "quantile"
)
sqrt(mean((Y_test - Predicts)^2))
plot(T1_result, X, Y)

#   In_sample_RMSE Out_of_sample_RMSE
# 1      0.7216729           1.051648

# Test 2 -----------------------------------------------------------------------

set.seed(1789)
X <- matrix(runif(500), 100, 5)
Y <- rnorm(100, -5, 3)
X_test <- matrix(runif(100), 20, 5)
Y_test <- rnorm(20, -5, 3)

T2_result <- AddiVortes(Y, X, 5,
  150, 50,
  6, 0.85, 3, 0.8, 3, 25,
  IntialSigma = "Linear"
)

T2_result[[8]]
Predicts <- predict(
  T2_result,
  X_test
)
sqrt(mean((Y_test - Predicts)^2))

#   In_sample_RMSE Out_of_sample_RMSE
# 1       2.642812           3.010738

# Test 3 -----------------------------------------------------------------------

set.seed(1234)
X <- matrix(rnorm(10000), 1000, 10)
Y <- runif(1000, -1, 3)
X_test <- matrix(rnorm(1000), 100, 10)
Y_test <- runif(100, -1, 3)

T3_result <- AddiVortes(Y, X, 10,
  200, 100,
  6, 0.85, 3, 0.8, 3, 25,
  IntialSigma = "Linear"
)

T3_result[[8]]
Predicts <- predict(
  T3_result,
  X_test
)
sqrt(mean((Y_test - Predicts)^2))

#   In_sample_RMSE Out_of_sample_RMSE
# 1       1.149055           1.093758

# Test 4 -----------------------------------------------------------------------
library(tictoc)
set.seed(5678)

X <- matrix(rnorm(10000), 1000, 10)
Y <- runif(1000, -1, 3)
X_test <- matrix(rnorm(1000), 100, 10)
Y_test <- runif(100, -1, 3)

tic("AddiVortes Test 4")
T4_result <- AddiVortes(Y, X, 10,
  2000, 100,
  6, 0.85, 3, 0.8, 3, 25,
  IntialSigma = "Linear",
  thinning = 3
)
toc()

T4_result[[8]]
Predicts <- predict(
  T4_result,
  X_test
)
sqrt(mean((Y_test - Predicts)^2))

#   In_sample_RMSE Out_of_sample_RMSE
# 1       1.140032           1.069564
