library(AddiVortes)

# Test 1 -----------------------------------------------------------------------

set.seed(111333)
X <- matrix(rnorm(100),10,10)
Y <- rnorm(10)
X_test <- matrix(rnorm(100),10,10)
Y_test <- rnorm(10)

T1_result <- AddiVortes(Y,X,10,
                        90,10,
                        6,0.85,3,0.8,3,25,
                        Y_test,X_test,
                        IntialSigma = "Linear")

T1_result[[2]]
PredictAddiVortes(T1_result[[1]],
                  X_test,
                  Y_test)[[1]]

#   In_sample_RMSE Out_of_sample_RMSE
# 1      0.7216729           1.051648

# Test 2 -----------------------------------------------------------------------

set.seed(1789)
X <- matrix(runif(500),100, 5)
Y <- rnorm(100, -5, 3)
X_test <- matrix(runif(100),20, 5)
Y_test <- rnorm(20, -5, 3)

T2_result <- AddiVortes(Y,X,5,
                        150,50,
                        6,0.85,3,0.8,3,25,
                        Y_test,X_test,
                        IntialSigma = "Linear")

T2_result[[2]]
PredictAddiVortes(T2_result[[1]],
                  X_test,
                  Y_test)[[1]]

#   In_sample_RMSE Out_of_sample_RMSE
# 1       2.642812           3.010738

# Test 3 -----------------------------------------------------------------------

set.seed(1234)
X <- matrix(rnorm(10000),1000, 10)
Y <- runif(1000, -1, 3)
X_test <- matrix(rnorm(1000),100, 10)
Y_test <- runif(100, -1, 3)

T3_result <- AddiVortes(Y,X,10,
                        200,100,
                        6,0.85,3,0.8,3,25,
                        Y_test,X_test,
                        IntialSigma = "Linear")

T3_result[[2]]
PredictAddiVortes(T3_result[[1]],
                  X_test,
                  Y_test)[[1]]

#   In_sample_RMSE Out_of_sample_RMSE
# 1       1.149055           1.093758

# Test 4 -----------------------------------------------------------------------

set.seed(5678)
X <- matrix(rnorm(10000),1000, 10)
Y <- runif(1000, -1, 3)
X_test <- matrix(rnorm(1000),100, 10)
Y_test <- runif(100, -1, 3)

T4_result <- AddiVortes(Y,X,10,
                        200,100,
                        6,0.85,3,0.8,3,25,
                        Y_test,X_test,
                        IntialSigma = "Linear",
                        thinning = 3)

T4_result[[2]]
PredictAddiVortes(T4_result[[1]],
                  X_test,
                  Y_test)[[1]]

#   In_sample_RMSE Out_of_sample_RMSE
# 1       1.144811           1.080191
