devtools::load_all()

sim_fried <- function(N,P,sigma) {
  X <- matrix(runif(N * P), nrow = N, ncol = P)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)**2 + 10 * X[,4] + 5 * X[,5]
  Y <- mu + sigma * rnorm(N)
  return(data.frame(X = X, Y = Y, mu = mu))
}
training_data <- sim_fried(2000, 5, 1)
test_data <- sim_fried(200,5, 1)
X_train <- model.matrix(Y ~ . - 1 - mu, data = training_data)
X_test <- model.matrix(Y ~ . - 1 - mu, data = test_data)

cat("Fitting local update algorithm...\n")
time_local <- system.time({
  set.seed(123)
  Model_AddiVortes_local <- AddiVortes(training_data$Y, X_train)
})

cat("Time taken for original:\n")
print(time_original)

cat("Time taken for local:\n")
print(time_local)

Model_AddiVortes<-AddiVortes(training_data$Y,X_train)
Model_AddiVortes_local<-AddiVortes_local(training_data$Y,X_train)
prediction_AddiVortes<-predict(Model_AddiVortes_local,X_test)
prediction_AddiVortes_local<-predict(Model_AddiVortes,X_test)

rmse_AddiVortes<-sqrt(mean((prediction_AddiVortes-test_data$Y)^2))
rmse_AddiVortes_local<-sqrt(mean((prediction_AddiVortes_local-test_data$Y)^2))
print(paste("RMSE for AddiVortes_local:", rmse_AddiVortes_local))
print(paste("RMSE for AddiVortes:", rmse_AddiVortes))

      