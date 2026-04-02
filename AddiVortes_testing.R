
devtools::load_all()

sim_fried <- function(N,P,sigma) {
  X <- matrix(runif(N * P), nrow = N, ncol = P)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)**2 + 10 * X[,4] + 5 * X[,5]
  Y <- mu + sigma * rnorm(N)
  return(data.frame(X = X, Y = Y, mu = mu))
}
training_data <- sim_fried(200, 5, 1)
test_data <- sim_fried(400,5, 1)
X_train <- model.matrix(Y ~ . - 1 - mu, data = training_data)
X_test <- model.matrix(Y ~ . - 1 - mu, data = test_data)

cat("Fitting local update algorithm...\n")
time_local <- system.time({
  set.seed(123)
  Model_AddiVortes_local <- AddiVortes(training_data$Y, X_train) #,totalMCMCIter = 200, mcmcBurnIn = 100,distancePower=2, p_rate=0.5, p_sd=0.05)
})

cat("Time taken for local:\n")
print(time_local)

prediction_AddiVortes<-predict(Model_AddiVortes_local,X_test)
rmse_AddiVortes<-sqrt(mean((prediction_AddiVortes-test_data$mu)^2))
print(paste("RMSE for AddiVortes:", rmse_AddiVortes))

cat("Fitting local update algorithm...\n")
time_soft <- system.time({
  set.seed(123)
  Model_softBART<-softbart(X_train,training_data$Y,X_test,opts = Opts(num_burn = 200, num_save = 1000))
})
print(time_soft)

print(paste("RMSE for SoftBART:", sqrt(mean((Model_softBART$y_hat_test_mean - test_data$mu)^2))))

### Graph of times ###

devtools::load_all()

sim_fried <- function(N, P, sigma) {
  X <- matrix(runif(N * P), nrow = N, ncol = P)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)**2 + 10 * X[,4] + 5 * X[,5]
  Y <- mu + sigma * rnorm(N)
  return(data.frame(X = X, Y = Y, mu = mu))
}

N_values <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000)
elapsed_times <- numeric(length(N_values))

for (i in seq_along(N_values)) {
  current_N <- N_values[i]
  cat("Fitting model for N =", current_N, "\n")
  
  training_data <- sim_fried(current_N, 5, 1)
  X_train <- model.matrix(Y ~ . - 1 - mu, data = training_data)
  
  recorded_time <- system.time({
    set.seed(123)
    Model_AddiVortes_local <- AddiVortes(training_data$Y, X_train)
  })
  
  elapsed_times[i] <- recorded_time["elapsed"]
}

plot(N_values, elapsed_times, type = "b", col = "blue", lwd = 2, pch = 16,
     xlab = "Number of Observations", 
     ylab = "Execution Time in Seconds",
     main = "Algorithm Execution Time Over Increasing Sample Sizes")

      