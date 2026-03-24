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

      