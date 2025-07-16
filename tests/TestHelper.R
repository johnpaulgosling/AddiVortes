require(AddiVortes)
require(tictoc)

Boston <- read.csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/BostonHousing_Data.csv")
X_Boston <- as.matrix(Boston[, 2:14])
Y_Boston <- as.numeric(as.matrix(Boston[, 15]))
rm(Boston)

n <- length(Y_Boston)

set.seed(1025)
TrainSet <- sort(sample.int(n, 5 * n / 6))
TestSet <- 1:n
TestSet <- TestSet[!TestSet %in% TrainSet]

tic("Start")
results <- AddiVortes(Y_Boston[TrainSet], X_Boston[TrainSet, ],
  200, 2000, 200, 6, 0.85, 3, 0.8, 3, 25,
  IntialSigma = "Linear"
)
preds <- predict(
  results,
  X_Boston[TestSet, ]
)
toc()

results[8]
sqrt(mean((Y_Boston[TestSet] - preds)^2))

# Plot predictions
plot(Y_Boston[TestSet],
  preds,
  xlab = "True Values",
  ylab = "Predicted Values",
  main = "AddiVortes Predictions vs True Values",
  xlim = c(min(Y_Boston[TestSet]) - 0.1, max(Y_Boston[TestSet]) + 0.1),
  ylim = c(min(preds) - 0.1, max(preds) + 0.1),
  pch = 19, col = "red"
)
# Add error lines
preds <- predict(
  results,
  X_Boston[TestSet, ],
  "quantile"
)
for (i in 1:nrow(preds)) {
  segments(Y_Boston[TestSet][i], preds[i, 1],
    Y_Boston[TestSet][i], preds[i, 2],
    col = "red", lwd = 1.5
  )
}
# Add in the equality line
abline(a = 0, b = 1, col = "blue", lwd = 2)

# Iteration 2000 out of 2000
# In_sample_RMSE Out_of_sample_RMSE
# 1       1.206556           3.209220
# In_sample_RMSE Out_of_sample_RMSE (vectorised sampling)
# 1       1.204092           3.237162
# In_sample_RMSE Out_of_sample_RMSE (computation speed-up)
# 1       1.172312           3.116504
# In_sample_RMSE Out_of_sample_RMSE (C++)
# 1       1.185435           3.060899

# Timings (Office comp)
# v0.0.14 = 87 sec elapsed
# v0.0.15 = 84 sec elapsed
# v0.0.16 = 82 sec elapsed
# v0.0.17 = 82 sec elapsed
# v0.0.20 = 115 sec elapsed
# v0.0.21 = 117 sec elapsed
# v0.0.29 = 123 sec elapsed
