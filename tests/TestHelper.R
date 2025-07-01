library(AddiVortes)
library(tictoc)

Boston <- read.csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/BostonHousing_Data.csv")
X_Boston <- as.matrix(Boston[,2:14])
Y_Boston <- as.numeric(as.matrix(Boston[,15]))
rm(Boston)

n <- length(Y_Boston)

set.seed(1025)
TrainSet <- sort(sample.int(n,5*n/6))
TestSet <- 1:n
TestSet <- TestSet[! TestSet %in% TrainSet]

tic("Start")
results <- AddiVortes(Y_Boston[TrainSet],X_Boston[TrainSet,],
           200,2000,200,6,0.85,3,0.8,3,25,
           Y_Boston[TestSet],X_Boston[TestSet,],
           IntialSigma = "Linear")
preds <- predictAddiVortes(results[[1]],
                           X_Boston[TestSet,],
                           Y_Boston[TestSet])
toc()

results[[2]][[1]]
preds[[1]]

# Plot predictions
plot(Y_Boston[TestSet],
     preds[[2]],
     xlab = "True Values",
     ylab = "Predicted Values",
     main = "AddiVortes Predictions vs True Values",
     xlim = c(min(Y_Boston[TestSet]) - 0.1, max(Y_Boston[TestSet]) + 0.1),
     ylim = c(min(preds[[3]]) - 0.1, max(preds[[3]]) + 0.1),
     pch = 19, col = "red")
# Add error lines
for (i in 1:ncol(preds[[3]])){
  segments(Y_Boston[TestSet][i], preds[[3]][1,i],
           Y_Boston[TestSet][i], preds[[3]][2,i],
           col = "red", lwd = 1.5)
}
# Add in the equality line
abline(a = 0, b = 1, col = "blue", lwd = 2)

# Iteration 2000 out of 2000
# In_sample_RMSE Out_of_sample_RMSE
# 1       1.206556           3.209220
# In_sample_RMSE Out_of_sample_RMSE (vectorised sampling)
# 1       1.204092           3.237162

# Timings (Office comp)
# v0.0.14 = 87 sec elapsed
# v0.0.15 = 84 sec elapsed
# v0.0.16 = 82 sec elapsed
# v0.0.17 = 82 sec elapsed
# v0.0.20 = 115 sec elapsed
# v0.0.21 = 117 sec elapsed
# v0.0.29 = 123 sec elapsed
