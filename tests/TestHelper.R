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
AddiVortes(Y_Boston[TrainSet],X_Boston[TrainSet,],
           200,2000,200,6,0.85,3,0.8,3,25,
           Y_Boston[TestSet],X_Boston[TestSet,],
           IntialSigma = "Linear")
toc()

# Iteration 2000 out of 2000
# In_sample_RMSE Out_of_sample_RMSE
# 1       1.206556            3.20922
# In_sample_RMSE Out_of_sample_RMSE (vectorised sampling)
# 1       1.204214           3.237258

# Timings (Office comp)
# v0.0.14 = 87 sec elapsed
# v0.0.15 = 84 sec elapsed
# v0.0.16 = 82 sec elapsed
# v0.0.17 = 82 sec elapsed
# v0.0.20 = 115 sec elapsed
# v0.0.21 = 118 sec elapsed
