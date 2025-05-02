library(AddiVortes)

# Code profiling
library(profvis)

Boston <- read.csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/BostonHousing_Data.csv")
X_Boston <- as.matrix(Boston[,2:14])
Y_Boston <- as.numeric(as.matrix(Boston[,15]))
rm(Boston)

n <- length(Y_Boston)

set.seed(1025)
TrainSet <- sort(sample.int(n,5*n/6))
TestSet <- 1:n
TestSet <- TestSet[! TestSet %in% TrainSet]

# Code profiling
profvis({
  # Run the algorithm
  result <- AddiVortes(X_Boston[TrainSet,], Y_Boston[TrainSet], X_Boston[TestSet,], Y_Boston[TestSet],
                       nfolds = 5, ntrees = 100, maxdepth = 10, minsize = 5,
                       alpha = 0.05, beta = 0.05, gamma = 0.05,
                       verbose = TRUE)
})
