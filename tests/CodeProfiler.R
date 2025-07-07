require(AddiVortes)

# Code profiling
require(profvis)

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
profvis::profvis({
  # Run the algorithm
  result <- AddiVortes(Y_Boston[TrainSet],X_Boston[TrainSet,],
                       200,2000,200,6,0.85,3,0.8,3,25,
                       Y_Boston[TestSet],X_Boston[TestSet,],
                       IntialSigma = "Linear")
})
