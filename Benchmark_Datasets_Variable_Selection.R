### Triazines ###
{# Import the dataset into a data frame
triazines_data <- read.csv("https://raw.githubusercontent.com/s3628730/datasets/master/triazines.csv")

# View the first few rows of the data
head(triazines_data)
X_triazines <- triazines_data[,1:60]
X_triazines <- X_triazines[,-c(43,44)]
Y_triazines <- triazines_data[,61]
}

### Boston Housing ###
{
  boston_data <- read.csv("https://raw.githubusercontent.com/s3628730/datasets/master/boston_housing.csv")
  
  # View the first few rows of the data
  head(boston_data)
  X_boston <- boston_data[,1:13]
  Y_boston <- boston_data[,14]
}


### Gisette ###
{
  gisette_data <- read.csv("C:/Users/adam6/Downloads/gisette.csv/gisette.csv")
  
  # View the first few rows of the data
  head(gisette_data)
  X_gisette <- gisette_data[,1:5000]
  Y_gisette <- gisette_data[,5001]
}


install.packages("hdi")

set.seed(123)

data <- data.frame(Y = Y_gisette, X_gisette)
# Determine the number of rows
n <- nrow(data)
train_indices <- sample(1:n, size = round(0.8 * n))
train_set <- data[train_indices, ]
test_set  <- data[-train_indices, ]
X_train <- model.matrix(Y ~ . -1, data = train_set)
X_test  <- model.matrix(Y ~ . -1, data = test_set)
Y_train <- train_set$Y
Y_test  <- test_set$Y
training_data<-data.frame(Y=Y_train,X_train)

{# Identify columns in X_train that have more than one unique value
non_constant_cols <- apply(X_train, 2, function(col) length(unique(col)) > 1)

# Subset both X_train and X_test using the same logical vector
X_train <- X_train[, non_constant_cols]
X_test  <- X_test[, non_constant_cols]

# Update the training_data data frame with the filtered predictors
training_data <- data.frame(Y = Y_train, X_train)}


cat("Fitting local update algorithm...\n")
time_local <- system.time({
  Model_AddiVortes_local <- AddiVortes(Y_train, X_train,m=20,totalMCMCIter = 2000, mcmcBurnIn = 1000,dirichletWarmup = 500,a_alpha = 0.5, b_alpha = 1,adaptBoost = 1, adaptPenalty = 1, momentumDecay = 0.6, kappa = 0.6,numChains = 1,IntialSigma = "Naive")
})

cat("Time taken for local:\n")
print(time_local)
##Graphs
{
  # Assuming you have run the algorithm and saved the output to 'results':
  # results <- Homo_AddiVortes_Algorithm(...)
  
  # 1. Extract the two metrics from the results list
  inclusion_probs <- Model_AddiVortes_local$ensembleInclusionProbabilities
  
  s_means <- Model_AddiVortes_local$posteriorDirichletWeightsMean
  s_lower <- Model_AddiVortes_local$posteriorDirichletWeightsLower
  s_upper <- Model_AddiVortes_local$posteriorDirichletWeightsUpper
  covariate_index <- 1:length(s_means)
  
  plot(covariate_index, s_means, type = "p", pch = 19, col = "darkorange",
       ylim = c(0, max(s_upper)+0.05),
       main = "Posterior Dirichlet Weights with 95% Credible Intervals",
       xlab = "Covariate Index", ylab = "Weight")
  
  arrows(x0 = covariate_index, y0 = s_lower, x1 = covariate_index, y1 = s_upper,
         angle = 90, code = 3, length = 0.05, col = "darkgray")
  # 4. Plot the Ensemble Inclusion Probabilities
  # We set ylim = c(0,1) because this is a strict probability bound
  inclusion_probs <- Model_AddiVortes_local$ensembleInclusionProbabilities
  
  barplot(inclusion_probs, 
          main = "Ensemble Inclusion Probabilities", 
          xlab = "Covariate Index", 
          ylab = "Probability (0 to 1)",
          col = "steelblue",
          border = NA,
          ylim = c(0, 1),
          names.arg = 1:length(inclusion_probs))
  
  #plot(Model_AddiVortes_local$posteriorAlphaTrace)
}
prediction_AddiVortes<-predict(Model_AddiVortes_local,X_test)
rmse_AddiVortes<-sqrt(mean((prediction_AddiVortes-Y_test)^2))
print(paste("RMSE for AddiVortes:", rmse_AddiVortes))
print(paste("Benchmark:", sqrt(mean((mean(Y_test)-Y_test)^2))))

cat("Fitting Orginal update algorithm...\n")
time_local <- system.time({
  Model_AddiVortes_Original <- AddiVortes(Y_train, X_train,varSelMode = 0,m=200,totalMCMCIter = 5000, mcmcBurnIn = 2000,dirichletWarmup = 1000,a_alpha = 0.5, b_alpha = 1,adaptBoost = 1, adaptPenalty = 1, momentumDecay = 0.9, kappa = 0.6,numChains = 1,IntialSigma = "LASSO")
})

cat("Time taken for local:\n")
print(time_local)
##Graphs
{
  # Assuming you have run the algorithm and saved the output to 'results':
  # results <- Homo_AddiVortes_Algorithm(...)
  
  # 1. Extract the two metrics from the results list
  inclusion_probs <- Model_AddiVortes_Original$ensembleInclusionProbabilities
  
  s_means <- Model_AddiVortes_Original$posteriorDirichletWeightsMean
  s_lower <- Model_AddiVortes_Original$posteriorDirichletWeightsLower
  s_upper <- Model_AddiVortes_Original$posteriorDirichletWeightsUpper
  covariate_index <- 1:length(s_means)
  
  plot(covariate_index, s_means, type = "p", pch = 19, col = "darkorange",
       ylim = c(0, max(s_upper)+0.05),
       main = "Posterior Dirichlet Weights with 95% Credible Intervals",
       xlab = "Covariate Index", ylab = "Weight")
  
  arrows(x0 = covariate_index, y0 = s_lower, x1 = covariate_index, y1 = s_upper,
         angle = 90, code = 3, length = 0.05, col = "darkgray")
  # 4. Plot the Ensemble Inclusion Probabilities
  # We set ylim = c(0,1) because this is a strict probability bound
  inclusion_probs <- Model_AddiVortes_Original$ensembleInclusionProbabilities
  
  barplot(inclusion_probs, 
          main = "Ensemble Inclusion Probabilities", 
          xlab = "Covariate Index", 
          ylab = "Probability (0 to 1)",
          col = "steelblue",
          border = NA,
          ylim = c(0, 1),
          names.arg = 1:length(inclusion_probs))
  
  #plot(Model_AddiVortes_Original$posteriorAlphaTrace)
}
prediction_AddiVortes_Original<-predict(Model_AddiVortes_Original,X_test)
rmse_AddiVortes_Original<-sqrt(mean((prediction_AddiVortes_Original-Y_test)^2))
print(paste("RMSE for AddiVortes:", rmse_AddiVortes_Original))

cat("Fitting Orginal update algorithm...\n")
time_local <- system.time({
  Model_AddiVortes_Dir <- AddiVortes(Y_train, X_train,varSelMode = 1,m=200,totalMCMCIter = 5000, mcmcBurnIn = 2000,dirichletWarmup = 1000,a_alpha = 0.5, b_alpha = 1,adaptBoost = 1, adaptPenalty = 1, momentumDecay = 0.9, kappa = 0.6,numChains = 1,IntialSigma = "LASSO")
})

cat("Time taken for local:\n")
print(time_local)
##Graphs
{
  # Assuming you have run the algorithm and saved the output to 'results':
  # results <- Homo_AddiVortes_Algorithm(...)
  
  # 1. Extract the two metrics from the results list
  inclusion_probs <- Model_AddiVortes_Dir$ensembleInclusionProbabilities
  
  s_means <- Model_AddiVortes_Dir$posteriorDirichletWeightsMean
  s_lower <- Model_AddiVortes_Dir$posteriorDirichletWeightsLower
  s_upper <- Model_AddiVortes_Dir$posteriorDirichletWeightsUpper
  covariate_index <- 1:length(s_means)
  
  plot(covariate_index, s_means, type = "p", pch = 19, col = "darkorange",
       ylim = c(0, max(s_upper)+0.05),
       main = "Posterior Dirichlet Weights with 95% Credible Intervals",
       xlab = "Covariate Index", ylab = "Weight")
  
  arrows(x0 = covariate_index, y0 = s_lower, x1 = covariate_index, y1 = s_upper,
         angle = 90, code = 3, length = 0.05, col = "darkgray")
  # 4. Plot the Ensemble Inclusion Probabilities
  # We set ylim = c(0,1) because this is a strict probability bound
  inclusion_probs <- Model_AddiVortes_Dir$ensembleInclusionProbabilities
  
  barplot(inclusion_probs, 
          main = "Ensemble Inclusion Probabilities", 
          xlab = "Covariate Index", 
          ylab = "Probability (0 to 1)",
          col = "steelblue",
          border = NA,
          ylim = c(0, 1),
          names.arg = 1:length(inclusion_probs))
  
  #plot(Model_AddiVortes_Dir$posteriorAlphaTrace)
}
prediction_AddiVortes_Dir<-predict(Model_AddiVortes_Dir,X_test)
rmse_AddiVortes_Dir<-sqrt(mean((prediction_AddiVortes_Dir-Y_test)^2))
print(paste("RMSE for AddiVortes:", rmse_AddiVortes_Dir))

task<- system.time({
  set.seed(123)
  Model_flexBART <- flexBART(Y~bart(.), train_data = training_data,M_vec=20,save_trees = TRUE,n.chains=1)
  prediction_flexBART<-predict.flexBART(Model_flexBART,test_set)
})
print(task)
rmse_flexBART<-sqrt(mean((colMeans(prediction_flexBART)-Y_test)^2))
#rmse_flexBART<-sqrt(mean((Model_flexBART$yhat.test.mean-test_data$mu)^2))
print(paste("RMSE for flexBART:", rmse_flexBART))

