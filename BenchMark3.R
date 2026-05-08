## Mice Genome dataset ###
{library(BGLR)
data(mice, package = "BGLR")

Y_mice <- as.numeric(mice.pheno$Obesity.BMI)
X_mice <- as.matrix(mice.X)

valid_indices <- !is.na(Y_mice)
Y_mice <- Y_mice[valid_indices]
X_mice <- X_mice[valid_indices, ]

data <- data.frame(Y = Y_mice, X_mice)

set.seed(123)
n <- nrow(data)
train_indices <- sample(1:n, size = round(0.8 * n))

train_set <- data[train_indices, ]
test_set  <- data[-train_indices, ]

X_train <- model.matrix(Y ~ . -1, data = train_set)
X_test  <- model.matrix(Y ~ . -1, data = test_set)
Y_train <- train_set$Y
Y_test  <- test_set$Y

training_data <- data.frame(Y = Y_train, X_train)
}

## Wheat Genome dataset ##
{library(BGLR)
  data(wheat, package = "BGLR")
  
  Y_maize <- as.numeric(wheat.Y[, 1])
  X_maize <- as.matrix(wheat.X)
  
  data <- data.frame(Y = Y_maize, X_maize)
  
  set.seed(123)
  n <- nrow(data)
  train_indices <- sample(1:n, size = round(0.8 * n))
  
  train_set <- data[train_indices, ]
  test_set  <- data[-train_indices, ]
  
  X_train <- model.matrix(Y ~ . -1, data = train_set)
  X_test  <- model.matrix(Y ~ . -1, data = test_set)
  Y_train <- train_set$Y
  Y_test  <- test_set$Y
  
  training_data <- data.frame(Y = Y_train, X_train)
  }

### cancer dataset ###
{
# Install the package if you have not already
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("TCGAbiolinks")

#library(TCGAbiolinks)
#library(SummarizedExperiment)

query_initial <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

results <- getResults(query_initial)

query_subset <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = results$cases[1:50]
)

GDCdownload(query_subset)
tcga_data <- GDCprepare(query_subset)

clinical_data <- colData(tcga_data)
Y_cancer <- as.numeric(clinical_data$age_at_diagnosis) / 365.25

X_cancer <- assay(tcga_data)
X_cancer <- t(X_cancer[1:1000, ])

valid_indices <- !is.na(Y_cancer)
Y_cancer <- Y_cancer[valid_indices]
X_cancer <- X_cancer[valid_indices, ]

set.seed(123)
n <- length(Y_cancer)
train_indices <- sample(1:n, size = round(0.8 * n))

X_train <- X_cancer[train_indices, ]
Y_train <- Y_cancer[train_indices]
X_test  <- X_cancer[-train_indices, ]
Y_test  <- Y_cancer[-train_indices]
}

cat("Fitting local update algorithm...\n")
time_local <- system.time({
  Model_AddiVortes_local <- AddiVortes(Y_train, X_train,m=20,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,a_alpha = 0.5, b_alpha = 1,adaptBoost = 1, adaptPenalty = 1, momentumDecay = 0.6, kappa = 0.6,numChains = 2,IntialSigma = "LASSO")
})

cat("Time taken for local:\n")
print(time_local)

##Graphs##
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
print(sqrt(mean((mean(Y_test)-Y_test)^2)))

task<- system.time({
  set.seed(123)
  Model_flexBART <- flexBART(Y~bart(.), train_data = training_data,M_vec=20,save_trees = TRUE,n.chains=2,sparse = TRUE )
  prediction_flexBART<-predict.flexBART(Model_flexBART,test_set)
})
print(task)
rmse_flexBART<-sqrt(mean((colMeans(prediction_flexBART)-Y_test)^2))
#rmse_flexBART<-sqrt(mean((Model_flexBART$yhat.test.mean-test_data$mu)^2))
print(paste("RMSE for flexBART:", rmse_flexBART))

