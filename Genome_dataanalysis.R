library(sommer)
library(vcfR)

genotypes <- read.delim("C:/Users/jd191252/Downloads/6281561/File_S3.txt")
phenotype_data <-read.delim("C:/Users/jd191252/Downloads/6281561/File_S1.txt")

common_ids <- intersect(genotypes$Id, phenotype_data$Id)
genotypes <- genotypes[genotypes$Id %in% common_ids, ]
phenotype_data <- phenotype_data[match(genotypes_subset$Id, phenotype_data$Id), ]

# Convert 'Cross' to a numeric column to include it as a feature
cross_effect <- as.numeric(as.factor(phenotype_data$Cross))

# 1. Calculate GRM (A.mat in sommer) from imputed genotypes
# Assuming 'genotypes' is your imputed 0,1,2 matrix
G <- A.mat(as.matrix(genotypes[,-1]))

# CRITICAL: Match the matrix names to the phenotype IDs
# File S3 'Id' column contains the names that G needs
rownames(G) <- colnames(G) <- genotypes$Id

# 2. Prepare Phenotype Data
# Ensure 'Id' matches the names in the matrix and convert 'Cross' to factor
phenotype_data$Id <- as.character(phenotype_data$Id)
phenotype_data$Cross <- as.factor(phenotype_data$Cross)

# 3. Fit the Animal Model
# We use the actual column names from your File_S1 (Weight, Cross, Id)
fit_gblup <- mmer(Weight ~ Cross,
                  random = ~ vsr(Id, Gu = G),
                  rcov = ~ units,
                  data = phenotype_data)

# 4. Extract Predicted Breeding Values (GEBVs)
# Note the capital 'Weight' here to match the formula above
gebv_gblup <- fit_gblup$U$`u:Id`$Weight


#install.packages("BiocManager")

# 2. Use the manager to install the 'impute' library from Bioconductor
#BiocManager::install("impute")
library(impute) # For KNN imputation

# Convert SNP genotypes to a numeric matrix (excluding the ID column)
X_snps <- as.matrix(genotypes[, -1])
storage.mode(X_snps) <- "double"

# Handle remaining NAs with KNN Imputation 
# (AddiVortes will fail if any NAs remain)
if(any(is.na(X_snps))){
  X_snps <- impute.knn(X_snps)$data
}

Y <- phenotype_data$Weight
X <- cbind(X_snps, cross_effect[match(genotypes_subset$Id, phenotype_data$Id)])
data_combined <- data.frame(Y = Y, X)

n <- nrow(data_combined)

train_indices <- sample(1:n, size = round(0.8 * n))
train_set <- data_combined[train_indices, ]
test_set  <- data_combined[-train_indices, ]
X_train <- model.matrix(Y ~ . - 1, data = train_set)
X_test  <- model.matrix(Y ~ . - 1, data = test_set)
Y_train <- train_set$Y
Y_test  <- test_set$Y
non_constant_cols <- apply(X_train, 2, function(col) length(unique(col)) > 1)
X_train <- X_train[, non_constant_cols, drop = FALSE]
X_test  <- X_test[, non_constant_cols, drop = FALSE]
training_data <- data.frame(Y = Y_train, X_train)

# Fit AddiVortes
fit_vortes <- AddiVortes( Y_train,X_train,m = 20,
                         totalMCMCIter = 5000, mcmcBurnIn = 2500, dirichletWarmup = 1000,
                         numChains = 1, IntialSigma = "LASSO", varSelMode = 2,splitMode = 0)

# Predict
predictions_vortes <- predict(fit_vortes, newdata = X_test)

cat("Fitting local update algorithm...\n")
time_soft <- system.time({
  set.seed(123)
  Model_softBART<-softbart(X_train,Y_train,X_test,hypers = Hypers(X_train, Y_train,num_tree = 20))#,opts = Opts(num_burn = 200, num_save = 1000))
})
print(paste("RMSE for SoftBART:", sqrt(mean((Model_softBART$y_hat_test_mean - Y_test)^2))))

# 1. Define RMSE function
calc_rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

# 2. Calculate metrics for G-BLUP (using your previous model results)
# Assuming gebv_gblup contains the predictions from mmer
rmse_gblup <- calc_rmse(phenotype_data$Weight, gebv_gblup)
cor_gblup  <- cor(phenotype_data$Weight, gebv_gblup)

# 3. Calculate metrics for AddiVortes
rmse_vortes <- calc_rmse(Y_test, predictions_vortes)
cor_vortes  <- cor(Y_test, predictions_vortes)

rmse_baseline<-calc_rmse(Y_test, mean(Y_test))
cor_baseline<-NA

# 4. Display Comparison Table
results <- data.frame(
  Metric = c("RMSE", "Correlation (r)"),
  G_BLUP = c(rmse_gblup, cor_gblup),
  AddiVortes = c(rmse_vortes, cor_vortes),
  baseline = c(rmse_baseline,cor_baseline)
)

print(results)

