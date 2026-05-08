results_summary_full <- data.frame(
  Dataset = character(),
  Metric_Type = character(),
  AddiVortes_Mode0_Score = numeric(),
  AddiVortes_Mode1_Score = numeric(),
  AddiVortes_Mode2_Score = numeric(),
  flexBART_Score = numeric(),
  softBART_Score = numeric(),
  wBART_Score = numeric(),
  RF_Score = numeric(),
  #SVM_Score = numeric(),
  XGB_Score = numeric(),
  baseline_score = numeric(),
  stringsAsFactors = FALSE
)

# ======================================================================
# GLOBAL SETTINGS & HYPERPARAMETERS
# ======================================================================

n_runs <- 20
num_cores <- 20

# Centralised MCMC Parameters
MCMC_ITER <- 5000
MCMC_BURNIN <- 2500
DIRICHLET_WARMUP <- 1000

# ======================================================================

library(doParallel)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  # If your script is running from the package root:
  devtools::load_all() 
  # Ensure all necessary packages are available to workers
  library(SoftBart)
  library(BART)
  library(ranger)
  library(e1071)
  library(xgboost)
})

### CV Hyperparameters ###
{
  # Base hyperparameters for Mode 0
  base_grid_0 <- expand.grid(
    m = c(20, 200),
    k = c(1, 3),
    Omega = c(1, 3),
    LambdaRate = c(5, 25)
  )
  
  # Specific (nu, q) pairs
  nu_q_cases <- data.frame(
    nu = c(3, 6),
    q = c(0.99, 0.85)
  )
  
  # Cross-join to create the final Mode 0 grid
  grid_mode0 <- merge(base_grid_0, nu_q_cases, by = NULL)
  
  hyperparam_cases_1 <- data.frame(
    alpha = c(1.0),
    a_alpha = c(0.5,1),
    b_alpha = c(1.0)
  )
  grid_mode1 <- merge(grid_mode0, hyperparam_cases_1, by = NULL)
  
  hyperparam_cases_2 <- data.frame(
    Case = c("Strict", "Balanced", "Permissive"),
    alpha = c(1,1,1),
    a_alpha = c(0.1, 0.5, 1),
    b_alpha = c(2.0, 1.0, 0.5),
    adaptBoost = c(1.0, 2.0, 5.0),
    adaptPenalty = c(8.0, 2.0, 0.1),
    momentumDecay = c(0.50, 0.90, 0.99),
    stringsAsFactors = FALSE
  )
  grid_mode2 <- merge(grid_mode0, hyperparam_cases_2, by = NULL)
  
  # --- EXPANDED SMART flexBART GRID ---
  # Base hyperparameters
  base_grid_flex <- expand.grid(
    M_vec = c(20, 200)
  )
  
  # (nu, q) pairs for the error variance prior
  # (3, 0.90) is standard BART, (3, 0.99) forces tighter fits, (10, 0.75) is highly constrained
  nu_q_cases <- data.frame(
    nu = c(3, 3, 10),
    q = c(0.90, 0.99, 0.75)
  )
  
  # (alpha, beta) pairs for the tree depth prior
  # (0.95, 2.0) is standard BART, (0.99, 3.0) allows for much deeper, complex trees
  alpha_beta_cases <- data.frame(
    alpha = c(0.95, 0.8),
    beta = c(2.0, 3.0)
  )
  
  sparsity_cases_flex<- data.frame(
    a_u = c( 0.1, 0.5, 1),
    b_u = c(1)
  )
  
  # Cross-join everything together
  grid_flexBART <- merge(base_grid_flex, nu_q_cases, by = NULL)
  grid_flexBART <- merge(grid_flexBART, alpha_beta_cases, by = NULL)
  grid_flexBART <- merge(grid_flexBART, sparsity_cases_flex, by = NULL)
  
  
  grid_softBART <- expand.grid(
    num_tree = c(20, 50, 100),
    alpha = c(0.95, 0.99), # 0.95 is standard BART, 0.99 allows deeper trees
    beta = c(2.0, 3.0)     # 2.0 is standard BART
  )
  grid_wBART <- expand.grid(
    ntree = c(20, 50, 100),
    base = c(0.95, 0.99),
    power = c(2.0, 3.0)
  )
  grid_rf <- expand.grid(
    mtry_prop = c(0.2, 0.4, 0.6, 0.8),
    min.node.size = c(1, 3, 5, 10),
    sample.fraction = c(0.632, 0.8),
    num.trees = c(250, 500),
    replace = c(TRUE, FALSE),
    use_extratrees = c(TRUE, FALSE)
  )
  # --- EXPANDED SMART SVM GRID ---
  # Radial Basis Function configurations
  grid_svm_radial <- expand.grid(
    kernel = "radial",
    cost = c(0.1, 1, 10, 50),
    gamma = c(0.001, 0.01, 0.1, 0.5),
    degree = 3,  # Placeholder (ignored by radial)
    coef0 = 0,   # Placeholder (ignored by radial)
    epsilon = c(0.01, 0.1)
  )
  
  # Polynomial configurations
  grid_svm_poly <- expand.grid(
    kernel = "polynomial",
    cost = c(0.1, 1, 10, 50),
    gamma = c(0.01, 0.1),
    degree = c(2, 3),
    coef0 = c(0, 1),
    epsilon = c(0.01, 0.1)
  )
  
  # Combine into one master grid
  grid_svm <- rbind(grid_svm_radial, grid_svm_poly)
  
  
  grid_xgb <- expand.grid(
    nrounds = c(100, 300), 
    max_depth = c(3, 5, 7), 
    learning_rate = c(0.01, 0.05, 0.1),
    subsample = c(0.7, 1.0),
    colsample_bytree = c(0.7, 1.0)
  )
}

createFolds <- function(y, k = 5, returnTrain = TRUE) {
  n <- length(y)
  shuffled_indices <- sample(1:n)
  folds <- split(shuffled_indices, cut(seq_along(shuffled_indices), breaks = k, labels = FALSE))
  if (returnTrain) {
    train_folds <- lapply(folds, function(test_idx) setdiff(1:n, test_idx))
    names(train_folds) <- paste0("Fold", 1:k)
    return(train_folds)
  } else {
    names(folds) <- paste0("Fold", 1:k)
    return(folds)
  }
}

# --- Helper to drop constant columns dynamically per CV fold ---
get_cv_data <- function(X_train, Y_train, train_idx) {
  X_cv_train <- X_train[train_idx, , drop = FALSE]
  X_cv_val   <- X_train[-train_idx, , drop = FALSE]
  Y_cv_train <- Y_train[train_idx]
  Y_cv_val   <- Y_train[-train_idx]
  
  non_const <- apply(X_cv_train, 2, function(col) length(unique(col)) > 1)
  
  list(
    X_train = X_cv_train[, non_const, drop = FALSE],
    X_val   = X_cv_val[, non_const, drop = FALSE],
    Y_train = Y_cv_train,
    Y_val   = Y_cv_val
  )
}

for (data_name in names(benchmark_datasets)) {
  cat("\nAnalysing dataset:", data_name, "\n")
  current_data <- benchmark_datasets[[data_name]]
  
  Y <- current_data$Y
  X <- current_data$X
  
  is_classification <- length(unique(na.omit(Y))) == 2
  metric_name <- ifelse(is_classification, "Accuracy", "RMSE")
  cat("Detected task:", ifelse(is_classification, "Binary Classification", "Regression"), "\n")
  cat("Target Metric:", metric_name, "\n\n")
  
  if (is_classification) {
    unique_vals <- sort(unique(na.omit(Y)))
    Y <- ifelse(Y == unique_vals[2], 1, 0)
  }
  
  data_combined <- data.frame(Y = Y, X)
  n <- nrow(data_combined)
  set.seed(123)
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
  cv_folds <- createFolds(Y_train, k = 5, returnTrain = TRUE)
  
  packages_list <- c("glmnet", "flexBART", "SoftBart", "BART", "ranger", "e1071", "xgboost")
  
  # --- TUNING AddiVortes Mode 0 ---
  cat("Tuning AddiVortes Mode 0...\n")
  cv_results_0 <- foreach(i = 1:nrow(grid_mode0), .combine = rbind, .packages = packages_list) %dopar% {
    params <- grid_mode0[i, , drop = FALSE]
    fold_score <- numeric(5)
    for (f in 1:5) {
      cv <- get_cv_data(X_train, Y_train, cv_folds[[f]])
      
      model <- AddiVortes(cv$Y_train, cv$X_train, 
                          m = params$m, 
                          k = params$k,
                          Omega = params$Omega,
                          LambdaRate = params$LambdaRate,
                          nu = params$nu,
                          q = params$q,
                          totalMCMCIter = MCMC_ITER, mcmcBurnIn = MCMC_BURNIN, 
                          numChains = 1, showProgress = FALSE, 
                          IntialSigma = "LASSO", varSelMode = 0)
      preds <- predict(model, cv$X_val, showProgress = FALSE)
      fold_score[f] <- if (is_classification) mean(round(preds) == cv$Y_val) else sqrt(mean((preds - cv$Y_val)^2))
    }
    data.frame(params, avg_score = mean(fold_score))
  }
  best_params_0 <- if (is_classification) cv_results_0[which.max(cv_results_0$avg_score), ] else cv_results_0[which.min(cv_results_0$avg_score), ]
  
  # --- TUNING AddiVortes Mode 1 ---
  cat("Tuning AddiVortes Mode 1...\n")
  cv_results_1 <- foreach(i = 1:nrow(grid_mode1), .combine = rbind, .packages = packages_list) %dopar% {
    params <- grid_mode1[i, ]
    fold_score <- numeric(5)
    for (f in 1:5) {
      cv <- get_cv_data(X_train, Y_train, cv_folds[[f]])
      
      model <- AddiVortes(cv$Y_train, cv$X_train, 
                          m = params$m, 
                          k = params$k,
                          Omega = params$Omega,
                          LambdaRate = params$LambdaRate,
                          nu = params$nu,
                          q = params$q,
                          a_alpha = params$a_alpha, b_alpha = params$b_alpha,
                          totalMCMCIter = MCMC_ITER, mcmcBurnIn = MCMC_BURNIN, dirichletWarmup = DIRICHLET_WARMUP,
                          numChains = 1, showProgress = FALSE, IntialSigma = "LASSO", varSelMode = 1)
      preds <- predict(model, cv$X_val, showProgress = FALSE)
      fold_score[f] <- if (is_classification) mean(round(preds) == cv$Y_val) else sqrt(mean((preds - cv$Y_val)^2))
    }
    data.frame(params, avg_score = mean(fold_score))
  }
  best_params_1 <- if (is_classification) cv_results_1[which.max(cv_results_1$avg_score), ] else cv_results_1[which.min(cv_results_1$avg_score), ]
  
  # --- TUNING AddiVortes Mode 2 ---
  cat("Tuning AddiVortes Mode 2...\n")
  cv_results_2 <- foreach(i = 1:nrow(grid_mode2), .combine = rbind, .packages = packages_list) %dopar% {
    params <- grid_mode2[i, ]
    fold_score <- numeric(5)
    for (f in 1:5) {
      cv <- get_cv_data(X_train, Y_train, cv_folds[[f]])
      
      model <- AddiVortes(cv$Y_train, cv$X_train, 
                          m = params$m, 
                          k = params$k,
                          Omega = params$Omega,
                          LambdaRate = params$LambdaRate,
                          nu = params$nu,
                          q = params$q,
                          a_alpha = params$a_alpha, b_alpha = params$b_alpha,
                          adaptBoost = params$adaptBoost, adaptPenalty = params$adaptPenalty, momentumDecay = params$momentumDecay,
                          totalMCMCIter = MCMC_ITER, mcmcBurnIn = MCMC_BURNIN, dirichletWarmup = DIRICHLET_WARMUP,
                          numChains = 1, showProgress = FALSE, IntialSigma = "LASSO", varSelMode = 2)
      preds <- predict(model, cv$X_val, showProgress = FALSE)
      fold_score[f] <- if (is_classification) mean(round(preds) == cv$Y_val) else sqrt(mean((preds - cv$Y_val)^2))
    }
    data.frame(params, avg_score = mean(fold_score))
  }
  best_params_2 <- if (is_classification) cv_results_2[which.max(cv_results_2$avg_score), ] else cv_results_2[which.min(cv_results_2$avg_score), ]
  
  # --- TUNING flexBART ---
  cat("Tuning flexBART...\n")
  cv_results_flex <- foreach(i = 1:nrow(grid_flexBART), .combine = rbind, .packages = packages_list) %dopar% {
    params <- grid_flexBART[i, ]
    fold_score <- numeric(5)
    for (f in 1:5) {
      cv <- get_cv_data(X_train, Y_train, cv_folds[[f]])
      
      train_cv_set <- data.frame(Y = cv$Y_train, cv$X_train)
      val_cv_set   <- data.frame(Y = cv$Y_val, cv$X_val)
      
      model<-flexBART(Y ~ bart(.), a = params$a_u, b = params$b_u, train_data = train_cv_set,
               M_vec = params$M_vec, save_trees = TRUE,alpha_vec=params$alpha,beta_vec=params$beta,nu=params$nu,sigquant=params$q,
               n.chains = 1, sparse = TRUE, burn = MCMC_BURNIN, nd = MCMC_ITER-MCMC_BURNIN)
      preds <- colMeans(predict.flexBART(model, val_cv_set))
      fold_score[f] <- if (is_classification) mean(ifelse(preds > 0.5, 1, 0) == cv$Y_val) else sqrt(mean((preds - cv$Y_val)^2))
    }
    data.frame(params, avg_score = mean(fold_score))
  }
  best_params_flex <- if (is_classification) cv_results_flex[which.max(cv_results_flex$avg_score), ] else cv_results_flex[which.min(cv_results_flex$avg_score), ]
  
  # --- TUNING softBART ---
  ##cat("Tuning softBART...\n")
  ##cv_results_soft <- foreach(i = 1:nrow(grid_softBART), .combine = rbind, .packages = packages_list) %dopar% {
  ##  params <- grid_softBART[i, ]
  ##  fold_score <- numeric(5)
  ##  for (f in 1:5) {
  ##    cv <- get_cv_data(X_train, Y_train, cv_folds[[f]])
  ##    
  ##    hypers <- Hypers(cv$X_train, cv$Y_train, 
  ##                     num_tree = params$num_tree,
  ##                     alpha = params$alpha,
  ##                     beta = params$beta)
  ##    opts <- Opts(num_burn = MCMC_BURNIN, num_save = MCMC_ITER - MCMC_BURNIN)
  ##    
  ##    if (is_classification) {
  ##      model <- softbart_probit(cv$X_train, cv$Y_train, cv$X_val, hypers = hypers, opts = opts)
  ##      preds <- model$p_mean
  ##      fold_score[f] <- mean(ifelse(preds > 0.5, 1, 0) == cv$Y_val)
  ##    } else {
  ##      model <- softbart(cv$X_train, cv$Y_train, cv$X_val, hypers = hypers, opts = opts)
  ##      preds <- model$y_hat_test_mean
  ##      fold_score[f] <- sqrt(mean((preds - cv$Y_val)^2))
  ##    }
  ##  }
  ##  data.frame(params, avg_score = mean(fold_score))
  ##}
  ##best_params_soft <- if (is_classification) cv_results_soft[which.max(cv_results_soft$avg_score), ] else cv_results_soft[which.min(cv_results_soft$avg_score), ]
  
  # --- TUNING wBART ---
  cat("Tuning wBART...\n")
  cv_results_wbart <- foreach(i = 1:nrow(grid_wBART), .combine = rbind, .packages = packages_list) %dopar% {
    params <- grid_wBART[i, ]
    fold_score <- numeric(5)
    for (f in 1:5) {
      cv <- get_cv_data(X_train, Y_train, cv_folds[[f]])
      
      if (is_classification) {
        model <- pbart(cv$X_train, cv$Y_train, cv$X_val, 
                       ntree=params$ntree, base=params$base, power=params$power,
                       ndpost=MCMC_ITER-MCMC_BURNIN, nskip=MCMC_BURNIN)
        fold_score[f] <- mean(ifelse(model$prob.test.mean > 0.5, 1, 0) == cv$Y_val)
      } else {
        model <- wbart(cv$X_train, cv$Y_train, cv$X_val, 
                       ntree=params$ntree, base=params$base, power=params$power,
                       ndpost=MCMC_ITER-MCMC_BURNIN, nskip=MCMC_BURNIN)
        fold_score[f] <- sqrt(mean((model$yhat.test.mean - cv$Y_val)^2))
      }
    }
    data.frame(params, avg_score = mean(fold_score))
  }
  best_params_wbart <- if (is_classification) cv_results_wbart[which.max(cv_results_wbart$avg_score), ] else cv_results_wbart[which.min(cv_results_wbart$avg_score), ]
  
  # --- TUNING Random Forest (ranger) ---
  cat("Tuning Random Forest...\n")
  cv_results_rf <- foreach(i = 1:nrow(grid_rf), .combine = rbind, .packages = packages_list) %dopar% {
    params <- grid_rf[i, ]
    fold_score <- numeric(5)
    
    # Dynamically assign the correct splitrule based on task and grid parameter
    if (is_classification) {
      current_splitrule <- ifelse(params$use_extratrees, "extratrees", "gini")
    } else {
      current_splitrule <- ifelse(params$use_extratrees, "extratrees", "variance")
    }
    
    for (f in 1:5) {
      cv <- get_cv_data(X_train, Y_train, cv_folds[[f]])
      mtry_val <- max(1, floor(params$mtry_prop * ncol(cv$X_train)))
      
      if (is_classification) {
        model <- ranger(x = cv$X_train, y = as.factor(cv$Y_train), 
                        num.trees = params$num.trees,
                        mtry = mtry_val, 
                        min.node.size = params$min.node.size,
                        sample.fraction = params$sample.fraction,
                        replace = params$replace,
                        splitrule = current_splitrule)
        preds <- predict(model, cv$X_val)$predictions
        fold_score[f] <- mean(as.numeric(as.character(preds)) == cv$Y_val)
      } else {
        model <- ranger(x = cv$X_train, y = cv$Y_train, 
                        num.trees = params$num.trees,
                        mtry = mtry_val, 
                        min.node.size = params$min.node.size,
                        sample.fraction = params$sample.fraction,
                        replace = params$replace,
                        splitrule = current_splitrule)
        preds <- predict(model, cv$X_val)$predictions
        fold_score[f] <- sqrt(mean((preds - cv$Y_val)^2))
      }
    }
    data.frame(params, avg_score = mean(fold_score))
  }
  best_params_rf <- if (is_classification) cv_results_rf[which.max(cv_results_rf$avg_score), ] else cv_results_rf[which.min(cv_results_rf$avg_score), ]
  
  # --- TUNING SVM ---
  cat("Tuning SVM...\n")
  cv_results_svm <- foreach(i = 1:nrow(grid_svm), .combine = rbind, .packages = packages_list) %dopar% {
    params <- grid_svm[i, ]
    fold_score <- numeric(5)
    for (f in 1:5) {
      cv <- get_cv_data(X_train, Y_train, cv_folds[[f]])
      
      if (is_classification) {
        model <- svm(x = cv$X_train, y = as.factor(cv$Y_train), 
                     kernel = as.character(params$kernel),
                     cost = params$cost, 
                     gamma = params$gamma,
                     degree = params$degree,
                     coef0 = params$coef0,
                     probability = FALSE)
        preds <- predict(model, cv$X_val)
        fold_score[f] <- mean(as.numeric(as.character(preds)) == cv$Y_val)
      } else {
        model <- svm(x = cv$X_train, y = cv$Y_train, 
                     kernel = as.character(params$kernel),
                     cost = params$cost, 
                     gamma = params$gamma,
                     degree = params$degree,
                     coef0 = params$coef0,
                     epsilon = params$epsilon)
        preds <- predict(model, cv$X_val)
        fold_score[f] <- sqrt(mean((preds - cv$Y_val)^2))
      }
    }
    data.frame(params, avg_score = mean(fold_score))
  }
  best_params_svm <- if (is_classification) cv_results_svm[which.max(cv_results_svm$avg_score), ] else cv_results_svm[which.min(cv_results_svm$avg_score), ]
  
  # --- TUNING XGBoost ---
  cat("Tuning XGBoost...\n")
  cv_results_xgb <- foreach(i = 1:nrow(grid_xgb), .combine = rbind, .packages = packages_list) %dopar% {
    params <- grid_xgb[i, ]
    fold_score <- numeric(5)
    obj_type <- ifelse(is_classification, "binary:logistic", "reg:squarederror")
    
    for (f in 1:5) {
      cv <- get_cv_data(X_train, Y_train, cv_folds[[f]])
      
      # Coerce to integer for classification to satisfy XGBoost's strict type checking
      target_label <- if (is_classification) as.factor(cv$Y_train) else cv$Y_train
      
      model <- xgboost(data = cv$X_train, 
                       label = target_label,
                       objective = obj_type,
                       nrounds = params$nrounds, 
                       max_depth = params$max_depth,
                       learning_rate = params$learning_rate,
                       subsample = params$subsample,
                       colsample_bytree = params$colsample_bytree,
                       verbose = 0)
      
      preds <- predict(model, cv$X_val)
      fold_score[f] <- if (is_classification) mean(ifelse(preds > 0.5, 1, 0) == cv$Y_val) else sqrt(mean((preds - cv$Y_val)^2))
    }
    data.frame(params, avg_score = mean(fold_score))
  }
  best_params_xgb <- if (is_classification) cv_results_xgb[which.max(cv_results_xgb$avg_score), ] else cv_results_xgb[which.min(cv_results_xgb$avg_score), ]
  
  # -------------------------------------------------------------
  cat("Fitting final models on full training set...\n")
  
  evaluation_results <- foreach(run = 1:n_runs, .combine = rbind, .packages = packages_list) %dopar% {
    
    set.seed(123 + run)
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
    
    # AddiVortes 0
    final_Addi_0 <- AddiVortes(Y_train, X_train, 
                               m = best_params_0$m, 
                               k = best_params_0$k,
                               Omega = best_params_0$Omega,
                               LambdaRate = best_params_0$LambdaRate,
                               nu = best_params_0$nu,
                               q = best_params_0$q,
                               totalMCMCIter = MCMC_ITER, mcmcBurnIn = MCMC_BURNIN, 
                               numChains = 1, IntialSigma = "LASSO", varSelMode = 0)
    preds_0 <- predict(final_Addi_0, X_test, showProgress = FALSE)
    score_addi_0 <- if (is_classification) mean(round(preds_0) == Y_test) else sqrt(mean((preds_0 - Y_test)^2))
    
    # AddiVortes 1
    final_Addi_1 <- AddiVortes(Y_train, X_train, 
                               m = best_params_1$m, 
                               k = best_params_1$k,
                               Omega = best_params_1$Omega,
                               LambdaRate = best_params_1$LambdaRate,
                               nu = best_params_1$nu,
                               q = best_params_1$q,
                               a_alpha = best_params_1$a_alpha, b_alpha = best_params_1$b_alpha, 
                               totalMCMCIter = MCMC_ITER, mcmcBurnIn = MCMC_BURNIN, dirichletWarmup = DIRICHLET_WARMUP, 
                               numChains = 1, IntialSigma = "LASSO", varSelMode = 1)
    preds_1 <- predict(final_Addi_1, X_test, showProgress = FALSE)
    score_addi_1 <- if (is_classification) mean(round(preds_1) == Y_test) else sqrt(mean((preds_1 - Y_test)^2))
    
    # AddiVortes 2
    final_Addi_2 <- AddiVortes(Y_train, X_train, 
                               m = best_params_2$m, 
                               k = best_params_2$k,
                               Omega = best_params_2$Omega,
                               LambdaRate = best_params_2$LambdaRate,
                               nu = best_params_2$nu,
                               q = best_params_2$q,
                               a_alpha = best_params_2$a_alpha, b_alpha = best_params_2$b_alpha, 
                               adaptBoost = best_params_2$adaptBoost, adaptPenalty = best_params_2$adaptPenalty, momentumDecay = best_params_2$momentumDecay,
                               totalMCMCIter = MCMC_ITER, mcmcBurnIn = MCMC_BURNIN, dirichletWarmup = DIRICHLET_WARMUP, 
                               numChains = 1, IntialSigma = "LASSO", varSelMode = 2)
    preds_2 <- predict(final_Addi_2, X_test, showProgress = FALSE)
    score_addi_2 <- if (is_classification) mean(round(preds_2) == Y_test) else sqrt(mean((preds_2 - Y_test)^2))
    
    # flexBART
    final_flex <- flexBART(Y ~ bart(.), a = best_params_flex$a, b = best_params_flex$b, train_data = training_data, M_vec = best_params_flex$M_vec, save_trees = TRUE, n.chains = 1, sparse = TRUE,alpha_vec=best_params_flex$alpha,beta_vec=best_params_flex$beta,nu=best_params_flex$nu,sigquant=best_params_flex$q, burn = MCMC_BURNIN, nd = MCMC_ITER-MCMC_BURNIN)
    preds_flex <- colMeans(predict.flexBART(final_flex, test_set[, c("Y", colnames(X_train)), drop=FALSE]))
    score_flex <- if (is_classification) mean(ifelse(preds_flex > 0.5, 1, 0) == Y_test) else sqrt(mean((preds_flex - Y_test)^2))
    
    # softBART
    ##hypers <- Hypers(X_train, Y_train, 
    ##                 num_tree = best_params_soft$num_tree,
    ##                 alpha = best_params_soft$alpha,
    ##                 beta = best_params_soft$beta)
    ##
    ##opts <- Opts(num_burn = MCMC_BURNIN, num_save = MCMC_ITER - MCMC_BURNIN)
    ##
    ##if (is_classification) {
    ##  final_soft <- softbart_probit(X_train, Y_train, X_test, hypers = hypers, opts = opts)
    ##  preds_soft <- final_soft$p_mean
    ##  score_soft <- mean(ifelse(preds_soft > 0.5, 1, 0) == Y_test)
    ##} else {
    ##  final_soft <- softbart(X_train, Y_train, X_test, hypers = hypers, opts = opts)
    ##  preds_soft <- final_soft$y_hat_test_mean
    ##  score_soft <- sqrt(mean((preds_soft - Y_test)^2))
    ##}
    
    # wBART
    if (is_classification) {
      final_wbart <- pbart(X_train, Y_train, X_test, 
                           ntree = best_params_wbart$ntree, 
                           base = best_params_wbart$base, 
                           power = best_params_wbart$power,
                           ndpost = MCMC_ITER - MCMC_BURNIN, 
                           nskip = MCMC_BURNIN)
      score_wbart <- mean(ifelse(final_wbart$prob.test.mean > 0.5, 1, 0) == Y_test)
    } else {
      final_wbart <- wbart(X_train, Y_train, X_test, 
                           ntree = best_params_wbart$ntree, 
                           base = best_params_wbart$base, 
                           power = best_params_wbart$power,
                           ndpost = MCMC_ITER - MCMC_BURNIN, 
                           nskip = MCMC_BURNIN)
      score_wbart <- sqrt(mean((final_wbart$yhat.test.mean - Y_test)^2))
    }
    
    # Random Forest
    rf_mtry <- max(1, floor(best_params_rf$mtry_prop * ncol(X_train)))
    
    if (is_classification) {
      current_splitrule <- ifelse(best_params_rf$use_extratrees, "extratrees", "gini")
      
      final_rf <- ranger(x = X_train, y = as.factor(Y_train), 
                         num.trees = best_params_rf$num.trees,
                         mtry = rf_mtry, 
                         min.node.size = best_params_rf$min.node.size,
                         sample.fraction = best_params_rf$sample.fraction,
                         replace = best_params_rf$replace,
                         splitrule = current_splitrule)
      
      preds_rf <- predict(final_rf, X_test)$predictions
      score_rf <- mean(as.numeric(as.character(preds_rf)) == Y_test)
      
    } else {
      current_splitrule <- ifelse(best_params_rf$use_extratrees, "extratrees", "variance")
      
      final_rf <- ranger(x = X_train, y = Y_train, 
                         num.trees = best_params_rf$num.trees,
                         mtry = rf_mtry, 
                         min.node.size = best_params_rf$min.node.size,
                         sample.fraction = best_params_rf$sample.fraction,
                         replace = best_params_rf$replace,
                         splitrule = current_splitrule)
      
      preds_rf <- predict(final_rf, X_test)$predictions
      score_rf <- sqrt(mean((preds_rf - Y_test)^2))
    }
    
    # SVM
    if (is_classification) {
      final_svm <- svm(x = X_train, y = as.factor(Y_train), 
                       kernel = as.character(best_params_svm$kernel),
                       cost = best_params_svm$cost, 
                       gamma = best_params_svm$gamma,
                       degree = best_params_svm$degree,
                       coef0 = best_params_svm$coef0,
                       probability = FALSE)
      preds_svm <- predict(final_svm, X_test)
      score_svm <- mean(as.numeric(as.character(preds_svm)) == Y_test)
    } else {
      final_svm <- svm(x = X_train, y = Y_train, 
                       kernel = as.character(best_params_svm$kernel),
                       cost = best_params_svm$cost, 
                       gamma = best_params_svm$gamma,
                       degree = best_params_svm$degree,
                       coef0 = best_params_svm$coef0,
                       epsilon = best_params_svm$epsilon)
      preds_svm <- predict(final_svm, X_test)
      score_svm <- sqrt(mean((preds_svm - Y_test)^2))
    }

    # XGBoost
    obj_type <- ifelse(is_classification, "binary:logistic", "reg:squarederror")
    target_label <- if (is_classification) as.factor(Y_train) else Y_train
    
    final_xgb <- xgboost(data = X_train, 
                         label = target_label,
                         objective = obj_type,
                         nrounds = best_params_xgb$nrounds, 
                         max_depth = best_params_xgb$max_depth,
                         learning_rate = best_params_xgb$learning_rate,
                         subsample = best_params_xgb$subsample,
                         colsample_bytree = best_params_xgb$colsample_bytree,
                         verbose = 0)
    
    preds_xgb <- predict(final_xgb, X_test)
    score_xgb <- if (is_classification) mean(ifelse(preds_xgb > 0.5, 1, 0) == Y_test) else sqrt(mean((preds_xgb - Y_test)^2))
    # Baseline
    if (is_classification) {
      baseline_score <- max(sum(Y_test==0), sum(Y_test==1)) / length(Y_test)
    } else {
      baseline_score <- sqrt(mean((mean(Y_test) - Y_test)^2))
    }
    
    data.frame(
      Run = run,
      AddiVortes_Mode0_Score = score_addi_0,
      AddiVortes_Mode1_Score = score_addi_1,
      AddiVortes_Mode2_Score = score_addi_2,
      flexBART_Score = score_flex,
      #softBART_Score = score_soft,
      wBART_Score = score_wbart,
      RF_Score = score_rf,
      SVM_Score = score_svm,
      XGB_Score = score_xgb,
      baseline = baseline_score
    )
  }
  
  cat(paste0("\n--- Full Results across 20 runs (", metric_name, ") ---\n"))
  print(evaluation_results)
  
  cat(paste0("\n--- Average ", metric_name, " across 20 runs ---\n"))
  print(colMeans(evaluation_results[, -1]))
  
  cat(paste0("\n--- Standard Deviation across 20 runs ---\n"))
  print(apply(evaluation_results[, -1], 2, sd))
  
  # Append to master summary dataframe
  results_summary_full <- rbind(results_summary_full, data.frame(
    Dataset = data_name,
    Metric_Type = metric_name,
    AddiVortes_Mode0_Score = mean(evaluation_results$AddiVortes_Mode0_Score),
    AddiVortes_Mode1_Score = mean(evaluation_results$AddiVortes_Mode1_Score),
    AddiVortes_Mode2_Score = mean(evaluation_results$AddiVortes_Mode2_Score),
    flexBART_Score = mean(evaluation_results$flexBART_Score),
    #softBART_Score = mean(evaluation_results$softBART_Score),
    wBART_Score = mean(evaluation_results$wBART_Score),
    RF_Score = mean(evaluation_results$RF_Score),
    SVM_Score = mean(evaluation_results$SVM_Score),
    XGB_Score = mean(evaluation_results$XGB_Score),
    baseline_score = mean(evaluation_results$baseline)
  ))
}

stopCluster(cl)
