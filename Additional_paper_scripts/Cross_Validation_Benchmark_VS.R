# ======================================================================
# RESULTS STORAGE
# ======================================================================

results_summary_full <- data.frame(
  Dataset                = character(),
  AddiVortes_Mode0_Score = numeric(),
  AddiVortes_Mode1_Score = numeric(),
  AddiVortes_Mode2_Score = numeric(),
  wBART_Score            = numeric(),
  RF_Score               = numeric(),
  SVM_Score              = numeric(),
  XGB_Score              = numeric(),
  baseline_score         = numeric(),
  stringsAsFactors = FALSE
)

# ======================================================================
# GLOBAL SETTINGS & HYPERPARAMETERS
# ======================================================================

n_runs    <- 20
num_cores <- 20

### PERFORMANCE FIX: SPLIT MCMC ITERATIONS ###
MCMC_ITER_CV        <- 3000
MCMC_BURNIN_CV      <- 1500
DIRICHLET_WARMUP_CV <- 400

MCMC_ITER_FINAL        <- 5000
MCMC_BURNIN_FINAL      <- 2500
DIRICHLET_WARMUP_FINAL <- 1000

# ======================================================================
# PARALLEL SETUP
# ======================================================================

library(doParallel)

set.seed(839)

# Detect OS and use FORK on Unix systems for Copy-on-Write shared memory.
# Fallback to PSOCK for Windows.
cluster_type <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
cl <- makeCluster(num_cores, type = cluster_type)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = 839)

clusterEvalQ(cl, {
  devtools::load_all()
  library(BART)
  library(ranger)
  library(e1071)
  library(xgboost)
})

# Export global MCMC constants once — never re-exported by foreach
clusterExport(cl,
              varlist = c("MCMC_ITER_CV",    "MCMC_BURNIN_CV",    "DIRICHLET_WARMUP_CV",
                          "MCMC_ITER_FINAL", "MCMC_BURNIN_FINAL", "DIRICHLET_WARMUP_FINAL"),
              envir = environment())

# ======================================================================
# CV HYPERPARAMETER GRIDS
# ======================================================================
{
  # --- AddiVortes Mode 0 ---
  base_grid_0 <- expand.grid(
    m          = c(200),
    k          = c(1, 3),
    Omega      = c(0.5, 1.5),
    LambdaRate = c(2, 4)
  )
  nu_q_cases  <- data.frame(nu = c(3, 6), q = c(0.99, 0.85))
  grid_mode0  <- merge(base_grid_0, nu_q_cases, by = NULL)
  
  # --- AddiVortes Mode 1 ---
  base_grid_1 <- expand.grid(
    m = c(200), k = c(1, 3), Omega = c(1.5, 3),
    LambdaRate = c(4), nu = c(6), q = c(0.85)
  )
  grid_mode_alpha <- data.frame(
    alpha = c(1, 20, 50), a_alpha = c(0.5, 20, 50), b_alpha = c(1.0)
  )
  grid_mode1 <- merge(base_grid_1, grid_mode_alpha, by = NULL)
  
  # --- AddiVortes Mode 2 ---
  grid_mode_usual <- expand.grid(
    m = c(200), k = c(2), Omega = c(1.5), LambdaRate = c(4),
    tau = c(20), momentumDecay = c(0.9)
  )
  grid_mode_comb       <- merge(grid_mode_usual, grid_mode_alpha, by = NULL)
  hyperparam_boost_penalty <- data.frame(
    adaptBoost = c(0.5, 2, 4, 10), adaptPenalty = c(0.5,2, 4, 10),
    stringsAsFactors = FALSE
  )
  grid_mode2 <- merge(grid_mode_comb, hyperparam_boost_penalty, by = NULL)
  
  # --- wBART ---
  grid_wBART <- expand.grid(
    ntree = c(20, 100), base = c(0.95, 0.99), power = c(2.0, 3.0),
    a = c(0.5, 1.0), b = 1.0, k = c(1.0, 2.0, 3.0)
  )
  
  # --- Random Forest ---
  grid_rf <- expand.grid(
    mtry_prop       = c(0.2, 0.4, 0.6, 0.8),
    min.node.size   = c(1, 3, 5, 10),
    sample.fraction = c(0.632, 0.8),
    num.trees       = c(250, 500),
    replace         = c(TRUE, FALSE),
    use_extratrees  = c(TRUE, FALSE)
  )
  
  # --- SVM ---
  grid_svm_radial <- expand.grid(
    kernel = "radial", cost = c(0.1, 1, 10, 50),
    gamma = c(0.001, 0.01, 0.1, 0.5), degree = 3, coef0 = 0,
    epsilon = c(0.01, 0.1)
  )
  grid_svm_poly <- expand.grid(
    kernel = "polynomial", cost = c(0.1, 1, 10, 50),
    gamma = c(0.01, 0.1), degree = c(2, 3), coef0 = c(0, 1),
    epsilon = c(0.01, 0.1)
  )
  grid_svm <- rbind(grid_svm_radial, grid_svm_poly)
  
  # --- XGBoost ---
  grid_xgb <- expand.grid(
    nrounds = c(100, 300), max_depth = c(3, 5, 7),
    learning_rate = c(0.01, 0.05, 0.1),
    subsample = c(0.7, 1.0), colsample_bytree = c(0.7, 1.0)
  )
}

# ======================================================================
# HELPER FUNCTIONS
# ======================================================================

createFolds <- function(y, k = 5, returnTrain = TRUE) {
  n                <- length(y)
  shuffled_indices <- sample(1:n)
  folds            <- split(shuffled_indices,
                            cut(seq_along(shuffled_indices), breaks = k, labels = FALSE))
  if (returnTrain) {
    train_folds        <- lapply(folds, function(test_idx) setdiff(1:n, test_idx))
    names(train_folds) <- paste0("Fold", 1:k)
    return(train_folds)
  } else {
    names(folds) <- paste0("Fold", 1:k)
    return(folds)
  }
}

get_cv_data <- function(X_train, Y_train, train_idx) {
  X_cv_train <- X_train[train_idx, , drop = FALSE]
  X_cv_val   <- X_train[-train_idx, , drop = FALSE]
  Y_cv_train <- Y_train[train_idx]
  Y_cv_val   <- Y_train[-train_idx]
  non_const  <- apply(X_cv_train, 2, function(col) length(unique(col)) > 1)
  list(
    X_train = X_cv_train[, non_const, drop = FALSE],
    X_val   = X_cv_val[,  non_const, drop = FALSE],
    Y_train = Y_cv_train,
    Y_val   = Y_cv_val
  )
}

# ======================================================================
# MAIN BENCHMARK LOOP
# ======================================================================

for (data_name in names(benchmark_datasets)[9]) {
  
  # Proactively force workers to clear memory from the previous dataset
  clusterEvalQ(cl, gc())
  
  cat("\nAnalysing dataset:", data_name, "\n")
  current_data <- benchmark_datasets[[data_name]]
  
  Y <- current_data$Y
  X <- current_data$X
  
  is_classification <- length(unique(na.omit(Y))) == 2
  metric_name       <- ifelse(is_classification, "Accuracy", "RMSE")
  cat("Detected task:  ", ifelse(is_classification, "Binary Classification", "Regression"), "\n")
  cat("Target Metric:  ", metric_name, "\n\n")
  
  if (is_classification) {
    unique_vals <- sort(unique(na.omit(Y)))
    Y <- as.integer(ifelse(Y == unique_vals[2], 1L, 0L))
  }
  
  data_combined <- data.frame(Y = Y, X)
  n             <- nrow(data_combined)
  train_indices <- sample(1:n, size = round(0.8 * n))
  
  train_set <- data_combined[train_indices, ]
  test_set  <- data_combined[-train_indices, ]
  
  X_train <- model.matrix(Y ~ . - 1, data = train_set)
  X_test  <- model.matrix(Y ~ . - 1, data = test_set)
  Y_train <- train_set$Y
  Y_test  <- test_set$Y
  
  non_constant_cols <- apply(X_train, 2, function(col) length(unique(col)) > 1)
  X_train <- X_train[, non_constant_cols, drop = FALSE]
  X_test  <- X_test[,  non_constant_cols, drop = FALSE]
  
  training_data       <- data.frame(Y = Y_train, X_train)
  cv_folds            <- createFolds(Y_train, k = 5, returnTrain = TRUE)
  precomputed_cv_data <- lapply(1:5, function(f) get_cv_data(X_train, Y_train, cv_folds[[f]]))
  
  packages_list <- c("glmnet", "BART", "ranger", "e1071", "xgboost")
  
  # Export per-dataset variables once so foreach never re-serialises them
  clusterExport(cl,
                varlist = c("precomputed_cv_data", "is_classification", "n", "data_combined"),
                envir   = environment())
  
  # CV blacklist — everything already on workers or too large to ship
  heavy_objects_cv <- c(
    "benchmark_datasets", "data_combined", "X_train", "X_test",
    "Y_train", "Y_test", "training_data", "train_set", "test_set",
    "X", "Y", "current_data",
    "precomputed_cv_data", "is_classification", "n",
    "MCMC_ITER_CV", "MCMC_BURNIN_CV", "DIRICHLET_WARMUP_CV"
  )
  
  # --------------------------------------------------------------------
  # TUNING: AddiVortes Mode 0  (param × fold flattened)
  # --------------------------------------------------------------------
  cat(sprintf("Tuning AddiVortes Mode 0 (%d combinations × 5 folds)...\n", nrow(grid_mode0)))
  pf_grid_0 <- expand.grid(param_idx = seq_len(nrow(grid_mode0)), fold = 1:5)
  
  raw_0 <- foreach(row = seq_len(nrow(pf_grid_0)), .combine = rbind,
                   .packages = packages_list, .noexport = heavy_objects_cv) %dopar% {
                     i      <- pf_grid_0$param_idx[row]; f <- pf_grid_0$fold[row]
                     params <- grid_mode0[i, , drop = FALSE]; cv <- precomputed_cv_data[[f]]
                     model  <- AddiVortes(cv$Y_train, cv$X_train,
                                          m = params$m, k = params$k, Omega = params$Omega,
                                          LambdaRate = params$LambdaRate, nu = params$nu, q = params$q,
                                          totalMCMCIter = MCMC_ITER_CV, mcmcBurnIn = MCMC_BURNIN_CV,
                                          numChains = 1, showProgress = FALSE,
                                          IntialSigma = "LASSO", varSelMode = 0)
                     preds <- predict(model, cv$X_val, showProgress = FALSE)
                     score <- if (is_classification) mean(round(preds) == cv$Y_val) else
                       sqrt(mean((preds - cv$Y_val)^2))
                     
                     # IMMEDIATE REMOVAL AND GC
                     rm(model, cv, preds)
                     gc()
                     
                     data.frame(param_idx = i, fold = f, score = score)
                   }
  cv_agg_0     <- aggregate(score ~ param_idx, data = raw_0, FUN = mean)
  cv_results_0 <- merge(cv_agg_0, cbind(param_idx = seq_len(nrow(grid_mode0)), grid_mode0), by = "param_idx")
  names(cv_results_0)[names(cv_results_0) == "score"] <- "avg_score"
  best_params_0 <- if (is_classification) cv_results_0[which.max(cv_results_0$avg_score), ] else
    cv_results_0[which.min(cv_results_0$avg_score), ]
  print(best_params_0)
  
  # --------------------------------------------------------------------
  # TUNING: AddiVortes Mode 1  (param × fold flattened)
  # --------------------------------------------------------------------
  cat(sprintf("Tuning AddiVortes Mode 1 (%d combinations × 5 folds)...\n", nrow(grid_mode1)))
  pf_grid_1 <- expand.grid(param_idx = seq_len(nrow(grid_mode1)), fold = 1:5)
  
  raw_1 <- foreach(row = seq_len(nrow(pf_grid_1)), .combine = rbind,
                   .packages = packages_list, .noexport = heavy_objects_cv) %dopar% {
                     i      <- pf_grid_1$param_idx[row]; f <- pf_grid_1$fold[row]
                     params <- grid_mode1[i, , drop = FALSE]; cv <- precomputed_cv_data[[f]]
                     model  <- AddiVortes(cv$Y_train, cv$X_train,
                                          m = params$m, k = params$k, Omega = params$Omega,
                                          LambdaRate = params$LambdaRate, nu = params$nu, q = params$q,
                                          a_alpha = params$a_alpha, b_alpha = params$b_alpha,
                                          totalMCMCIter = MCMC_ITER_CV, mcmcBurnIn = MCMC_BURNIN_CV,
                                          dirichletWarmup = DIRICHLET_WARMUP_CV,
                                          numChains = 1, showProgress = FALSE,
                                          IntialSigma = "LASSO", varSelMode = 1)
                     preds <- predict(model, cv$X_val, showProgress = FALSE)
                     score <- if (is_classification) mean(round(preds) == cv$Y_val) else
                       sqrt(mean((preds - cv$Y_val)^2))
                     
                     # IMMEDIATE REMOVAL AND GC
                     rm(model, cv, preds)
                     gc()
                     
                     data.frame(param_idx = i, fold = f, score = score)
                   }
  cv_agg_1     <- aggregate(score ~ param_idx, data = raw_1, FUN = mean)
  cv_results_1 <- merge(cv_agg_1, cbind(param_idx = seq_len(nrow(grid_mode1)), grid_mode1), by = "param_idx")
  names(cv_results_1)[names(cv_results_1) == "score"] <- "avg_score"
  best_params_1 <- if (is_classification) cv_results_1[which.max(cv_results_1$avg_score), ] else
    cv_results_1[which.min(cv_results_1$avg_score), ]
  
  print(best_params_1)
  # --------------------------------------------------------------------
  # TUNING: AddiVortes Mode 2  (param × fold flattened)
  # --------------------------------------------------------------------
  cat(sprintf("Tuning AddiVortes Mode 2 (%d combinations × 5 folds)...\n", nrow(grid_mode2)))
  pf_grid_2 <- expand.grid(param_idx = seq_len(nrow(grid_mode2)), fold = 1:5)
  
  raw_2 <- foreach(row = seq_len(nrow(pf_grid_2)), .combine = rbind,
                   .packages = packages_list, .noexport = heavy_objects_cv) %dopar% {
                     i      <- pf_grid_2$param_idx[row]; f <- pf_grid_2$fold[row]
                     params <- grid_mode2[i, ]; cv <- precomputed_cv_data[[f]]
                     model  <- AddiVortes(cv$Y_train, cv$X_train,
                                          m = params$m, k = params$k, Omega = params$Omega,
                                          LambdaRate = params$LambdaRate,
                                          a_alpha = params$a_alpha, b_alpha = params$b_alpha,
                                          adaptBoost = params$adaptBoost, adaptPenalty = params$adaptPenalty,
                                          momentumDecay = params$momentumDecay, tau = params$tau,
                                          totalMCMCIter = MCMC_ITER_CV, mcmcBurnIn = MCMC_BURNIN_CV,
                                          dirichletWarmup = DIRICHLET_WARMUP_CV,
                                          numChains = 1, showProgress = FALSE,
                                          IntialSigma = "LASSO", varSelMode = 2)
                     preds <- predict(model, cv$X_val, showProgress = FALSE)
                     score <- if (is_classification) mean(round(preds) == cv$Y_val) else
                       sqrt(mean((preds - cv$Y_val)^2))
                     
                     # IMMEDIATE REMOVAL AND GC
                     rm(model, cv, preds)
                     gc()
                     
                     data.frame(param_idx = i, fold = f, score = score)
                   }
  cv_agg_2     <- aggregate(score ~ param_idx, data = raw_2, FUN = mean)
  cv_results_2 <- merge(cv_agg_2, cbind(param_idx = seq_len(nrow(grid_mode2)), grid_mode2), by = "param_idx")
  names(cv_results_2)[names(cv_results_2) == "score"] <- "avg_score"
  best_params_2 <- if (is_classification) cv_results_2[which.max(cv_results_2$avg_score), ] else
    cv_results_2[which.min(cv_results_2$avg_score), ]
  print(best_params_2)
  # --------------------------------------------------------------------
  # TUNING: wBART  (param × fold flattened)
  # --------------------------------------------------------------------
  cat(sprintf("Tuning wBART (%d combinations × 5 folds)...\n", nrow(grid_wBART)))
  pf_grid_wbart <- expand.grid(param_idx = seq_len(nrow(grid_wBART)), fold = 1:5)
  
  raw_wbart <- foreach(row = seq_len(nrow(pf_grid_wbart)), .combine = rbind,
                       .packages = packages_list, .noexport = heavy_objects_cv) %dopar% {
                         i      <- pf_grid_wbart$param_idx[row]; f <- pf_grid_wbart$fold[row]
                         params <- grid_wBART[i, ]; cv <- precomputed_cv_data[[f]]
                         if (is_classification) {
                           model <- pbart(cv$X_train, cv$Y_train, cv$X_val,
                                          ntree = params$ntree, base = params$base, power = params$power,
                                          sparse = TRUE, a = params$a, k = params$k,
                                          ndpost = MCMC_ITER_CV - MCMC_BURNIN_CV, nskip = MCMC_BURNIN_CV,
                                          nkeeptreedraws = 0, printevery = .Machine$integer.max)
                           score <- mean(ifelse(model$prob.test.mean > 0.5, 1, 0) == cv$Y_val)
                         } else {
                           model <- wbart(cv$X_train, cv$Y_train, cv$X_val,
                                          ntree = params$ntree, base = params$base, power = params$power,
                                          sparse = TRUE, a = params$a, k = params$k,
                                          ndpost = MCMC_ITER_CV - MCMC_BURNIN_CV, nskip = MCMC_BURNIN_CV,
                                          nkeeptreedraws = 0, printevery = .Machine$integer.max)
                           score <- sqrt(mean((model$yhat.test.mean - cv$Y_val)^2))
                         }
                         
                         # IMMEDIATE REMOVAL AND GC
                         rm(model, cv)
                         gc()
                         
                         data.frame(param_idx = i, fold = f, score = score)
                       }
  cv_agg_wbart     <- aggregate(score ~ param_idx, data = raw_wbart, FUN = mean)
  cv_results_wbart <- merge(cv_agg_wbart, cbind(param_idx = seq_len(nrow(grid_wBART)), grid_wBART), by = "param_idx")
  names(cv_results_wbart)[names(cv_results_wbart) == "score"] <- "avg_score"
  best_params_wbart <- if (is_classification) cv_results_wbart[which.max(cv_results_wbart$avg_score), ] else
    cv_results_wbart[which.min(cv_results_wbart$avg_score), ]
  
  # --------------------------------------------------------------------
  # TUNING: Random Forest  (param × fold flattened)
  # --------------------------------------------------------------------
  cat(sprintf("Tuning Random Forest (%d combinations × 5 folds)...\n", nrow(grid_rf)))
  pf_grid_rf <- expand.grid(param_idx = seq_len(nrow(grid_rf)), fold = 1:5)
  
  raw_rf <- foreach(row = seq_len(nrow(pf_grid_rf)), .combine = rbind,
                    .packages = packages_list, .noexport = heavy_objects_cv) %dopar% {
                      i      <- pf_grid_rf$param_idx[row]; f <- pf_grid_rf$fold[row]
                      params <- grid_rf[i, ]; cv <- precomputed_cv_data[[f]]
                      current_splitrule <- if (is_classification) ifelse(params$use_extratrees, "extratrees", "gini") else
                        ifelse(params$use_extratrees, "extratrees", "variance")
                      mtry_val <- max(1, floor(params$mtry_prop * ncol(cv$X_train)))
                      
                      if (is_classification) {
                        model <- ranger(x = cv$X_train, y = as.factor(cv$Y_train),
                                        num.trees = params$num.trees, mtry = mtry_val,
                                        min.node.size = params$min.node.size, sample.fraction = params$sample.fraction,
                                        replace = params$replace, splitrule = current_splitrule,
                                        importance = "none", num.threads = 1)
                        preds <- predict(model, cv$X_val)$predictions
                        score <- mean(as.numeric(as.character(preds)) == cv$Y_val)
                      } else {
                        model <- ranger(x = cv$X_train, y = cv$Y_train,
                                        num.trees = params$num.trees, mtry = mtry_val,
                                        min.node.size = params$min.node.size, sample.fraction = params$sample.fraction,
                                        replace = params$replace, splitrule = current_splitrule,
                                        importance = "none", num.threads = 1)
                        preds <- predict(model, cv$X_val)$predictions
                        score <- sqrt(mean((preds - cv$Y_val)^2))
                      }
                      
                      # IMMEDIATE REMOVAL AND GC
                      rm(model, cv, preds)
                      gc()
                      
                      data.frame(param_idx = i, fold = f, score = score)
                    }
  cv_agg_rf     <- aggregate(score ~ param_idx, data = raw_rf, FUN = mean)
  cv_results_rf <- merge(cv_agg_rf, cbind(param_idx = seq_len(nrow(grid_rf)), grid_rf), by = "param_idx")
  names(cv_results_rf)[names(cv_results_rf) == "score"] <- "avg_score"
  best_params_rf <- if (is_classification) cv_results_rf[which.max(cv_results_rf$avg_score), ] else
    cv_results_rf[which.min(cv_results_rf$avg_score), ]
  
  # --------------------------------------------------------------------
  # TUNING: SVM  (param × fold flattened)
  # --------------------------------------------------------------------
  cat(sprintf("Tuning SVM (%d combinations × 5 folds)...\n", nrow(grid_svm)))
  pf_grid_svm <- expand.grid(param_idx = seq_len(nrow(grid_svm)), fold = 1:5)
  
  raw_svm <- foreach(row = seq_len(nrow(pf_grid_svm)), .combine = rbind,
                     .packages = packages_list, .noexport = heavy_objects_cv) %dopar% {
                       i      <- pf_grid_svm$param_idx[row]; f <- pf_grid_svm$fold[row]
                       params <- grid_svm[i, ]; cv <- precomputed_cv_data[[f]]
                       if (is_classification) {
                         model <- svm(x = cv$X_train, y = as.factor(cv$Y_train),
                                      kernel = as.character(params$kernel), cost = params$cost,
                                      gamma = params$gamma, degree = params$degree, coef0 = params$coef0,
                                      probability = FALSE)
                         preds <- predict(model, cv$X_val)
                         score <- mean(as.numeric(as.character(preds)) == cv$Y_val)
                       } else {
                         model <- svm(x = cv$X_train, y = cv$Y_train,
                                      kernel = as.character(params$kernel), cost = params$cost,
                                      gamma = params$gamma, degree = params$degree,
                                      coef0 = params$coef0, epsilon = params$epsilon)
                         preds <- predict(model, cv$X_val)
                         score <- sqrt(mean((preds - cv$Y_val)^2))
                       }
                       
                       # IMMEDIATE REMOVAL AND GC
                       rm(model, cv, preds)
                       gc()
                       
                       data.frame(param_idx = i, fold = f, score = score)
                     }
  cv_agg_svm     <- aggregate(score ~ param_idx, data = raw_svm, FUN = mean)
  cv_results_svm <- merge(cv_agg_svm, cbind(param_idx = seq_len(nrow(grid_svm)), grid_svm), by = "param_idx")
  names(cv_results_svm)[names(cv_results_svm) == "score"] <- "avg_score"
  best_params_svm <- if (is_classification) cv_results_svm[which.max(cv_results_svm$avg_score), ] else
    cv_results_svm[which.min(cv_results_svm$avg_score), ]
  
  # --------------------------------------------------------------------
  # TUNING: XGBoost  (param × fold flattened, device = "cpu")
  # --------------------------------------------------------------------
  cat(sprintf("Tuning XGBoost (%d combinations × 5 folds)...\n", nrow(grid_xgb)))
  pf_grid_xgb <- expand.grid(param_idx = seq_len(nrow(grid_xgb)), fold = 1:5)
  
  raw_xgb <- foreach(row = seq_len(nrow(pf_grid_xgb)), .combine = rbind,
                     .packages = packages_list, .noexport = heavy_objects_cv) %dopar% {
                       i        <- pf_grid_xgb$param_idx[row]; f <- pf_grid_xgb$fold[row]
                       params   <- grid_xgb[i, ]; cv <- precomputed_cv_data[[f]]
                       obj_type <- ifelse(is_classification, "binary:logistic", "reg:squarederror")
                       model    <- xgboost(data = cv$X_train,
                                           label = if (is_classification) as.integer(cv$Y_train) else cv$Y_train,
                                           objective = obj_type, nrounds = params$nrounds,
                                           max_depth = params$max_depth, learning_rate = params$learning_rate,
                                           subsample = params$subsample, colsample_bytree = params$colsample_bytree,
                                           tree_method = "hist", device = "cpu", nthread = 1, verbose = 0)
                       preds <- predict(model, cv$X_val)
                       score <- if (is_classification) mean(ifelse(preds > 0.5, 1, 0) == cv$Y_val) else
                         sqrt(mean((preds - cv$Y_val)^2))
                       
                       # IMMEDIATE REMOVAL AND GC
                       rm(model, cv, preds)
                       gc()
                       
                       data.frame(param_idx = i, fold = f, score = score)
                     }
  cv_agg_xgb     <- aggregate(score ~ param_idx, data = raw_xgb, FUN = mean)
  cv_results_xgb <- merge(cv_agg_xgb, cbind(param_idx = seq_len(nrow(grid_xgb)), grid_xgb), by = "param_idx")
  names(cv_results_xgb)[names(cv_results_xgb) == "score"] <- "avg_score"
  best_params_xgb <- if (is_classification) cv_results_xgb[which.max(cv_results_xgb$avg_score), ] else
    cv_results_xgb[which.min(cv_results_xgb$avg_score), ]
  
  # ====================================================================
  # FINAL EVALUATION
  # ====================================================================
  cat("Running final evaluation (20 runs × 8 models = 160 atomic jobs)...\n")
  
  # Step 1: Pre-generate all train/test split indices
  final_splits <- lapply(seq_len(n_runs), function(r) sample(1:n, size = round(0.8 * n)))
  
  # Step 2: Export best params + splits to workers once
  clusterExport(cl,
                varlist = c("final_splits",
                            "best_params_0", "best_params_1", "best_params_2",
                            "best_params_wbart", "best_params_rf",
                            "best_params_svm", "best_params_xgb"),
                envir = environment())
  
  # Step 3: Build (run × model) job grid
  model_names    <- c("addi0", "addi1", "addi2", "wbart", "rf", "svm", "xgb", "baseline")
  job_grid_final <- expand.grid(run = seq_len(n_runs), model = model_names,
                                stringsAsFactors = FALSE)
  
  # Comprehensive blacklist — everything already on workers via clusterExport
  heavy_objects_final <- c(
    "benchmark_datasets", "precomputed_cv_data", "cv_folds",
    "X_train", "X_test", "Y_train", "Y_test", "training_data",
    "train_set", "test_set", "X", "Y", "current_data",
    "data_combined", "n", "is_classification",
    "MCMC_ITER_FINAL", "MCMC_BURNIN_FINAL", "DIRICHLET_WARMUP_FINAL",
    "final_splits",
    "best_params_0", "best_params_1", "best_params_2",
    "best_params_wbart", "best_params_rf", "best_params_svm", "best_params_xgb"
  )
  
  # Step 4: Run — each job = one model fit, one score returned
  raw_final <- foreach(
    job       = seq_len(nrow(job_grid_final)),
    .combine  = rbind,
    .packages = packages_list,
    .noexport = heavy_objects_final
  ) %dopar% {
    
    run_id     <- job_grid_final$run[job]
    model_name <- job_grid_final$model[job]
    
    # Reconstruct this run's train/test split
    train_idx <- final_splits[[run_id]]
    train_set <- data_combined[train_idx, ]
    test_set  <- data_combined[-train_idx, ]
    X_tr      <- model.matrix(Y ~ . - 1, data = train_set)
    X_te      <- model.matrix(Y ~ . - 1, data = test_set)
    Y_tr      <- train_set$Y
    Y_te      <- test_set$Y
    
    non_const <- apply(X_tr, 2, function(col) length(unique(col)) > 1)
    X_tr <- X_tr[, non_const, drop = FALSE]
    X_te <- X_te[, non_const, drop = FALSE]
    
    # Fit the single assigned model, compute score, and IMMEDIATELY REMOVE object
    score_val <- if (model_name == "addi0") {
      
      m <- AddiVortes(Y_tr, X_tr,
                      m = best_params_0$m, k = best_params_0$k,
                      Omega = best_params_0$Omega, LambdaRate = best_params_0$LambdaRate,
                      nu = best_params_0$nu, q = best_params_0$q,
                      totalMCMCIter = MCMC_ITER_FINAL, mcmcBurnIn = MCMC_BURNIN_FINAL,
                      numChains = 1, showProgress = FALSE, IntialSigma = "LASSO", varSelMode = 0)
      p <- predict(m, X_te, showProgress = FALSE)
      
      # IMMEDIATE REMOVAL
      rm(m)
      gc()
      
      if (is_classification) mean(round(p) == Y_te) else sqrt(mean((p - Y_te)^2))
      
    } else if (model_name == "addi1") {
      
      m <- AddiVortes(Y_tr, X_tr,
                      m = best_params_1$m, k = best_params_1$k,
                      Omega = best_params_1$Omega, LambdaRate = best_params_1$LambdaRate,
                      nu = best_params_1$nu, q = best_params_1$q,
                      a_alpha = best_params_1$a_alpha, b_alpha = best_params_1$b_alpha,
                      totalMCMCIter = MCMC_ITER_FINAL, mcmcBurnIn = MCMC_BURNIN_FINAL,
                      dirichletWarmup = DIRICHLET_WARMUP_FINAL,
                      numChains = 1, showProgress = FALSE, IntialSigma = "LASSO", varSelMode = 1)
      p <- predict(m, X_te, showProgress = FALSE)
      
      # IMMEDIATE REMOVAL
      rm(m)
      gc()
      
      if (is_classification) mean(round(p) == Y_te) else sqrt(mean((p - Y_te)^2))
      
    } else if (model_name == "addi2") {
      
      m <- AddiVortes(Y_tr, X_tr,
                      m = best_params_2$m, k = best_params_2$k,
                      Omega = best_params_2$Omega, LambdaRate = best_params_2$LambdaRate,
                      a_alpha = best_params_2$a_alpha, b_alpha = best_params_2$b_alpha,
                      tau = best_params_2$tau,
                      adaptBoost = best_params_2$adaptBoost, adaptPenalty = best_params_2$adaptPenalty,
                      momentumDecay = best_params_2$momentumDecay,
                      totalMCMCIter = MCMC_ITER_FINAL, mcmcBurnIn = MCMC_BURNIN_FINAL,
                      dirichletWarmup = DIRICHLET_WARMUP_FINAL,
                      numChains = 1, showProgress = FALSE, IntialSigma = "LASSO", varSelMode = 2)
      p <- predict(m, X_te, showProgress = FALSE)
      
      # IMMEDIATE REMOVAL
      rm(m)
      gc()
      
      if (is_classification) mean(round(p) == Y_te) else sqrt(mean((p - Y_te)^2))
      
    } else if (model_name == "wbart") {
      
      if (is_classification) {
        m <- pbart(X_tr, Y_tr, X_te,
                   ntree = best_params_wbart$ntree, base = best_params_wbart$base,
                   power = best_params_wbart$power, a = best_params_wbart$a,
                   k = best_params_wbart$k, sparse = TRUE,
                   ndpost = MCMC_ITER_FINAL - MCMC_BURNIN_FINAL, nskip = MCMC_BURNIN_FINAL,
                   nkeeptreedraws = 0, printevery = .Machine$integer.max)
        val <- mean(ifelse(m$prob.test.mean > 0.5, 1, 0) == Y_te)
        
        # IMMEDIATE REMOVAL
        rm(m)
        gc()
        val
      } else {
        m <- wbart(X_tr, Y_tr, X_te,
                   ntree = best_params_wbart$ntree, base = best_params_wbart$base,
                   power = best_params_wbart$power, a = best_params_wbart$a,
                   k = best_params_wbart$k, sparse = TRUE,
                   ndpost = MCMC_ITER_FINAL - MCMC_BURNIN_FINAL, nskip = MCMC_BURNIN_FINAL,
                   nkeeptreedraws = 0, printevery = .Machine$integer.max)
        val <- sqrt(mean((m$yhat.test.mean - Y_te)^2))
        
        # IMMEDIATE REMOVAL
        rm(m)
        gc()
        val
      }
      
    } else if (model_name == "rf") {
      
      rf_mtry <- max(1, floor(best_params_rf$mtry_prop * ncol(X_tr)))
      if (is_classification) {
        current_splitrule <- ifelse(best_params_rf$use_extratrees, "extratrees", "gini")
        m <- ranger(x = X_tr, y = as.factor(Y_tr),
                    num.trees = best_params_rf$num.trees, mtry = rf_mtry,
                    min.node.size = best_params_rf$min.node.size,
                    sample.fraction = best_params_rf$sample.fraction,
                    replace = best_params_rf$replace, splitrule = current_splitrule,
                    importance = "none", num.threads = 1)
        p <- predict(m, X_te)$predictions
        
        # IMMEDIATE REMOVAL
        rm(m)
        gc()
        
        mean(as.numeric(as.character(p)) == Y_te)
      } else {
        current_splitrule <- ifelse(best_params_rf$use_extratrees, "extratrees", "variance")
        m <- ranger(x = X_tr, y = Y_tr,
                    num.trees = best_params_rf$num.trees, mtry = rf_mtry,
                    min.node.size = best_params_rf$min.node.size,
                    sample.fraction = best_params_rf$sample.fraction,
                    replace = best_params_rf$replace, splitrule = current_splitrule,
                    importance = "none", num.threads = 1)
        p <- predict(m, X_te)$predictions
        
        # IMMEDIATE REMOVAL
        rm(m)
        gc()
        
        sqrt(mean((p - Y_te)^2))
      }
      
    } else if (model_name == "svm") {
      
      if (is_classification) {
        m <- svm(x = X_tr, y = as.factor(Y_tr),
                 kernel = as.character(best_params_svm$kernel), cost = best_params_svm$cost,
                 gamma = best_params_svm$gamma, degree = best_params_svm$degree,
                 coef0 = best_params_svm$coef0, probability = FALSE)
        p <- predict(m, X_te)
        
        # IMMEDIATE REMOVAL
        rm(m)
        gc()
        
        mean(as.numeric(as.character(p)) == Y_te)
      } else {
        m <- svm(x = X_tr, y = Y_tr,
                 kernel = as.character(best_params_svm$kernel), cost = best_params_svm$cost,
                 gamma = best_params_svm$gamma, degree = best_params_svm$degree,
                 coef0 = best_params_svm$coef0, epsilon = best_params_svm$epsilon)
        p <- predict(m, X_te)
        
        # IMMEDIATE REMOVAL
        rm(m)
        gc()
        
        sqrt(mean((p - Y_te)^2))
      }
      
    } else if (model_name == "xgb") {
      
      obj_type <- ifelse(is_classification, "binary:logistic", "reg:squarederror")
      m <- xgboost(data  = X_tr,
                   label = if (is_classification) as.integer(Y_tr) else Y_tr,
                   objective = obj_type,
                   nrounds = best_params_xgb$nrounds, max_depth = best_params_xgb$max_depth,
                   learning_rate = best_params_xgb$learning_rate,
                   subsample = best_params_xgb$subsample,
                   colsample_bytree = best_params_xgb$colsample_bytree,
                   tree_method = "hist", device = "cpu", nthread = 1, verbose = 0)
      p <- predict(m, X_te)
      
      # IMMEDIATE REMOVAL
      rm(m)
      gc()
      
      if (is_classification) mean(ifelse(p > 0.5, 1, 0) == Y_te) else sqrt(mean((p - Y_te)^2))
      
    } else { # baseline
      
      if (is_classification) max(sum(Y_te == 0), sum(Y_te == 1)) / length(Y_te) else
        sqrt(mean((mean(Y_te) - Y_te)^2))
      
    }
    
    # Return the clean scalar score
    data.frame(run = run_id, model = model_name, score = score_val)
  }
  
  # Step 5: Reshape long -> wide (one row per run)
  evaluation_results <- reshape(raw_final, idvar = "run", timevar = "model", direction = "wide")
  # Rename score.* columns to match original conventions
  col_map <- c(
    "score.addi0"    = "AddiVortes_Mode0_Score",
    "score.addi1"    = "AddiVortes_Mode1_Score",
    "score.addi2"    = "AddiVortes_Mode2_Score",
    "score.wbart"    = "wBART_Score",
    "score.rf"       = "RF_Score",
    "score.svm"      = "SVM_Score",
    "score.xgb"      = "XGB_Score",
    "score.baseline" = "baseline"
  )
  names(evaluation_results) <- c("Run", col_map[names(evaluation_results)[-1]])
  
  cat(paste0("\n--- Average ", metric_name, " across 20 runs ---\n"))
  print(colMeans(evaluation_results[, -1]))
  
  cat(paste0("\n--- Standard Deviation across 20 runs ---\n"))
  print(apply(evaluation_results[, -1], 2, sd))
  
  results_summary_full <- rbind(results_summary_full, data.frame(
    Dataset                = data_name,
    AddiVortes_Mode0_Score = mean(evaluation_results$AddiVortes_Mode0_Score),
    AddiVortes_Mode1_Score = mean(evaluation_results$AddiVortes_Mode1_Score),
    AddiVortes_Mode2_Score = mean(evaluation_results$AddiVortes_Mode2_Score),
    wBART_Score            = mean(evaluation_results$wBART_Score),
    RF_Score               = mean(evaluation_results$RF_Score),
    SVM_Score              = mean(evaluation_results$SVM_Score),
    XGB_Score              = mean(evaluation_results$XGB_Score),
    baseline_score         = mean(evaluation_results$baseline)
  ))
}

stopCluster(cl)

# ======================================================================
# RANKING LOGIC
# ======================================================================
algorithm_names <- c("AddiVortes_Mode0", "AddiVortes_Mode1", "AddiVortes_Mode2",
                     "wBART", "RF", "SVM", "XGBoost", "Baseline")
average_spot    <- setNames(rep(0, 8), algorithm_names)

for (i in 1:nrow(results_summary_full)) {
  scores <- as.numeric(results_summary_full[i, 2:9])
  current_ranks <- rank(scores, ties.method = "average")
  average_spot  <- average_spot + current_ranks
  print(current_ranks)
}

final_avg_rank <- average_spot / nrow(results_summary_full)

cat("\n--- Final Average Algorithm Ranks (1 = best) ---\n")
print(final_avg_rank)


