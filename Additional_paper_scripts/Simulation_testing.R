# ==============================================================================
# AddiVortes Comprehensive Simulation Study
# ==============================================================================
library(MASS)
library(dbarts)
library(randomForest)
library(glmnet)

# 1. Data Generation Function
# ------------------------------------------------------------------------------
generate_sim_data <- function(n, p, dgm = "friedman", rho = 0.5, snr = 3.0) {
  
  # Construct a Toeplitz auto-regressive covariance matrix for covariates
  Sigma <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] <- rho^abs(i - j)
    }
  }
  
  # Generate covariates from Multivariate Normal, then transform to Uniform(0,1)
  X_raw <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X <- pnorm(X_raw)
  
  # Initialize true target surface
  f_true <- rep(0, n)
  active_vars <- integer(0)
  
  # Define Data Generating Mechanisms (DGMs)
  if (dgm == "linear") {
    # Sparse Linear (Favours LASSO)
    beta <- c(3, 2, 1.5, 1, 0.5, rep(0, p - 5))
    f_true <- X %*% beta
    active_vars <- 1:5
    
  } else if (dgm == "friedman") {
    # Friedman 1 (Standard Non-Linear Benchmark)
    f_true <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
    active_vars <- 1:5
    
  } else if (dgm == "oblique") {
    # Oblique Decision Boundary (Favours Voronoi over axis-aligned trees)
    f_true <- ifelse(X[,1] + X[,2] > 1.0, 5, -5) + ifelse(X[,3] - X[,4] > 0, 3, -3)
    active_vars <- 1:4
    
  } else if (dgm == "checkerboard") {
    # Complex Spatial Discontinuity 
    f_true <- 5 * sign(sin(2 * pi * X[,1]) * sin(2 * pi * X[,2])) * X[,3]
    active_vars <- 1:3
  }
  
  # Scale noise to achieve exact Signal-to-Noise Ratio (SNR)
  var_f <- var(f_true)
  sigma_noise <- sqrt(var_f / snr)
  y <- f_true + rnorm(n, 0, sigma_noise)
  
  return(list(X = X, y = y, f_true = f_true, active_vars = active_vars, sigma = sigma_noise))
}

# 2. Simulation Execution Loop
# ------------------------------------------------------------------------------
run_simulation <- function(n_train = 500, n_test = 500, p = 100, 
                           dgm_list = c("linear", "friedman", "oblique", "checkerboard"),
                           rho_list = c(0.0, 0.8), snr_list = c(1.0, 5.0), n_reps = 10) {
  
  results <- data.frame()
  
  for (dgm in dgm_list) {
    for (rho in rho_list) {
      for (snr in snr_list) {
        for (rep in 1:n_reps) {
          
          cat(sprintf("Running: DGM=%s, Rho=%.1f, SNR=%.1f, Rep=%d\n", dgm, rho, snr, rep))
          
          # Generate Train and Test sets
          train_data <- generate_sim_data(n_train, p, dgm, rho, snr)
          test_data  <- generate_sim_data(n_test, p, dgm, rho, snr)
          
          X_tr <- train_data$X; y_tr <- train_data$y
          X_te <- test_data$X; y_te <- test_data$y
          
          # --- Model 1: AddiVortes ---
          # Note: Tune hyperparameters to match your paper's exact specifications
          fit_av <- AddiVortes(y = y_tr, x = X_tr, m = 200, totalMCMCIter = 3000, 
                               mcmcBurnIn = 1500, varSelMode = 2, showProgress = FALSE,numChains=1,IntialSigma = "LASSO")
          # Assuming you have a predict_AddiVortes function:
           pred_av <- predict(fit_av, X_te)
           rmse_av <- sqrt(mean((y_te - pred_av)^2))
          #rmse_av <- fit_av$inSampleRmse # Placeholder if predict is unavailable
          
          # Variable Selection Metrics (AddiVortes)
          inc_probs <- fit_av$ensembleInclusionProbabilities
          selected_av <- which(inc_probs > 0.5)
          tpr_av <- sum(selected_av %in% train_data$active_vars) / length(train_data$active_vars)
          fpr_av <- sum(!(selected_av %in% train_data$active_vars)) / (p - length(train_data$active_vars))
          
          # --- Model 2: BART (dbarts) ---
          fit_bart <- bart(X_tr, y_tr, x.test = X_te, ndpost = 1500, nskip = 1500, verbose = FALSE,sigest)
          rmse_bart <- sqrt(mean((y_te - fit_bart$yhat.test.mean)^2))
          
          # --- Model 3: Random Forest ---
          fit_rf <- randomForest(x = X_tr, y = y_tr, xtest = X_te)
          rmse_rf <- sqrt(mean((y_te - fit_rf$test$predicted)^2))
          
          # --- Model 4: LASSO (glmnet) ---
          cv_lasso <- cv.glmnet(X_tr, y_tr, alpha = 1)
          pred_lasso <- predict(cv_lasso, newx = X_te, s = "lambda.min")
          rmse_lasso <- sqrt(mean((y_te - pred_lasso)^2))
          
          # Store Results
          res_row <- data.frame(
            DGM = dgm, Rho = rho, SNR = snr, Rep = rep,
            RMSE_AV = rmse_av, RMSE_BART = rmse_bart, RMSE_RF = rmse_rf, RMSE_LASSO = rmse_lasso,
            TPR_AV = tpr_av, FPR_AV = fpr_av
          )
          results <- rbind(results, res_row)
        }
      }
    }
  }
  return(results)
}

# Execute (Warning: This will take significant time. Run in parallel on a cluster for the actual paper).
sim_results <- run_simulation(n_train = 100,n_test=100, p = 10, n_reps = 2)

colSum()
# write.csv(sim_results, "AddiVortes_Sim_Results.csv", row.names = FALSE)


# ==============================================================================
# AddiVortes Data Generating Mechanism (Regression & Classification)
# ==============================================================================
library(MASS)

generate_sim_data <- function(n = 500, p = 100, 
                              dgm = "friedman_reg", 
                              rho = 0.5, snr = 3.0, 
                              prop_binary_features = 0.2) {
  
  # 1. Construct Auto-regressive Covariance Matrix
  Sigma <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] <- rho^abs(i - j)
    }
  }
  
  # 2. Generate Covariates (Gaussian Copula to Uniform)
  X_raw <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X <- pnorm(X_raw) 
  
  # 3. Convert a proportion of features to Binary Classifiers
  # We select the last (p * prop_binary_features) columns
  n_bin <- max(1, floor(p * prop_binary_features))
  bin_cols <- (p - n_bin + 1):p
  X[, bin_cols] <- ifelse(X[, bin_cols] > 0.5, 1, 0)
  
  # Extract indices for easy assignment in DGMs
  c1 <- 1; c2 <- 2; c3 <- 3; c4 <- 4; c5 <- 5 # Continuous indices
  b1 <- bin_cols[1]; b2 <- bin_cols[2]        # Binary indices
  
  f_true <- rep(0, n)
  active_vars <- integer(0)
  is_classification <- FALSE
  
  # ============================================================================
  # CONTINUOUS REGRESSION DGMs
  # ============================================================================
  if (dgm == "linear_reg") {
    # Sparse Linear with mixed continuous and binary effects
    f_true <- 3 * X[,c1] + 2 * X[,c2] - 1.5 * X[,c3] + 2.5 * X[,b1] - 2 * X[,b2]
    active_vars <- c(c1, c2, c3, b1, b2)
    
  } else if (dgm == "friedman_reg") {
    # Standard Friedman augmented with a binary shift
    f_true <- 10 * sin(pi * X[,c1] * X[,c2]) + 20 * (X[,c3] - 0.5)^2 + 
      10 * X[,c4] + 5 * X[,c5] + 5 * X[,b1]
    active_vars <- c(c1, c2, c3, c4, c5, b1)
  }
  
  # ============================================================================
  # BINARY CLASSIFICATION DGMs (Using Probit Link)
  # ============================================================================
  else if (dgm == "linear_class") {
    # Linear decision boundary with mixed features
    f_true <- 1.5 * X[,c1] - 1.5 * X[,c2] + 2.0 * X[,b1] - 2.0 * X[,b2]
    active_vars <- c(c1, c2, b1, b2)
    is_classification <- TRUE
    
  } else if (dgm == "rings_class") {
    # Non-linear Concentric Rings (Highly favourable to Voronoi topologies)
    # The binary feature acts as a spatial shift
    radius <- sqrt((X[,c1] - 0.5)^2 + (X[,c2] - 0.5)^2)
    f_true <- 5 - 15 * radius + 2.0 * X[,b1]
    active_vars <- c(c1, c2, b1)
    is_classification <- TRUE
  }
  
  # ============================================================================
  # RESPONSE GENERATION
  # ============================================================================
  if (is_classification) {
    # For classification, f_true acts as the latent continuous variable Z
    # We pass it through the standard normal CDF (Probit link)
    probs <- pnorm(f_true)
    y <- rbinom(n, size = 1, prob = probs)
    sigma_noise <- NA # Not applicable for categorical outcomes
    
  } else {
    # For regression, we scale the noise to match the exact desired SNR
    var_f <- var(f_true)
    sigma_noise <- sqrt(var_f / snr)
    y <- f_true + rnorm(n, 0, sigma_noise)
  }
  
  return(list(
    X = X, 
    y = y, 
    f_true = f_true, 
    active_vars = active_vars, 
    sigma = sigma_noise,
    is_classification = is_classification
  ))
}
