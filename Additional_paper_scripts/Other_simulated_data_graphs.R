#Friedman dataset
{
  
  sim_fried <- function(N,P,sigma) {
    X <- matrix(runif(N * P), nrow = N, ncol = P)
    mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)**2 + 10 * X[,4] + 5 * X[,5]
    Y <- mu + sigma * rnorm(N)
    
    return(list(X = X, Y = Y,mu=mu))
  }
  
  set.seed(45)
  training_data <- sim_fried(250, 200, sqrt(1))
  test_data <- sim_fried(250,200, sqrt(1))
  X_train <- training_data$X
  X_test <-test_data$X
  Y_test_mu<-test_data$mu
}

#Second Simulated study (works except Soft)
{
  sim_spam_continuous <- function(N, P, sigma) {
    X <- matrix(runif(N * P, min = -2.5, max = 2.5), nrow = N, ncol = P)
    
    f1 <- -sin(2 * X[, 1])
    f2 <- X[, 2]^2 - (25 / 12)
    f3 <- X[, 3]
    f4 <- exp(-X[, 4]) - (2 / 5) * sinh(2.5)
    
    mu <- f1 + f2 + f3 + f4
    Y <- mu + sigma * rnorm(N)
    
    return(data.frame(X = X, Y = Y,mu=mu))
  }
  set.seed(789)
  training_data <- sim_spam_continuous(1000, 500,sqrt(10))
  test_data <- sim_spam_continuous(1000,500, sqrt(10))
  X_train <- training_data[,1:500]
  X_test <- test_data[,1:500]
  Y_test_mu<-test_data$mu
}

#Binary Input (VerSelMode = 1 eventually converges if 500,000 iterations)
{
  sim_regime_binary <- function(N, P_cont, P_bin, sigma) {
    # Generate continuous predictors uniformly between -1 and 1
    X <- matrix(runif(N * P_cont, -1, 1), nrow = N, ncol = P_cont)
    
    # Generate binary predictors as fair coin flips
    B <- matrix(rbinom(N * P_bin, 1, 0.5), nrow = N, ncol = P_bin)
    
    # Initialize the mean response vector
    mu <- numeric(N)
    
    # Regime 1: B1 == 0 & B2 == 0
    idx1 <- B[, 1] == 0 & B[, 2] == 0
    mu[idx1] <- 15 * sin(pi * X[idx1, 1] * X[idx1, 2]) + 
      10 * (X[idx1, 3])^2 - 
      5 * X[idx1, 4]
    
    # Regime 2: B1 == 1 & B2 == 0
    idx2 <- B[, 1] == 1 & B[, 2] == 0
    mu[idx2] <- 15 * cos(pi * X[idx2, 1] * X[idx2, 2]) - 
      10 * (X[idx2, 3])^2 + 
      5 * X[idx2, 5]
    
    # Regime 3: B1 == 0 & B2 == 1
    idx3 <- B[, 1] == 0 & B[, 2] == 1
    mu[idx3] <- 10 * X[idx3, 1] + 
      10 * X[idx3, 2] + 
      20 * (X[idx3, 3])^2
    
    # Regime 4: B1 == 1 & B2 == 1
    idx4 <- B[, 1] == 1 & B[, 2] == 1
    mu[idx4] <- -10 * X[idx4, 1] - 
      10 * X[idx4, 2] - 
      20 * (X[idx4, 3])^2
    
    # Generate final response with Gaussian noise
    Y <- mu + sigma * rnorm(N)
    
    # Return data and true active indices for validation
    return(list(
      X = X, 
      B = B, 
      Y = Y,
      active_cont = 1:5,
      active_bin = 1:2,
      mu = mu
    ))
  }
  
  set.seed(74)
  training_data <- sim_regime_binary(500, 500, 50,sqrt(10))
  test_data <- sim_regime_binary(500,500,50, sqrt(10))
  X_train <- cbind(training_data$X,training_data$B)
  X_test <- cbind(test_data$X,test_data$B)
  Y_test_mu<-test_data$mu
}

#simulation study CM1 (works except Soft)
{
  sim_CM1 <- function(n, p, sigma) {
    # Calculate the number of binary and continuous predictors
    p_bin <- ceiling(p / 2)
    p_cont <- p - p_bin
    
    # Generate predictors
    X_bin <- matrix(rbinom(n * p_bin, size = 1, prob = 0.5), nrow = n, ncol = p_bin)
    X_cont <- matrix(runif(n * p_cont, min = 0, max = 1), nrow = n, ncol = p_cont)
    X <- cbind(X_bin, X_cont)
    
    # Calculate the true non-linear and non-additive signal
    f0 <- 10 * sin(pi * X[, p_bin + 1] * X[, p_bin + 2]) + 
      20 * (X[, p_bin + 3] - 0.5)^2 + 
      10 * X[, 1] + 
      5 * X[, 2]
    
    # Generate response with Gaussian noise
    y <- rnorm(n, mean = f0, sd = sigma)
    
    return(list(X = data.frame(X), Y = y,mu=f0))
  }
  set.seed(0394)
  training_data <- sim_CM1 (500, 500,sqrt(10))
  test_data <- sim_CM1 (500,500, sqrt(10))
  X_train <- training_data$X
  X_test <- test_data$X
  Y_test_mu<-test_data$mu
}

#simulation study CM2 (perofrms badly atm)
{
  sim_CM2 <- function(n, sigma) {
    # Generate binary predictors with varying probabilities
    X1_20 <- matrix(rbinom(n * 20, size = 1, prob = 0.2), nrow = n, ncol = 20)
    X21_40 <- matrix(rbinom(n * 20, size = 1, prob = 0.5), nrow = n, ncol = 20)
    
    # Generate continuous predictors with a block correlation of 0.3
    p_mvn <- 44 
    Sigma <- matrix(0.3, nrow = p_mvn, ncol = p_mvn)
    diag(Sigma) <- 1
    X41_84 <- MASS::mvrnorm(n, mu = rep(0, p_mvn), Sigma = Sigma)
    
    # Combine into the final design matrix
    X <- cbind(X1_20, X21_40, X41_84)
    
    # Calculate the true complex functional form
    f0 <- -4 + 
      2 * X[, 1] + 
      sin(X[, 12] * X[, 44]) * X[, 21] + 
      0.6 * X[, 41] * X[, 42] + 
      exp(-2 * (X[, 42] + 1)^2) - 
      X[, 43] + 
      0.5 * X[, 44]
    
    # Generate response with Gaussian noise
    y <- rnorm(n, mean = f0, sd = sigma)
    
    return(list(X = data.frame(X), Y = y,mu=f0))
  }
  
  set.seed(9834)
  training_data <- sim_CM2 (1500,sqrt(10))
  test_data <- sim_CM2 (1500, sqrt(10))
  X_train <- training_data$X
  X_test <- test_data$X
  Y_test_mu<-test_data$mu
}

#classification case (Binary classification is wrong in code)
{
  sim_fried_class <- function(N, P, sigma) {
    X <- matrix(runif(N * P), nrow = N, ncol = P)
    mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
    
    latent <- mu + sigma * rnorm(N)
    latent_scaled <- as.numeric(scale(latent))
    prob <- 1 / (1 + exp(-latent_scaled))
    
    Y <- rbinom(N, size = 1, prob = prob)
    
    return(data.frame(X = X, Y = Y))#, prob = prob, mu = mu))
  }
  
  set.seed(74)
  training_data <- sim_fried_class(500, 100, sqrt(1))
  test_data <- sim_fried_class(500,100, sqrt(1))
  X_train <- model.matrix(Y ~ . - 1, data = training_data)
  X_test <- model.matrix(Y ~ . - 1, data = test_data)
}

#Non-sparse huge dataset
{
  simulated_study_data<-function(n=500,p_active=500,p_inactive=500,sigma_true=1.0 ){
    set.seed(42)
    
    p <- p_active + p_inactive
    
    # 2. Generate independent uniform covariates [-1, 1]
    X <- matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p)
    colnames(X) <- paste0("V", 1:p)
    
    # 3. Generate the true Additive Response (Y) using only the first 500 variables
    #    - First 250 variables have a linear effect
    #    - Next 250 variables have a non-linear (quadratic/sine) effect
    linear_part <- rowSums(X[, 1:(p/2)])
    nonlinear_part <- rowSums(X[, round(p/2):round(3*p/4)]^2) + rowSums(sin(pi * X[, round(3*p/4):p]))
    
    Y_true <- linear_part + nonlinear_part
    
    # 4. Add Gaussian noise to create the final observation vector
    Y <- Y_true + rnorm(n, mean = 0, sd = sigma_true)
    
    return(list(X = X, Y = Y,mu=Y_true))
  }
  
  set.seed(74)
  training_data <- simulated_study_data(500, 100,50, sqrt(1))
  test_data <- simulated_study_data(500,100,50, sqrt(1))
  X_train <- training_data$X
  X_test <-test_data$X
  Y_test_mu<-test_data$mu
}

#correlated dataset
{
# Generates Correlated Data (Interaction-based, AR(1) dependence)
generate_dataset_2 <- function(n, p, sigma,rho) {
  Sigma <- outer(1:p, 1:p, function(j, k) rho^abs(j - k))
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:p)
  f_x <- 2 * X[, 1] * X[, 4] + 2 * X[, 7] * X[, 10]
  y <- f_x + rnorm(n, 0, sigma)
  return(list(X = X, Y = y,mu=f_x, true_vars = paste0("X", c(1, 4, 7, 10))))
}
  set.seed(74)
  training_data <- generate_dataset_2(250, 75, sqrt(1),0.9)
  test_data <- generate_dataset_2(250,75, sqrt(1),0.9)
  X_train <- training_data$X
  X_test <-test_data$X
  Y_test_mu<-test_data$mu
  
}
set.seed(186)

tau<-20
boost_val<-5
penalty_val<-5
decay_val<-0.6
alpha_val<-1
a_alpha_val <- 0.5
b_alpha_val <- 1
kappa_val <- 0.8

Model_AddiVortes_local <- AddiVortes(training_data$Y, X_train,m=200,Omega=1.5,thinning = 1,varSelMode = 2,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = TRUE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = kappa_val,numChains = 1,IntialSigma = "LASSO",tau=tau,splitMode =1)#,rho_alpha = 1)#,power = p_init_val, p_shape = p_shape_val, p_rate = p_rate_val, p_sd = p_sd_val)

Model_AddiVortes_original <- AddiVortes(training_data$Y, X_train,m=200,Omega=1.5,thinning = 1,varSelMode = 0,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = FALSE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = kappa_val,numChains = 1,IntialSigma = "LASSO",tau=20,splitMode =1)#,rho_alpha = 1)#,power = p_init_val, p_shape = p_shape_val, p_rate = p_rate_val, p_sd = p_sd_val)

Model_AddiVortes_dir <- AddiVortes(training_data$Y, X_train,m=200,Omega=1.5,thinning = 1,varSelMode = 1,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = TRUE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = kappa_val,numChains = 1,IntialSigma = "LASSO",tau=tau,splitMode =1)#,rho_alpha = 1)#,power = p_init_val, p_shape = p_shape_val, p_rate = p_rate_val, p_sd = p_sd_val)

Model_DART<-wbart(X_train,training_data$Y,X_test,ntree = 200,ndpost = 2500, nskip = 2500,sparse = TRUE)

#Model_DART<-wbart(X_train,training_data$Y,X_test,ntree = 200,ndpost = 2500, nskip = 2500,sparse = TRUE,augment =TRUE)

#Model_BART<-wbart(X_train,training_data$Y,X_test,ntree = 200,ndpost = 2500, nskip = 2500,sparse = FALSE)


col_adaptive <- "darkorange"
col_original <- "firebrick"
col_dart <- "forestgreen"
col_dir <- "darkorchid"
col_bart <- "royalblue"



### Figure 3 Augmented counts ###
{
  old_par <- par(no.readonly = TRUE)
  
  # Configure 1x2 layout with tight spacing and closer axis labels
  par(mfrow = c(1, 2), mgp = c(2.2, 0.4, 0))
  
  # Helper function to render horizontal barplots with 95% Credible Intervals
  plot_aug_counts <- function(model, bar_col) {
    aug_counts <- model$posteriorAugmentedCounts
    mean_counts <- rowMeans(aug_counts)
    ci_lower    <- apply(aug_counts, 1, quantile, probs = 0.025)
    ci_upper    <- apply(aug_counts, 1, quantile, probs = 0.975)
    
    cov_names <- if (!is.null(names(model$xCentres))) names(model$xCentres) else paste0("X", seq_along(mean_counts))
    ord <- order(mean_counts, decreasing = FALSE)
    if (length(ord) > 20) ord <- tail(ord, 20)
    
    # Dynamically adjust left margin for variable names; minimise top and right margins
    name_margin <- max(4.1, max(nchar(cov_names[ord])) * 0.55)
    par(mar = c(3.5, name_margin, 0.5, 0.5))
    
    bp <- barplot(mean_counts[ord], horiz = TRUE, names.arg = cov_names[ord], las = 1,
                  col = bar_col, border = NA, xlab = "Failure Count", cex.names = 0.8,
                  xlim = c(0, max(ci_upper[ord]) * 1.05))
    
    # 95% Credible Intervals with caps
    segments(x0 = ci_lower[ord], y0 = bp, x1 = ci_upper[ord], y1 = bp, lwd = 1.5)
    segments(x0 = ci_lower[ord], y0 = bp - 0.15, x1 = ci_lower[ord], y1 = bp + 0.15, lwd = 1.5)
    segments(x0 = ci_upper[ord], y0 = bp - 0.15, x1 = ci_upper[ord], y1 = bp + 0.15, lwd = 1.5)
  }
  
  # Render plots
  plot_aug_counts(Model_AddiVortes_dir, col_dir)
  plot_aug_counts(Model_AddiVortes_local, col_adaptive)
  
  # Reset plotting parameters
  par(old_par)
}

#### Figure 2 Sigma trace ####
{
  old_par <- par(no.readonly = TRUE)
  
  # Configure a 2x2 layout with zero inner margins (mar) and shared outer margins (oma)
  # This makes the plots perfectly flush and identical in size.
  par(mfrow = c(2, 2), mar = c(0.4, 0.4, 0.4, 0.4), oma = c(4, 5, 2, 1), mgp = c(2.5, 1, 0))
  
  # List of models and their respective data
  # Note: Extracting data to variables first keeps the plotting loop clean
  plot_data <- list(
    list(data = Model_AddiVortes_local$posteriorSigma, col = col_adaptive),
    list(data = Model_AddiVortes_dir$posteriorSigma,    col = col_dir),
    list(data = Model_DART$sigma[2500:5000]^2,         col = col_dart),
    list(data = Model_AddiVortes_original$posteriorSigma, col = col_original)
  )
  
  # Helper to draw each flush panel
  for (i in 1:4) {
    # Determine position (rows/cols)
    row <- if (i <= 2) 1 else 2
    col <- if (i %% 2 != 0) 1 else 2
    
    # Plot without axes or labels
    plot(plot_data[[i]]$data, type = "l", col = plot_data[[i]]$col, 
         lwd = 1.5, ylim = c(2, 19), axes = FALSE, xlab = "", ylab = "", main = "")
    
    # Add horizontal reference line
    abline(h = 10, col = "black", lwd = 1.5)
    
    # Draw border box
    box()
    
    # Add y-axis ticks only to the first column (plots 1 and 3)
    if (col == 1) axis(2, las = 1, xpd = NA)
    
    # Add x-axis ticks only to the bottom row (plots 3 and 4)
    if (row == 2) axis(1, xpd = NA)
  }
  
  # Add shared outer labels
  mtext("MCMC Iteration", side = 1, outer = TRUE, line = 2.5)
  mtext(expression(sigma^2 ~ "sample"), side = 2, outer = TRUE, line = 2)
  
  # Reset plotting parameters
  par(old_par)
}

### Figure 4 TRUE vs Predicted ###
{
  # 1. Compute posterior predictive means and 95% intervals for each model
  pred_addivortes_local    <- predict(Model_AddiVortes_local, as.matrix(X_test), showProgress = FALSE)
  ci_addivortes_local      <- predict(Model_AddiVortes_local, as.matrix(X_test), type = "quantile", quantiles = c(0.025, 0.975), showProgress = FALSE)
  
  pred_addivortes_original <- predict(Model_AddiVortes_original, as.matrix(X_test), showProgress = FALSE)
  ci_addivortes_original   <- predict(Model_AddiVortes_original, as.matrix(X_test), type = "quantile", quantiles = c(0.025, 0.975), showProgress = FALSE)
  
  pred_addivortes_dir      <- predict(Model_AddiVortes_dir, as.matrix(X_test), showProgress = FALSE)
  ci_addivortes_dir        <- predict(Model_AddiVortes_dir, as.matrix(X_test), type = "quantile", quantiles = c(0.025, 0.975), showProgress = FALSE)
  
  pred_dart <- colMeans(Model_DART$yhat.test)
  ci_dart   <- t(apply(Model_DART$yhat.test, 2, quantile, probs = c(0.025, 0.975)))
  
  # 2. Configure a tight 2x2 plotting grid
  # oma sets outer margins to hold the shared axes and labels, allowing inner margins to be tiny
  par(mfrow = c(2, 2), oma = c(4, 4, 1, 1), mar = c(0.5, 0.5, 0.5, 0.5), mgp = c(2.5, 1, 0))
  
  # Determine common axis limits using both point predictions and interval boundaries
  all_values <- c(Y_test_mu, ci_addivortes_local, ci_addivortes_original, ci_addivortes_dir, ci_dart)
  plot_limits <- range(all_values, na.rm = TRUE)
  
  # Helper function updated to manage axes manually and use light grey bars
  plot_with_ci <- function(x_vals, y_vals, ci_matrix, pt_col, draw_x_axis, draw_y_axis) {
    # Set up an empty plot with correct limits, removing internal axes completely
    plot(x_vals, y_vals, type = "n", main = "", xlab = "", ylab = "", 
         xlim = plot_limits, ylim = plot_limits, xaxt = "n", yaxt = "n")
    
    # Draw light grey confidence intervals first so they sit behind the points
    segments(x0 = x_vals, y0 = ci_matrix[, 1], x1 = x_vals, y1 = ci_matrix[, 2], 
             col = "lightgrey", lwd = 2.5)
    
    # Add the solid point predictions and the 1:1 reference line
    points(x_vals, y_vals, col = pt_col, pch = 16, cex = 0.8)
    abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = 2)
    
    # Conditionally draw the axes in the outer margin space
    if(draw_x_axis) axis(1)
    if(draw_y_axis) axis(2)
  }
  
  # 3. Plot 1: Adaptive AddiVortes (Top-Left)
  plot_with_ci(Y_test_mu, pred_addivortes_local, ci_addivortes_local, 
               col_adaptive, draw_x_axis = FALSE, draw_y_axis = TRUE)
  
  # 4. Plot 2: AddiVortes Dirichlet (Top-Right)
  plot_with_ci(Y_test_mu, pred_addivortes_dir, ci_addivortes_dir, 
               col_dir, draw_x_axis = FALSE, draw_y_axis = FALSE)
  
  # 5. Plot 3: DART (Bottom-Left)
  plot_with_ci(Y_test_mu, pred_dart, ci_dart, 
               col_dart, draw_x_axis = TRUE, draw_y_axis = TRUE)
  
  # 6. Plot 4: Original AddiVortes (Bottom-Right)
  plot_with_ci(Y_test_mu, pred_addivortes_original, ci_addivortes_original, 
               col_original, draw_x_axis = TRUE, draw_y_axis = FALSE)
  
  # Add shared axis labels using the outer margin
  mtext("True values", side = 1, outer = TRUE, line = 2.5)
  mtext("Predicted values", side = 2, outer = TRUE, line = 2.5)
  
  # Reset plotting parameters to default, including outer margins
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0), oma = c(0, 0, 0, 0))
  
  # Output the RMSE scores
  print(paste("Adaptive: ", sqrt(mean((Y_test_mu - pred_addivortes_local)^2))))
  print(paste("Uniform: ", sqrt(mean((Y_test_mu - pred_addivortes_original)^2))))
  print(paste("Dirichlet: ", sqrt(mean((Y_test_mu - pred_addivortes_dir)^2))))
  print(paste("DART: ", sqrt(mean((Y_test_mu - pred_dart)^2))))
  print(paste("Benchmark: ", sqrt(mean((Y_test_mu - mean(Y_test_mu))^2))))
}

### Figure 9 Number of Active Dimensions ###
{
  active_adapt <- colSums(Model_AddiVortes_local$posteriorVariableSelection > 0)
  active_orig  <- colSums(Model_AddiVortes_original$posteriorVariableSelection > 0)
  active_dir <- colSums(Model_AddiVortes_dir$posteriorVariableSelection > 0)
  
  # DART: varcount is (ndpost x p)
  active_dart  <- rowSums(Model_DART$varcount > 0)
  active_bart  <- rowSums(Model_BART$varcount > 0)
  
  # Plot
  plot(active_adapt, type = "l", col = col_adaptive, lwd = 2,
       ylim = c(0, max(active_adapt, active_orig, active_dart)),
       xlab = "MCMC Iteration", ylab = "Number of Active Covariates")
  lines(active_orig, col = col_original, lwd = 2)
  lines(active_dart, col = col_dart, lwd = 2)
  lines(active_dir, col = col_dir, lwd = 2)
  lines(active_bart, col = col_bart, lwd = 2)
  
  
  # Add ground truth line (5 active variables)
  abline(h = 5, col = "black", lwd = 2, lty = 2)
  
  plot(active_adapt, type = "l", col = col_adaptive, lwd = 2,
       ylim = c(0, max(active_adapt, active_dart)),
       xlab = "MCMC Iteration", ylab = "Number of Active Covariates")
  lines(active_dart, col = col_dart, lwd = 2)
  lines(active_dir, col = col_dir, lwd = 2)
  
  # Add ground truth line (5 active variables)
  abline(h = 5, col = "black", lwd = 2, lty = 2)
  legend("topright", legend = c("Adaptive AddiVortes","Dirichlet AddiVortes", "DART" ,"True Sparsity (5)"),
         col = c(col_adaptive, col_dir, col_dart, "black"), 
         lwd = c(2, 2,2, 2), lty = c(1, 1,1,2), bty = "n", cex = 0.8)
}

### Figure 6 Posterior Dirichlet Weights ###
{
  old_par <- par(no.readonly = TRUE)
  
  # Define independent text scaling sizes
  label_size <- 1.4
  y_axis_size <- 1.2
  x_axis_size <- 0.7
  
  # 1. Prepare DART matrix early to calculate global y-axis limits
  split_counts <- Model_DART$varcount
  proportions_matrix <- split_counts / rowSums(split_counts)
  dart_upper <- apply(proportions_matrix, 2, quantile, probs = 0.975)
  
  # 2. Calculate the global maximum upper bound for a unified y-axis
  max_local <- max(Model_AddiVortes_local$posteriorDirichletWeightsUpper, na.rm = TRUE)
  max_dir   <- max(Model_AddiVortes_dir$posteriorDirichletWeightsUpper, na.rm = TRUE)
  max_dart  <- max(dart_upper, na.rm = TRUE)
  global_y_max <- max(0.8, max_local, max_dir, max_dart) * 1.05
  
  # Configure 1x3 layout. 
  # 'oma' creates outer margins to hold the shared x-axis label at the bottom and y-axis on the left.
  par(mfrow = c(1, 3), mgp = c(2.5, 1, 0), oma = c(4, 5, 1, 1), 
      cex.lab = label_size)
  
  # Helper function to extract, sort, and plot the inclusion proportions
  plot_inclusion <- function(means, lower, upper, bar_col, draw_y_axis) {
    
    number_points <- 6
    ordered_rest <- order(means[-c(1:5)], decreasing = TRUE)[1:number_points] + 5
    
    y_vals <- c(means[1:5], means[ordered_rest])
    y_lower <- c(lower[1:5], lower[ordered_rest])
    y_upper <- c(upper[1:5], upper[ordered_rest])
    x_names <- c(1:5, ordered_rest)
    
    # Draw the bar chart WITHOUT internal axes or labels to allow global control
    bp <- barplot(y_vals, 
                  col = bar_col,
                  main = "",
                  xlab = "", 
                  names.arg = rep("", length(y_vals)), 
                  ylab = "", 
                  ylim = c(0, global_y_max),
                  axes = FALSE)
    
    # Force horizontal x-axis labels to render at every bar location using the smaller x_axis_size
    mtext(text = x_names, side = 1, at = bp, line = 0.8, cex = x_axis_size, las = 1)
    
    # Conditionally draw the y-axis in the global outer margin space using the larger y_axis_size
    if(draw_y_axis) {
      axis(2, outer = TRUE, las = 1, cex.axis = y_axis_size)
    }
    
    # Add the error bars
    arrows(x0 = bp, y0 = y_lower, x1 = bp, y1 = y_upper, 
           angle = 90, code = 3, length = 0.05, col = "black")
    
    # Add a vertical dashed red line perfectly between the 5th and 6th bars
    abline(v = (bp[5] + bp[6]) / 2, col = "red", lty = 2)
  }
  
  # Set identical inner margins for ALL plots to guarantee exact same width.
  # bottom = 3, left = 0 (y-axis pushed outside), top = 1, right = 0.5 (for tight spacing).
  uniform_mar <- c(3, 0, 1, 0.5)
  
  # --- Plot 1: Adaptive AddiVortes (Left) ---
  par(mar = uniform_mar) 
  plot_inclusion(Model_AddiVortes_local$posteriorDirichletWeightsMean,
                 Model_AddiVortes_local$posteriorDirichletWeightsLower,
                 Model_AddiVortes_local$posteriorDirichletWeightsUpper,
                 col_adaptive, draw_y_axis = TRUE)
  
  # Add the y-axis label to the outer margin
  mtext(expression(s[j]), side = 2, outer = TRUE, line = 3, cex = label_size)
  
  # --- Plot 2: Dirichlet AddiVortes (Middle) ---
  par(mar = uniform_mar) 
  plot_inclusion(Model_AddiVortes_dir$posteriorDirichletWeightsMean,
                 Model_AddiVortes_dir$posteriorDirichletWeightsLower,
                 Model_AddiVortes_dir$posteriorDirichletWeightsUpper,
                 col_dir, draw_y_axis = FALSE)
  
  # --- Plot 3: DART (Right) ---
  par(mar = uniform_mar) 
  plot_inclusion(colMeans(proportions_matrix),
                 apply(proportions_matrix, 2, quantile, probs = 0.025),
                 dart_upper,
                 col_dart, draw_y_axis = FALSE)
  
  # Add the single overarching x-axis label at the very bottom
  mtext("Covariate Index", side = 1, outer = TRUE, line = 2.5, cex = label_size)
  
  # Reset plotting parameters
  par(old_par)
}

### Figure 7 Ensemble inclusion probabilities ###
{
  old_par <- par(no.readonly = TRUE)
  
  # Define text scaling sizes to make the font bigger
  label_size <- 1.4
  axis_size <- 1.2
  
  # Configure a tight 2x2 layout
  # oma (outer margins) holds the shared axes and labels.
  par(mfrow = c(2, 2), oma = c(4, 5, 1, 1), mgp = c(2.5, 1, 0))
  
  # Set identical inner margins for ALL plots to guarantee exact same width.
  # bottom = 2.5 (room for x-values), left = 0 (y-axis pushed outside), top = 0.5, right = 1 (small spacing).
  par(mar = c(2.5, 0, 0.5, 1))
  
  # Helper function to extract, sort, and plot the ensemble inclusion probabilities
  plot_ensemble_inclusion <- function(probs, bar_col, draw_y_axis) {
    bars_displayed <- 10
    
    # Identify the original indices of the top noise variables (covariates 6 onwards)
    noise_indices <- 6:length(probs)
    noise_probs <- probs[noise_indices]
    
    # Get the indices of the highest noise variables and map back to original indices
    top_noise_order <- order(noise_probs, decreasing = TRUE)[1:bars_displayed]
    top_noise_actual_indices <- noise_indices[top_noise_order]
    
    # Combine the values and names for the plot
    plot_values <- c(probs[1:5], probs[top_noise_actual_indices])
    plot_names <- c(1:5, top_noise_actual_indices)
    
    # Generate the barplot. 
    # names.arg ensures the x-values print on ALL graphs.
    # axes = FALSE ensures NO solid lines are drawn for the axes.
    bp <- barplot(plot_values, 
                  names.arg = plot_names,
                  cex.names = axis_size,
                  main = "", xlab = "", ylab = "", 
                  col = bar_col, border = NA, ylim = c(0, 1),
                  axes = FALSE)
    
    # Conditionally draw the y-axis in the global outer margin space for the left plots
    if(draw_y_axis) {
      axis(2, outer = TRUE, las = 1, cex.axis = axis_size)
    }
  }
  
  # --- Plot 1: Adaptive AddiVortes (Top-Left) ---
  plot_ensemble_inclusion(Model_AddiVortes_local$ensembleInclusionProbabilities, 
                          col_adaptive, draw_y_axis = TRUE)
  
  # --- Plot 2: Dirichlet AddiVortes (Top-Right) ---
  plot_ensemble_inclusion(Model_AddiVortes_dir$ensembleInclusionProbabilities, 
                          col_dir, draw_y_axis = FALSE)
  
  # --- Plot 3: DART (Bottom-Left) ---
  # Calculate probability of inclusion in the ensemble (count > 0)
  prob_in_ensemble_dart <- colMeans(Model_DART$varcount > 0)
  
  plot_ensemble_inclusion(prob_in_ensemble_dart, 
                          col_dart, draw_y_axis = TRUE)
  
  # --- Plot 4: Original AddiVortes (Bottom-Right) ---
  plot_ensemble_inclusion(Model_AddiVortes_original$ensembleInclusionProbabilities, 
                          col_original, draw_y_axis = FALSE)
  
  # Add shared axis labels using the outer margin with bigger font
  mtext("Covariate Index", side = 1, outer = TRUE, line = 2, cex = label_size)
  mtext("Proportion of inclusion", side = 2, outer = TRUE, line = 3.5, cex = label_size)
  
  # Reset plotting parameters
  par(old_par)
}

### Figure 10 trace of variable count ###
{
  par(mfrow=c(2,5))
  for (j in c(1:8,378,254)) {
    # Pre-allocating the vector size speeds up the loop significantly
    number_of_inclusions <- numeric(2500) 
    for(i in 1:2500) {
      number_of_inclusions[i] <- sum(unlist(Model_AddiVortes_local$posteriorDim[[i]]) == j)
    }
    
    # Conditionals for axis limits and labels
    y_lim <- if (j <= 5) c(0, 200) else c(0, 10)
    y_lab <- if (j == 1 || j == 6) "Number of inclusions" else ""
    
    plot(number_of_inclusions, type = 'l', col = col_adaptive, lwd = 2, 
         ylim = y_lim, xlab = "MCMC iteration", ylab = y_lab)
  }
  
  # --- Dirichlet AddiVortes ---
  for (j in 1:10) {
    number_of_inclusions <- numeric(2500)
    for(i in 1:2500) {
      number_of_inclusions[i] <- sum(unlist(Model_AddiVortes_dir$posteriorDim[[i]]) == j)
    }
    
    y_lim <- if (j <= 5) c(0, 200) else c(0, 10)
    y_lab <- if (j == 1 || j == 6) "Number of inclusions" else ""
    
    plot(number_of_inclusions, type = 'l', col = col_dir, lwd = 2, 
         ylim = y_lim, xlab = "MCMC iteration", ylab = y_lab)
  }
  
  # --- DART Model ---
  for (j in 1:10) {
    y_lim <- if (j <= 5) c(0, 200) else c(0, 10)
    y_lab <- if (j == 1 || j == 6) "Number of inclusions" else ""
    
    plot(Model_DART$varcount[, j], type = 'l', col = col_dart, lwd = 2, 
         ylim = y_lim, xlab = "MCMC iteration", ylab = y_lab)
  }
  par(mfrow=c(1,1))
  
  {
    old_par <- par(no.readonly = TRUE)
    
    # 1. Add a small inner margin (mar) for gaps, keeping oma for shared outer axes
    par(mfrow = c(3, 5), mar = c(0.4, 0.4, 0.4, 0.4), oma = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
    
    var_indices <- 1:5
    n_iter <- 2500
    
    # 2. Pre-allocate and extract inclusion counts for speed
    inc_adaptive <- matrix(0, nrow = n_iter, ncol = 5)
    inc_dir      <- matrix(0, nrow = n_iter, ncol = 5)
    
    for (i in 1:n_iter) {
      for (j in var_indices) {
        inc_adaptive[i, j] <- sum(unlist(Model_AddiVortes_local$posteriorDim[[i]]) == j)
        inc_dir[i, j]      <- sum(unlist(Model_AddiVortes_dir$posteriorDim[[i]]) == j)
      }
    }
    
    # Helper to draw each panel with a small gap
    plot_panel <- function(data, bar_col, row, col) {
      # Draw the plot without default axes or labels
      plot(data, type = 'l', col = bar_col, lwd = 2, ylim = c(0, 200),
           axes = FALSE, main = "", xlab = "", ylab = "")
      
      # Draw a clean border box around each plot
      box()
      
      # Add y-axis ticks only to the first column
      if (col == 1) {
        axis(2, las = 1, xpd = NA)
      }
      
      # Add x-axis ticks only to the bottom row
      if (row == 3) {
        axis(1, xpd = NA)
      }
    }
    
    # --- Row 1: Adaptive AddiVortes ---
    for (j in var_indices) {
      plot_panel(inc_adaptive[, j], col_adaptive, row = 1, col = j)
    }
    
    # --- Row 2: Dirichlet AddiVortes ---
    for (j in var_indices) {
      plot_panel(inc_dir[, j], col_dir, row = 2, col = j)
    }
    
    # --- Row 3: DART Model ---
    for (j in var_indices) {
      plot_panel(Model_DART$varcount[, j], col_dart, row = 3, col = j)
    }
    
    # 3. Add the shared outer labels
    mtext("MCMC iteration", side = 1, outer = TRUE, line = 2.5)
    mtext("Number of inclusions", side = 2, outer = TRUE, line = 2.5)
    
    # Reset plotting parameters
    par(old_par)
  }
}

### Figure 9 Alpha trace ###
{
  old_par <- par(no.readonly = TRUE)
  
  # Configure a tight 1x2 layout
  par(mfrow = c(1, 2), mgp = c(2.5, 1, 0))
  
  # --- Left Plot: Adaptive Alpha Trace ---
  # Standard left/bottom margins, minimised top and right margins
  par(mar = c(4, 4, 0.5, 0.5))
  plot(1:length(Model_AddiVortes_local$posteriorAlpha), Model_AddiVortes_local$posteriorAlpha, 
       type = "l", xlab = "MCMC Iteration", ylab = expression(alpha ~ "sample"), main = "",
       col = col_adaptive, lwd = 1.5, ylim = c(0, 6.5))
  abline(h = mean(Model_AddiVortes_local$posteriorAlpha, na.rm = TRUE), col = "red", lty = 2)
  
  # --- Right Plot: Dirichlet Alpha Trace ---
  # Minimized left margin to pull it close, no y-axis label
  par(mar = c(4, 1.5, 0.5, 1))
  plot(1:length(Model_AddiVortes_dir$posteriorAlpha), Model_AddiVortes_dir$posteriorAlpha, 
       type = "l", xlab = "MCMC Iteration", ylab = "", main = "",
       col = col_dir, lwd = 1.5, ylim = c(0, 6.5))
  abline(h = mean(Model_AddiVortes_dir$posteriorAlpha, na.rm = TRUE), col = "red", lty = 2)
  
  # Reset plotting parameters
  par(old_par)
  
}

### Figure 7 Probability of included active variable ### 
{
  old_par <- par(no.readonly = TRUE)
  
  # Define a common size for all axis labels so they match perfectly
  common_label_size <- 1.1
  
  # Configure a tight 1x3 layout.
  # oma (outer margins): bottom, left, top, right. 
  # Increased the left margin from 4 to 6 to prevent the y-axis label from being cut off.
  par(mfrow = c(1, 3), oma = c(3, 6, 0.5, 0.5), mgp = c(2.5, 1, 0), cex.lab = common_label_size)
  
  num_covariates <- length(X_test[1, ])
  
  # Set identical inner margins for ALL plots to guarantee exact same width.
  # (bottom = 2, left = 0, top = 0.5, right = 1). 
  # The right margin of 1 creates a small, consistent spacing between the graphs.
  par(mar = c(2, 0, 0.5, 1))
  
  # --- Left Plot: Adaptive AddiVortes Weights ---
  sum_active_local <- colSums(Model_AddiVortes_local$posteriorDirichletWeights[1:5, ])
  sum_noise_local <- colSums(Model_AddiVortes_local$posteriorDirichletWeights[6:num_covariates, ])
  
  # yaxt = "n" prevents standard axis generation inside the plot
  plot(sum_active_local, type = "l", col = col_adaptive,
       xlab = "", ylab = "", yaxt = "n", 
       main = "", ylim = c(0, 1))
  lines(sum_noise_local, col = 'red')
  abline(h = mean(sum_active_local), lty = 2)
  abline(h = mean(sum_noise_local), col = 'red', lty = 2)
  
  # Draw the y-axis values in the left outer margin so they don't compress the plot width
  axis(2, outer = TRUE)
  
  # Draw the y-axis label in the left outer margin, pushed slightly further out (line = 3)
  mtext(expression("Sum of active and inactive "~ "s"[j]), side = 2, outer = TRUE, line = 3, cex = common_label_size)
  
  # --- Middle Plot: Dirichlet AddiVortes Weights ---
  sum_active_dir <- colSums(Model_AddiVortes_dir$posteriorDirichletWeights[1:5, ])
  sum_noise_dir <- colSums(Model_AddiVortes_dir$posteriorDirichletWeights[6:num_covariates, ])
  
  plot(sum_active_dir, type = "l", col = col_dir,
       xlab = "", ylab = "", yaxt = "n", 
       main = "", ylim = c(0, 1))
  lines(sum_noise_dir, col = 'red')
  abline(h = mean(sum_active_dir), lty = 2)
  abline(h = mean(sum_noise_dir), col = 'red', lty = 2)
  
  # --- Right Plot: DART Weights ---
  sum_active_dart <- rowSums(Model_DART$varprob[, 1:5])
  sum_noise_dart <- rowSums(Model_DART$varprob[, 6:num_covariates])
  
  plot(sum_active_dart, type = "l", col = col_dart, 
       xlab = "", ylab = "", yaxt = "n", 
       main = "", ylim = c(0, 1))
  lines(sum_noise_dart, col = 'red')
  abline(h = mean(sum_active_dart), lty = 2)
  abline(h = mean(sum_noise_dart), col = 'red', lty = 2)
  
  # Add the single shared x-axis label in the bottom outer margin
  mtext("MCMC Iteration", side = 1, outer = TRUE, line = 1, cex = common_label_size)
  
  # Reset plotting parameters
  par(old_par)
}

### Figure 11 Momentum barchart ###
{
  old_par <- par(no.readonly = TRUE)
  
  # 1. Calculate the mean, lower (2.5%), and upper (97.5%) bounds
  mean_momentum <- rowMeans(Model_AddiVortes_local$posteriorMomentum)
  lower_momentum <- apply(Model_AddiVortes_local$posteriorMomentum, 1, quantile, probs = 0.025)
  upper_momentum <- apply(Model_AddiVortes_local$posteriorMomentum, 1, quantile, probs = 0.975)
  
  cov_names <- if(!is.null(names(Model_AddiVortes_local$xCentres))) names(Model_AddiVortes_local$xCentres) else paste0("X", seq_along(mean_momentum))
  
  # 2. Sort based on the mean
  ord <- order(mean_momentum, decreasing = FALSE)
  
  max_vars_to_plot <- 20
  if (length(mean_momentum) > max_vars_to_plot) {
    ord <- tail(ord, max_vars_to_plot)
  }
  
  # Dynamically adjust left margin for variable names; minimise top and right margins
  name_margin <- max(4.1, max(nchar(cov_names[ord])) * 0.6)
  par(mar = c(4, name_margin, 0.5, 1) + 0.1)
  
  # 3. Calculate dynamic xlim to ensure error bars aren't cut off
  x_min <- min(0, min(lower_momentum[ord])) 
  x_max <- max(upper_momentum[ord])
  
  # 4. Save the bar midpoints into 'bp' (Title removed)
  bp <- barplot(mean_momentum[ord], horiz = TRUE, names.arg = cov_names[ord], las = 1,
                col = col_adaptive, border = NA, xlab = "Mean Momentum",
                main = "", cex.names = 0.8,
                xlim = c(x_min, x_max * 1.05))
  
  # 5. Draw the main horizontal error bar lines
  segments(x0 = lower_momentum[ord], y0 = bp, 
           x1 = upper_momentum[ord], y1 = bp, 
           col = "black", lwd = 1.5)
  
  # 6. Draw the vertical caps for the error bars
  epsilon <- 0.15 # Height of the cap
  segments(x0 = lower_momentum[ord], y0 = bp - epsilon, 
           x1 = lower_momentum[ord], y1 = bp + epsilon, 
           col = "black", lwd = 1.5)
  segments(x0 = upper_momentum[ord], y0 = bp - epsilon, 
           x1 = upper_momentum[ord], y1 = bp + epsilon, 
           col = "black", lwd = 1.5)
  
  # Reset plotting parameters
  par(old_par)
}

### Figure 12 Momentum traces ###
{
  old_par <- par(no.readonly = TRUE)
  
  # Configure a 2x5 grid with small inner gaps and shared outer margins
  par(mfrow = c(2, 5), mar = c(0.4, 0.4, 0.4, 0.4), oma = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
  
  trace_indices <- c(1, 2, 3, 4, 5, 8, 168, 362, 341, 78)
  
  shared_ylim <- range(Model_AddiVortes_local$posteriorMomentum[trace_indices, ])
  
  # Helper function to draw each trace panel cleanly
  plot_trace_panel <- function(idx, row, col) {
    is_top_row <- row == 1
    show_y <- col == 1
    show_x <- row == 2
    
    # Draw the plot without default axes or labels
    plot(Model_AddiVortes_local$posteriorMomentum[idx, ], 
         type = 'l', col = col_adaptive, lwd = 1.5, ylim = shared_ylim,
         axes = FALSE, main = "", xlab = "", ylab = "")
    
    # Draw a clean border box around each plot
    box()
    
    # Add y-axis ticks only to the first column, allowing them to bleed outside
    if (show_y) {
      axis(2, las = 1, xpd = NA)
    }
    
    # Add x-axis ticks only to the bottom row, allowing them to bleed outside
    if (show_x) {
      axis(1, xpd = NA)
    }
    
    # Add the dashed red line for the burn-in threshold
    abline(v = 1500, col = "red", lty = 2)
  }
  
  # Plot all 10 traces across the 2x5 grid
  for (i in seq_along(trace_indices)) {
    idx <- trace_indices[i]
    current_row <- if (i <= 5) 1 else 2
    current_col <- if (i %% 5 == 0) 5 else i %% 5
    
    plot_trace_panel(idx, row = current_row, col = current_col)
  }
  
  # Add the shared outer labels
  mtext("MCMC Iteration", side = 1, outer = TRUE, line = 2.5)
  mtext("Posterior Momentum", side = 2, outer = TRUE, line = 2.5)
  
  # Reset plotting parameters
  par(old_par)
}

### Figure RMSE 8 values comparison ### Time ~ 4 minutes
{
  library(parallel)
  library(pbapply)
  library(BART)

  
  # Detect available cores and leave one free for system stability
  n_cores <- 10 #max(1, detectCores() - 3)
  cl <- makeCluster(n_cores)
  
  # Export all current variables and functions from your global environment to the workers
  #clusterExport(cl, ls())
  
  library(doRNG)
  
  # Ensure the BART package is loaded on all worker nodes
  clusterEvalQ(cl, { 
    library(BART)
    devtools::load_all() 
  })
  registerDoRNG(seed = 539)
  
  n_iterations <- 10
  model_names <- c("Original", "Dirichlet", "Adaptive", "DART")
  
  cat("Running parallel simulations across", n_cores, "cores...\n")
  
  # Run the parallelised loop with a progress bar
  results_list <- pblapply(1:n_iterations, function(i) {
    
    sim_fried <- function(N,P,sigma) {
      X <- matrix(runif(N * P), nrow = N, ncol = P)
      mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)**2 + 10 * X[,4] + 5 * X[,5]
      Y <- mu + sigma * rnorm(N)
      
      return(list(X = X, Y = Y,mu=mu))
    }
    
    # 1. Generate new data for this specific iteration
    iter_train <- sim_fried(200, 400, sqrt(10))
    iter_test <- sim_fried(200, 400, sqrt(10))
    
    tau<-20
    boost_val<-3
    penalty_val<-3
    decay_val<-0.9
    alpha_val<-1
    a_alpha_val <- 0.5
    b_alpha_val <- 1
    kappa_val <- 0.8
    
    # 2. Fit the models using your predefined hyperparameters
    fit_orig <- AddiVortes(iter_train$Y, iter_train$X, m=200,Omega=1.5,thinning = 1,varSelMode = 0,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = TRUE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = kappa_val,numChains = 1,IntialSigma = "LASSO",tau=tau,splitMode =1, showProgress=FALSE)
    
    fit_dir <- AddiVortes(iter_train$Y, iter_train$X, m=200,Omega=1.5,thinning = 1,varSelMode = 1,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = TRUE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = kappa_val,numChains = 1,IntialSigma = "LASSO",tau=tau,splitMode =1, showProgress=FALSE)
    
    fit_adapt <- AddiVortes(iter_train$Y, iter_train$X, m=200,Omega=1.5,thinning = 1,varSelMode = 2,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = TRUE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = kappa_val,numChains = 1,IntialSigma = "LASSO",tau=tau,splitMode =1, showProgress=FALSE)
    
    fit_dart <- wbart(iter_train$X, iter_train$Y, iter_test$X, ntree=200, 
                      ndpost=2500, nskip=2500, sparse=TRUE)
    
    # 3. Calculate point predictions (posterior means) for the test set
    pred_orig <- predict(fit_orig, as.matrix(iter_test$X), showProgress=FALSE)
    pred_dir <- predict(fit_dir, as.matrix(iter_test$X), showProgress=FALSE)
    pred_adapt <- predict(fit_adapt, as.matrix(iter_test$X), showProgress=FALSE)
    pred_dart <- colMeans(fit_dart$yhat.test)
    
    # 4. Calculate and return the RMSE against the true test mu
    c(
      Original = sqrt(mean((pred_orig - iter_test$mu)^2)),
      Dirichlet = sqrt(mean((pred_dir - iter_test$mu)^2)),
      Adaptive = sqrt(mean((pred_adapt - iter_test$mu)^2)),
      DART = sqrt(mean((pred_dart - iter_test$mu)^2))
    )
  }, cl = cl)
  
  stopCluster(cl)
  
  # Combine the individual list results into the final matrix
  rmse_results <- do.call(rbind, results_list)
  rownames(rmse_results) <- NULL
  
  
  old_par <- par(no.readonly = TRUE)
  
  # Configure a tight 1x1 layout with a minimised top margin
  par(mfrow = c(1, 1), mar = c(4, 4, 1, 1) + 0.1)
  
  # 1. Align the custom colours with the column order of rmse_results
  simulation_colours <- c(col_original, col_dir, col_adaptive, col_dart)
  
  # 2. Generate the boxplot for the independent runs with no title
  boxplot(rmse_results,
          main = "",
          xlab = "RMSE",
          col = simulation_colours,
          pch = 16,
          outcol = rgb(0, 0, 0, alpha = 0.3),
          horizontal = TRUE
  )
  
  # Reset plotting parameters
  par(old_par)
  
}


### Additional Graphs ###

### Average #Centres/#Dimensions per tessellation
{
  average_dimensions_adaptive <- sapply(Model_AddiVortes_local$posteriorDim, function(tessellations) {
    length(unlist(tessellations))/200
  })
  average_dimensions_original <- sapply(Model_AddiVortes_original$posteriorDim, function(tessellations) {
    length(unlist(tessellations))/200
  })
  average_dimensions_dir <- sapply(Model_AddiVortes_dir$posteriorDim, function(tessellations) {
    length(unlist(tessellations))/200
  })
  
  average_centres_adaptive <- sapply(Model_AddiVortes_local$posteriorPred, function(tessellations) {
    length(unlist(tessellations))/200
  })
  average_centres_dir <- sapply(Model_AddiVortes_dir$posteriorPred, function(tessellations) {
    length(unlist(tessellations))/200
  })
  average_centres_original <- sapply(Model_AddiVortes_original$posteriorPred, function(tessellations) {
    length(unlist(tessellations))/200
  })
  
  par(mfrow=c(1,2))
  plot(average_dimensions_adaptive, type = "l", col = col_adaptive, lwd = 2,
       xlab = "MCMC Iteration", ylab = "Average number of dimensions",ylim=c(2,5))
  lines(average_dimensions_dir,col= col_dir)
  lines(average_dimensions_original,col= col_original)
  
  plot(average_centres_adaptive, type = "l", col = col_adaptive, lwd = 2,
       xlab = "MCMC Iteration", ylab = "Average number of centres",ylim=c(4,6))
  lines(average_centres_dir,col= col_dir)
  lines(average_centres_original,col= col_original)
  legend("bottomright", legend = c("Adaptive AddiVortes","Dirichlet AddiVortes", "Original AddiVortes"),
         col = c(col_adaptive, col_dir,col_original, col_dart,col_bart, "black"), 
         lwd = c(2, 2, 2), lty = c(1, 1, 1), bty = "n", cex = 0.8)
  
  par(mfrow=c(1,1))
}

### Trace of augment counts ##
{
  {
    old_par <- par(no.readonly = TRUE)
    
    # Configure a 2x5 grid with small inner gaps and shared outer margins
    par(mfrow = c(2, 5), mar = c(0.4, 0.4, 0.4, 0.4), oma = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
    
    var_indices <- 1:5
    
    # Calculate a shared y-axis limit for a fair comparison across all ten traces
    max_count <- max(
      Model_AddiVortes_local$posteriorAugmentedCounts[var_indices, ],
      Model_AddiVortes_dir$posteriorAugmentedCounts[var_indices, ],
      na.rm = TRUE
    )
    shared_ylim <- c(0, max_count * 1.05)
    
    # Helper function to draw each trace panel cleanly
    plot_trace_panel <- function(data, trace_col, row, col) {
      # Draw the plot without default axes or labels
      plot(data, type = 'l', col = trace_col, lwd = 1.5, ylim = shared_ylim,
           axes = FALSE, main = "", xlab = "", ylab = "")
      
      # Draw a clean border box around each plot
      box()
      
      # Add y-axis ticks only to the first column, allowing them to bleed outside
      if (col == 1) {
        axis(2, las = 1, xpd = NA)
      }
      
      # Add x-axis ticks only to the bottom row, allowing them to bleed outside
      if (row == 2) {
        axis(1, xpd = NA)
      }
    }
    
    # --- Row 1: Adaptive AddiVortes ---
    for (j in var_indices) {
      plot_trace_panel(Model_AddiVortes_local$posteriorAugmentedCounts[j, ], col_adaptive, row = 1, col = j)
    }
    
    # --- Row 2: Dirichlet AddiVortes ---
    for (j in var_indices) {
      plot_trace_panel(Model_AddiVortes_dir$posteriorAugmentedCounts[j, ], col_dir, row = 2, col = j)
    }
    
    # Add the shared outer labels
    mtext("MCMC iteration", side = 1, outer = TRUE, line = 2.5)
    mtext("Augmented Counts", side = 2, outer = TRUE, line = 2.5)
    
    var_indices <- 6:10
    
    # Calculate a shared y-axis limit for covariates 6-10
    max_count <- max(
      Model_AddiVortes_local$posteriorAugmentedCounts[var_indices, ],
      Model_AddiVortes_dir$posteriorAugmentedCounts[var_indices, ],
      na.rm = TRUE
    )
    shared_ylim <- c(0, max_count * 1.05)
    
    # Helper function to draw each trace panel cleanly
    plot_trace_panel <- function(data, trace_col, row, col) {
      # Draw the plot without default axes or labels
      plot(data, type = 'l', col = trace_col, lwd = 1.5, ylim = shared_ylim,
           axes = FALSE, main = "", xlab = "", ylab = "")
      
      # Draw a clean border box around each plot
      box()
      
      # Add y-axis ticks only to the first column, allowing them to bleed outside
      if (col == 1) {
        axis(2, las = 1, xpd = NA)
      }
      
      # Add x-axis ticks only to the bottom row, allowing them to bleed outside
      if (row == 2) {
        axis(1, xpd = NA)
      }
    }
    
    # --- Row 1: Adaptive AddiVortes ---
    for (idx in seq_along(var_indices)) {
      j <- var_indices[idx]
      plot_trace_panel(Model_AddiVortes_local$posteriorAugmentedCounts[j, ], col_adaptive, row = 1, col = idx)
    }
    
    # --- Row 2: Dirichlet AddiVortes ---
    for (idx in seq_along(var_indices)) {
      j <- var_indices[idx]
      plot_trace_panel(Model_AddiVortes_dir$posteriorAugmentedCounts[j, ], col_dir, row = 2, col = idx)
    }
    
    # Add the shared outer labels
    mtext("MCMC iteration", side = 1, outer = TRUE, line = 2.5)
    mtext("Augmented Counts", side = 2, outer = TRUE, line = 2.5)
    
    # Reset plotting parameters
    par(old_par)
  }
}

### Graph 6.5 Trace of s values ###
{
  par(mfrow=c(2,5))
  plot(Model_AddiVortes_local$posteriorDirichletWeights[1,], type = 'l', col = col_adaptive, lwd = 2, xlab = "MCMC iteration")
  plot(Model_AddiVortes_local$posteriorDirichletWeights[2,], type = 'l', col = col_adaptive, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_local$posteriorDirichletWeights[3,], type = 'l', col = col_adaptive, lwd = 2, xlab = "MCMC iteration")
  plot(Model_AddiVortes_local$posteriorDirichletWeights[4,], type = 'l', col = col_adaptive, lwd = 2, xlab = "MCMC iteration")
  plot(Model_AddiVortes_local$posteriorDirichletWeights[5,], type = 'l', col = col_adaptive, lwd = 2, xlab = "MCMC iteration")
  
  plot(Model_AddiVortes_local$posteriorDirichletWeights[10,], type = 'l', col = col_adaptive, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_local$posteriorDirichletWeights[150,], type = 'l', col = col_adaptive, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_local$posteriorDirichletWeights[400,], type = 'l', col = col_adaptive, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_local$posteriorDirichletWeights[327,], type = 'l', col = col_adaptive, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_local$posteriorDirichletWeights[27,], type = 'l', col = col_adaptive, lwd = 2, xlab = "MCMC iteration",ylab="")
  
  
  plot(Model_AddiVortes_dir$posteriorDirichletWeights[1,], type = 'l', col = col_dir, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_dir$posteriorDirichletWeights[2,], type = 'l', col = col_dir, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_dir$posteriorDirichletWeights[3,], type = 'l', col = col_dir, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_dir$posteriorDirichletWeights[4,], type = 'l', col = col_dir, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_dir$posteriorDirichletWeights[5,], type = 'l', col = col_dir, lwd = 2, xlab = "MCMC iteration",ylab="")
  
  plot(Model_AddiVortes_dir$posteriorDirichletWeights[10,], type = 'l',col = col_dir, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_dir$posteriorDirichletWeights[150,], type = 'l',col = col_dir, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_dir$posteriorDirichletWeights[400,], type = 'l',col = col_dir, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_dir$posteriorDirichletWeights[327,], type = 'l',col = col_dir, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_AddiVortes_dir$posteriorDirichletWeights[27,], type = 'l',col = col_dir, lwd = 2, xlab = "MCMC iteration",ylab="")
  
  
  plot(Model_DART$varprob[,1], type = 'l', col = col_dart, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_DART$varprob[,2], type = 'l', col = col_dart, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_DART$varprob[,3], type = 'l', col = col_dart, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_DART$varprob[,4], type = 'l', col = col_dart, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_DART$varprob[,5], type = 'l', col = col_dart, lwd = 2, xlab = "MCMC iteration",ylab="")
  
  plot(Model_DART$varprob[,10], type = 'l', col = col_dart, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_DART$varprob[,150], type = 'l', col = col_dart, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_DART$varprob[,400], type = 'l', col = col_dart, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_DART$varprob[,327], type = 'l', col = col_dart, lwd = 2, xlab = "MCMC iteration",ylab="")
  plot(Model_DART$varprob[,27], type = 'l', col = col_dart, lwd = 2, xlab = "MCMC iteration",ylab="")
  
  par(mfrow=c(1,1))
}



##F1_score
{# Define the total number of covariates and the true active set
  p <- 75
  true_active <- c(1, 4, 7, 10)
  
  # Create the ground truth binary vector
  y_true <- rep(0, p)
  y_true[true_active] <- 1
  
  # Define a function to calculate F1 for a single iteration's prediction vector
  calc_f1_vector <- function(y_pred) {
    tp <- sum(y_true == 1 & y_pred == 1)
    fp <- sum(y_true == 0 & y_pred == 1)
    fn <- sum(y_true == 1 & y_pred == 0)
    
    precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall    <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    
    f1 <- ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
    return(f1)
  }
  
  # Define a tiny tolerance to account for floating-point noise in continuous weights
  tol <- 1e-8
  
  # Binarise AddiVortes weights (transposing so iterations are rows, matching DART)
  inc_local <- t(Model_AddiVortes_local$posteriorVariableSelection) > tol
  inc_orig  <- t(Model_AddiVortes_original$posteriorVariableSelection) > tol
  inc_dir   <- t(Model_AddiVortes_dir$posteriorVariableSelection) > tol
  
  # Binarise DART split counts (iterations are already rows)
  inc_dart <- Model_DART$varcount > 0
  
  # Calculate the mean F1 score across all posterior iterations
  f1_local <- mean(apply(inc_local, 1, calc_f1_vector))
  f1_orig  <- mean(apply(inc_orig, 1, calc_f1_vector))
  f1_dir   <- mean(apply(inc_dir, 1, calc_f1_vector))
  f1_dart  <- mean(apply(inc_dart, 1, calc_f1_vector))
  
  # Compile and print the final F1 scores
  f1_results <- data.frame(
    Model = c("AV-Adaptive", "AV-Original", "AV-Dirichlet", "DART"),
    Mean_Posterior_F1 = c(f1_local, f1_orig, f1_dir, f1_dart)
  )
  
  print(f1_results)}








### Graph 8 RMSE values comparison ###
{
  n_iterations <- 5
  model_names <- c("Original", "Dirichlet", "Adaptive", "DART")
  rmse_results <- matrix(NA, nrow = n_iterations, ncol = length(model_names))
  colnames(rmse_results) <- model_names
  
  # Initialise a text progress bar to monitor the loop
  pb <- txtProgressBar(min = 0, max = n_iterations, style = 3)
  
  for(i in 1:n_iterations) {
    
    # 1. Generate new data for this specific iteration
    iter_train <- sim_fried(500, 500, sqrt(10))
    iter_test <- sim_fried(500, 500, sqrt(10))
    
    # 2. Fit the models using your predefined hyperparameters
    fit_orig <- AddiVortes(iter_train$Y, iter_train$X, m=200, thinning=1, varSelMode=0, 
                           totalMCMCIter=5000, mcmcBurnIn=2500, dirichletWarmup=1000, nu=6, q=0.9, 
                           updateAlpha=FALSE, alpha=alpha_val, a_alpha=a_alpha_val, b_alpha=b_alpha_val, 
                           adaptBoost=boost_val, adaptPenalty=penalty_val, momentumDecay=decay_val, 
                           kappa=0.6, numChains=1, IntialSigma="LASSO", tau=100, splitMode=1)
    
    fit_dir <- AddiVortes(iter_train$Y, iter_train$X, m=200, thinning=1, varSelMode=1, 
                          totalMCMCIter=5000, mcmcBurnIn=2500, dirichletWarmup=1000, nu=6, q=0.9, 
                          updateAlpha=FALSE, alpha=alpha_val, a_alpha=a_alpha_val, b_alpha=b_alpha_val, 
                          adaptBoost=boost_val, adaptPenalty=penalty_val, momentumDecay=decay_val, 
                          kappa=0.6, numChains=1, IntialSigma="LASSO", tau=100, splitMode=1)
    
    fit_adapt <- AddiVortes(iter_train$Y, iter_train$X, m=200, thinning=1, varSelMode=2, 
                            totalMCMCIter=5000, mcmcBurnIn=2500, dirichletWarmup=1000, nu=6, q=0.9, 
                            updateAlpha=FALSE, alpha=alpha_val, a_alpha=a_alpha_val, b_alpha=b_alpha_val, 
                            adaptBoost=boost_val, adaptPenalty=penalty_val, momentumDecay=decay_val, 
                            kappa=0.6, numChains=1, IntialSigma="LASSO", tau=100, splitMode=1)
    
    fit_dart <- wbart(iter_train$X, iter_train$Y, iter_test$X, ntree=200, 
                      ndpost=2500, nskip=2500, sparse=TRUE)
    
    # 3. Calculate point predictions (posterior means) for the test set
    pred_orig <- predict(fit_orig,as.matrix(iter_test$X))
    pred_dir <- predict(fit_dir,as.matrix(iter_test$X))
    pred_adapt <- predict(fit_adapt,as.matrix(iter_test$X))
    pred_dart <- colMeans(fit_dart$yhat.test)
    
    # 4. Calculate and store the RMSE against the true test mu
    rmse_results[i, "Original"] <- sqrt(mean((pred_orig - iter_test$mu)^2))
    rmse_results[i, "Dirichlet"] <- sqrt(mean((pred_dir - iter_test$mu)^2))
    rmse_results[i, "Adaptive"] <- sqrt(mean((pred_adapt - iter_test$mu)^2))
    rmse_results[i, "DART"] <- sqrt(mean((pred_dart - iter_test$mu)^2))
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  { par(mfrow=c(1,1))
    # 1. Align the custom colours with the column order of rmse_results
    # ("Original", "Dirichlet", "Adaptive", "DART")
    simulation_colours <- c(col_original, col_dir, col_adaptive, col_dart)
    
    # 2. Generate the boxplot for the 100 independent runs
    boxplot(rmse_results,
            main = "Out-of-Sample RMSE Across 100 Simulations (Friedman Dataset)",
            xlab = "RMSE",
            col = simulation_colours,
            pch = 16,
            outcol = rgb(0, 0, 0, alpha = 0.3),
            horizontal=TRUE
    ) 
  }
}

### Graph 9 Hyperparmeter robustness ###
{
  # Define the grids for each hyperparameter to be tested
  param_grids <- list(
    adaptBoost = c(0.5, 2.0,10,50),
    adaptPenalty = c(0.5, 2.0, 10,50),
    tau = c(10, 100, 1000),
    m = c(50, 100, 200),
    alpha = c(0.5, 1.0, 2.0, 5.0)
  )
  
  n_reps <- 2
  
  # Set up a 2x3 plotting layout to accommodate the five graphs
  par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))
  
  # Loop through each hyperparameter
  for(param_name in names(param_grids)) {
    grid_values <- param_grids[[param_name]]
    robustness_results <- matrix(NA, nrow = n_reps, ncol = length(grid_values))
    colnames(robustness_results) <- paste0(param_name, " = ", grid_values)
    
    cat("\nTesting robustness for:", param_name, "\n")
    pb <- txtProgressBar(min = 0, max = length(grid_values) * n_reps, style = 3)
    counter <- 0
    
    ## Loop through each hyperparameter
    for(param_name in names(param_grids)) {
      grid_values <- param_grids[[param_name]]
      robustness_results <- matrix(NA, nrow = n_reps, ncol = length(grid_values))
      colnames(robustness_results) <- paste0(param_name, " = ", grid_values)
      
      cat("\nTesting robustness for:", param_name, "\n")
      pb <- txtProgressBar(min = 0, max = length(grid_values) * n_reps, style = 3)
      counter <- 0
      
      # Loop through the grid values and repetitions
      for(j in seq_along(grid_values)) {
        for(i in 1:n_reps) {
          # Generate fresh training and test sets
          iter_train <- sim_fried(100, 100, sqrt(10))
          iter_test <- sim_fried(100, 100, sqrt(10))
          
          # Define the base arguments for the model using lowercase y and x
          model_args <- list(
            y = iter_train$Y, x = iter_train$X, m = 200, thinning = 1, varSelMode = 2, 
            totalMCMCIter = 5000, mcmcBurnIn = 2500, dirichletWarmup = 1000, nu = 6, q = 0.9, 
            updateAlpha = FALSE, alpha = alpha_val, a_alpha = a_alpha_val, b_alpha = b_alpha_val, 
            adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, 
            kappa = 0.6, numChains = 1, IntialSigma = "LASSO", tau = 100, splitMode = 1
          )
          
          # Dynamically overwrite the specific hyperparameter being tested
          model_args[[param_name]] <- grid_values[j]
          
          # Fit the model using do.call to pass the updated argument list
          fit_adapt <- do.call(AddiVortes, model_args)
          
          # Calculate predictions and evaluate against true mu
          pred_adapt <- predict(fit_adapt,iter_test$X)
          robustness_results[i, j] <- sqrt(mean((pred_adapt - iter_test$mu)^2))
          
          counter <- counter + 1
          setTxtProgressBar(pb, counter)
        }
      }
      
      close(pb)
      
      # Generate the boxplot for the current parameter
      boxplot(robustness_results,
              main = paste("Robustness:", param_name),
              ylab = "Out-of-Sample RMSE",
              col = col_adaptive,
              pch = 16,
              outcol = rgb(0, 0, 0, alpha = 0.3),
              las = 1,
              cex.axis = 0.8)
    }
    
    # Reset the graphical parameters to default
    par(mfrow = c(1, 1))
  }
}

### Acceptence probability ## 


t<-1:4000
kappa<-0.9
tau<-c(1,10,50,100,1000)
col<-c('red','blue','green','purple','orange')

plot(((1+tau[1])/(t+tau[1]))^kappa,type ='l',col=col[1])
for(i in 2:length(tau)){
  values<-((1+tau[i])/(t+tau[i]))^kappa
  
  lines(values,col=col[i])
}

abline(v=1500)

t<-1:4000
kappa<-c(0.3,0.51,0.8,0.9,0.99)
col<-c('red','blue','green','purple','orange')

plot(((1)/(t))^kappa[1],type ='l',col=col[1],ylim = c(0,1))
for(i in 2:length(kappa)){
  values<-((1)/(t))^kappa[i]
  
  lines(values,col=col[i])
}

abline(v=1500)