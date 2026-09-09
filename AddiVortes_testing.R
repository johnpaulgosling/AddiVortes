#devtools::install()
devtools::load_all()

#Friedman dataset (Works except Soft)
{

sim_fried <- function(N,P,sigma) {
  X <- matrix(runif(N * P), nrow = N, ncol = P)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)**2 + 10 * X[,4] + 5 * X[,5]
  Y <- mu + sigma * rnorm(N)
  
  return(list(X = X, Y = Y,mu=mu))
}

set.seed(74)
training_data <- sim_fried(500, 500, sqrt(10))
test_data <- sim_fried(500,500, sqrt(10))
X_train <- training_data$X
X_test <-test_data$X
Y_test_mu<-test_data$mu
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
  set.seed(74)
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
  training_data <- sim_CM2 (1500,sqrt(10))
  test_data <- sim_CM2 (1500, sqrt(10))
  X_train <- training_data$X
  X_test <- test_data$X
  Y_test_mu<-test_data$mu
}

#simulation case study (perofrms badly atm)
{
  library(MASS)
   simulated_case_study <- function(N,P) {

     Sigma <- diag(P)
     
     block_size <- 20
     num_blocks <- 5
     rho <- 0.95
     
     for (b in 1:num_blocks) {
       start_idx <- (b - 1) * block_size + 1
       end_idx <- b * block_size
       
       Sigma[start_idx:end_idx, start_idx:end_idx] <- rho
       diag(Sigma[start_idx:end_idx, start_idx:end_idx]) <- 1
     }
     
     X <- mvrnorm(n = N, mu = rep(0, P), Sigma = Sigma)
     
     x1  <- X[, 1]
     x21 <- X[, 21]
     x41 <- X[, 41]
     x61 <- X[, 61]
     
     f_x <- 10 * sin(pi * x1 * x21) + 20 * (x41 - 0.5)^2 + 10 * x61
     
     noise_sd <- 1.0
     Y <- f_x + rnorm(N, mean = 0, sd = noise_sd)
     
     simulated_data <- list(Y = Y, X=X,mu=f_x)
    return(simulated_data)
   }
   training_data <- simulated_case_study (300,2000)
   test_data <- simulated_case_study (300, 2000)
   X_train <- training_data$X
   X_test <- test_data$X
   Y_test_mu<-test_data$mu
}

#Smoothness Cases:
{
  # --- CASE 1: THE SMOOTHER (Low p / Forces gentle gradients) ---
  {
  p_init_val  <- 1.0
  p_shape_val <- 2.0
  p_rate_val  <- 5.0
  p_sd_val    <- 0.25
}
  # --- CASE 2: THE EXPLORER (Balanced / Data-driven adjustments) ---
  {
   p_init_val  <- 2.0
   p_shape_val <- 2.0
   p_rate_val  <- 0.5
   p_sd_val    <- 0.5
  }
  # --- CASE 3: THE SHARPENER (High p / Approximates hard Voronoi boundaries) ---
  { p_init_val  <- 5.0
   p_shape_val <- 2.0
   p_rate_val  <- 0.05
   p_sd_val    <- 1.0}
}

#Sparsity cases:
{
# --- CASE 1: EXTREME (Strict Sparsity / Ruthless Pruning) ---
{
alpha_val   <- 1
a_alpha_val <- 0.5
b_alpha_val <- 1
boost_val   <- 10.0
penalty_val <- 10.0
decay_val   <- 0.75
}
# --- CASE 2: DEFAULT (Balanced / Let the Data Decide) ---
{
 alpha_val   <- 1.0
 a_alpha_val <- 0.5
 b_alpha_val <- 1.0
 boost_val   <- 2.0
 penalty_val <- 10.0
 decay_val   <- 0.90
}
# --- CASE 3: CONSERVATIVE (Permissive / Retains Many Variables) ---
{ alpha_val   <- 5.0
 a_alpha_val <- 5.0
 b_alpha_val <- 0.5
 boost_val   <- 5.0
 penalty_val <- 0.1
 decay_val   <- 0.9
}
}

boost_val<-10
penalty_val<-10
alpha_val<-1

cat("Fitting local update algorithm...\n")
time_local <- system.time({
  set.seed(78)
  Model_AddiVortes_local<- AddiVortes(training_data$Y, X_train,m=200,thinning = 1,varSelMode = 2,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = TRUE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = 0.6,numChains = 1,IntialSigma = "LASSO",tau=50,splitMode =1)#,rho_alpha = 1)#,power = p_init_val, p_shape = p_shape_val, p_rate = p_rate_val, p_sd = p_sd_val)
})

time_local <- system.time({
  set.seed(78)
  Model_AddiVortes_original <- AddiVortes:::AddiVortes(training_data$Y, X_train,m=200,thinning = 1,varSelMode = 0,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = FALSE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = 0.6,numChains = 1,IntialSigma = "LASSO",tau=100,splitMode =1)#,rho_alpha = 1)#,power = p_init_val, p_shape = p_shape_val, p_rate = p_rate_val, p_sd = p_sd_val)
})

cat("Time taken for local:\n")
print(time_local)

##Graphs
{
  # Assuming you have run the algorithm and saved the output:
  # Model_AddiVortes_local <- Homo_AddiVortes_Algorithm(...)
  
  # --- 1. Extract the Dirichlet Weight metrics ---
  s_means <- Model_AddiVortes_local$posteriorDirichletWeightsMean
  s_lower <- Model_AddiVortes_local$posteriorDirichletWeightsLower
  s_upper <- Model_AddiVortes_local$posteriorDirichletWeightsUpper
  covariate_index <- 1:length(s_means)
  
  # --- 2. Plot the Posterior Dirichlet Weights ---
  plot(covariate_index, s_means, type = "p", pch = 19, col = "darkorange",
       ylim = c(0, max(s_upper) + 0.05),
       main = "Posterior Dirichlet Weights with 95% Credible Intervals",
       xlab = "Covariate Index", ylab = "Weight")
  
  # FIX: Create a logical mask to only draw arrows where the credible interval length > 0.
  # This completely prevents the "zero-length arrow" warning.
  valid_arrows <- (s_upper > s_lower)
  
  suppressWarnings(
    arrows(x0 = covariate_index, 
           y0 = s_lower, 
           x1 = covariate_index, 
           y1 = s_upper,
           angle = 90, code = 3, length = 0.05, col = "darkgray")
  )
  
    # --- 2. Plot the Posterior Dirichlet Weights ---
  plot(covariate_index, s_means, type = "p", pch = 19, col = "darkorange",
       ylim = c(0, max(s_upper) + 0.05),
       main = "Posterior Dirichlet Weights with 95% Credible Intervals",
       xlab = "Covariate Index", ylab = "Weight")
  
  # FIX: Create a logical mask to only draw arrows where the credible interval length > 0.
  # This completely prevents the "zero-length arrow" warning.
  valid_arrows <- (s_upper > s_lower)
  
  suppressWarnings(
    arrows(x0 = covariate_index, 
           y0 = s_lower, 
           x1 = covariate_index, 
           y1 = s_upper,
           angle = 90, code = 3, length = 0.05, col = "darkgray")
  )
  
  number_points<- 10
  ordered_indexes_mean_values<-order(s_means[-c(1:5)],decreasing = TRUE)[1:number_points]+5
  
  # 1. Define your data vectors
  x_vals <- c(1:5, covariate_index[ordered_indexes_mean_values])
  y_vals <- c(s_means[1:5], s_means[ordered_indexes_mean_values])
  y_lower <- c(s_lower[1:5], s_lower[ordered_indexes_mean_values])
  y_upper <- c(s_upper[1:5], s_upper[ordered_indexes_mean_values])
  
  # 2. Draw the bar chart and save the x-coordinates of the bars
  bp <- barplot(y_vals, 
                names.arg = x_vals, 
                col = "darkorange",
                ylim = c(0, max(s_upper) + 0.01),
                main = "Posterior Dirichlet Weights with 95% Credible Intervals",
                xlab = "Covariate Index", 
                ylab = "Weight")
  
  # 3. Add the error bars using the arrows() function
  arrows(x0 = bp, y0 = y_lower, x1 = bp, y1 = y_upper, 
         angle = 90, code = 3, length = 0.05, col = "black")
  
  # --- 3. Plot the Ensemble Inclusion Probabilities ---
  inclusion_probs <- Model_AddiVortes_local$ensembleInclusionProbabilities
  
  barplot(inclusion_probs, 
          main = "Ensemble Inclusion Probabilities", 
          xlab = "Covariate Index", 
          ylab = "Probability (0 to 1)",
          col = "steelblue",
          border = NA,
          ylim = c(0, 1),
          names.arg = covariate_index)

  bars_displayed<-10
  barplot(c(inclusion_probs[1:5],sort(inclusion_probs[-c(1:5)],decreasing = TRUE)[1:bars_displayed]), 
          main = "Ensemble Inclusion Probabilities", 
          xlab = "Covariate Index", 
          ylab = "Probability (0 to 1)",
          col = "steelblue",
          border = NA,
          ylim = c(0, 1),
          names.arg = 5+bars_displayed)
  
  # --- 4. MCMC Trace Plots ---
  plot(Model_AddiVortes_local$posteriorAlpha, type = "l", col = "purple",
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Alpha Value")
  
  plot(Model_AddiVortes_local,which=7)
  
  plot(Model_AddiVortes_local$posteriorSigma, type = "l", col = "blue", lwd = 2,
       main = "Posterior Sigma squared Trace Plot",
       xlab = "MCMC Iteration", ylab = "Sigma squared Value")
  
  # Replace '10' with the true sigma value if known for simulated data
  abline(h = 10, col = "red", lwd = 2, lty = 2) 
  
  acf(Model_AddiVortes_local$posteriorSigma, 
      main = expression(paste("ACF: Error Variance (", sigma^2, ")")),
      lwd = 2, col = "darkgreen", ylab = "Autocorrelation", xlab = "Lag")
}

prediction_AddiVortes<-predict(Model_AddiVortes_local,as.matrix(X_test))
rmse_AddiVortes<-sqrt(mean((prediction_AddiVortes-Y_test_mu)^2))
Accuracy_AddiVortes<-sum(round(pnorm(prediction_AddiVortes))==(as.numeric(test_data$Y)))/length(test_data$Y)
print(paste("RMSE for AddiVortes:", rmse_AddiVortes))
print(paste("Accuracy for AddiVortes:", Accuracy_AddiVortes))
sqrt(mean((mean(test_data$Y)-Y_test_mu)^2))
1-sum(as.numeric(test_data$Y))/length(test_data$Y)
plot(sort(pnorm(prediction_AddiVortes)),col = ifelse(test_data$Y[order(prediction_AddiVortes)]==1,'red','black'))
abline(h=0.5)

plot.AddiVortesFit(Model_AddiVortes_local_test,as.matrix(X_test),test_data$mu,which=c(1:3))

cat("Fitting local update algorithTRUEm...\n")
time_dart <- system.time({
  set.seed(123)
  Model_DART<-wbart(X_train,training_data$Y,X_test,ntree = 200,ndpost = 2500, nskip = 2500,sparse = TRUE)
})
print(time_dart)
plot(Model_DART$sigma[2500:5000]^2, type = "l", col = "blue", lwd = 2)
abline(h = 10, col = "red")
print(paste("RMSE for DART:", sqrt(mean((Model_DART$yhat.test.mean - Y_test_mu)^2))))
split_counts <- Model_DART$varcount
inclusion_proportions <- colMeans(split_counts / rowSums(split_counts))
barplot(inclusion_proportions, main = "dbarts Splitting Proportions (p > n)",ylab = "Posterior Splitting Proportion", xlab = "Covariate Index")
acf(Model_DART$sigma[2501:5000], 
    main = expression(paste("ACF: Error Variance (", sigma^2, ")")),
    lwd = 2, col = "darkgreen", ylab = "Autocorrelation", xlab = "Lag")

cat("Fitting local update algorithm...\n")
time_soft <- system.time({
  set.seed(123)
  Model_softBART<-softbart(X_train,training_data$Y,X_test,hypers = Hypers(X_train, training_data$Y,num_tree = 20))#,opts = Opts(num_burn = 200, num_save = 1000))
})
print(time_soft)
plot(Model_softBART$sigma^2, type = "l", col = "blue", lwd = 2)
abline(h = 10, col = "red")
print(Model_softBART$alpha)
barplot(colMeans(Model_softBART$s))
bars_displayed<-10
barplot(c(colMeans(Model_softBART$s)[1:5],sort(colMeans(Model_softBART$s)[-c(1:5)],decreasing = TRUE)[1:bars_displayed]), 
        main = "Ensemble Inclusion Probabilities", 
        xlab = "Covariate Index", 
        ylab = "Probability (0 to 1)",
        col = "darkorange",
        border = NA,
        names.arg = 5+bars_displayed)
print(paste("RMSE for SoftBART:", sqrt(mean((Model_softBART$y_hat_test_mean - Y_test_mu)^2))))
plot(Model_softBART$alpha)

task<- system.time({
  set.seed(123)
  Model_flexBART <- flexBART(Y~bart(.), train_data = training_data,M_vec=20,save_trees = TRUE,n.chains=1)#,family = binomial(link = "probit"))
  prediction_flexBART<-predict.flexBART(Model_flexBART, test_data)
})
print(task)
rmse_flexBART<-sqrt(mean((colMeans(prediction_flexBART)-test_data$Y)^2))
#rmse_flexBART<-sqrt(mean((Model_flexBART$yhat.test.mean-Y_test_mu)^2))
Accuracy_flexBART<-sum(round(colMeans(prediction_flexBART))==(as.numeric(test_data$Y)))/length(test_data$Y)
print(paste("Accuracy for flexBART:", Accuracy_flexBART))
print(paste("RMSE for flexBART:", rmse_flexBART))
plot(Model_flexBART$sigma^2, type = "l", col = "blue", lwd = 2)
abline(h = 10, col = "red")

plot(sort(colMeans(prediction_flexBART)),col = ifelse(test_data$Y[order(prediction_AddiVortes)]==1,'red','black'))
abline(h=0.5)

      