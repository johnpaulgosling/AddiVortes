
devtools::load_all()

true_value<- function(X){
  mu<-10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)**2 + 10 * X[,4] + 5 * X[,5]
  return(mu)
}

sim_fried <- function(N,P,sigma) {
  X <- matrix(runif(N * P), nrow = N, ncol = P)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)**2 + 10 * X[,4] + 5 * X[,5]
  Y <- mu + sigma * rnorm(N)
  
  return(data.frame(X = X, Y = Y))
}

sim_fried_class <- function(N, P, sigma) {
  X <- matrix(runif(N * P), nrow = N, ncol = P)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
  
  latent <- mu + sigma * rnorm(N)
  latent_scaled <- as.numeric(scale(latent))
  prob <- 1 / (1 + exp(-latent_scaled))
  
  Y <- rbinom(N, size = 1, prob = prob)
  
  return(data.frame(X = X, Y = Y))#, prob = prob, mu = mu))
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
alpha_val   <- 0.1
a_alpha_val <- 0.1
b_alpha_val <- 10.0
boost_val   <- 1.0
penalty_val <- 8.0
decay_val   <- 0.50
}
# --- CASE 2: DEFAULT (Balanced / Let the Data Decide) ---
{
 alpha_val   <- 1.0
 a_alpha_val <- 0.5
 b_alpha_val <- 1.0
 boost_val   <- 2.0
 penalty_val <- 2.0
 decay_val   <- 0.90
}
# --- CASE 3: CONSERVATIVE (Permissive / Retains Many Variables) ---
{ alpha_val   <- 5.0
 a_alpha_val <- 5.0
 b_alpha_val <- 0.5
 boost_val   <- 5.0
 penalty_val <- 0.1
 decay_val   <- 0.99
}
}

set.seed(74)
training_data <- sim_fried(100, 100, sqrt(1))
test_data <- sim_fried(100,100, sqrt(1))
X_train <- model.matrix(Y ~ . - 1 , data = training_data)
X_test <- model.matrix(Y ~ . - 1, data = test_data)
Y_test_mu<-true_value(X_test)

set.seed(74)
training_data <- sim_fried_class(100, 100, sqrt(10))
test_data <- sim_fried_class(100,100, sqrt(10))
X_train <- model.matrix(Y ~ . - 1, data = training_data)
X_test <- model.matrix(Y ~ . - 1, data = test_data)

cat("Fitting local update algorithm...\n")
time_local <- system.time({
  Model_AddiVortes_local <- AddiVortes(training_data$Y, X_train,m=200,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 500,nu=6,q=0.9,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = 0.6,numChains = 1,IntialSigma = "LASSO",varSelMode = 1,splitMode =1)#,power = p_init_val, p_shape = p_shape_val, p_rate = p_rate_val, p_sd = p_sd_val)
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
  
  # --- 4. MCMC Trace Plots ---
  plot(Model_AddiVortes_local$posteriorAlpha, type = "l", col = "purple",
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Alpha Value")
  
  plot(Model_AddiVortes_local$posteriorSigma, type = "l", col = "blue", lwd = 2,
       main = "Posterior Sigma squared Trace Plot",
       xlab = "MCMC Iteration", ylab = "Sigma Value")
  
  # Replace '10' with the true sigma value if known for simulated data
  abline(h = 10, col = "red", lwd = 2, lty = 2) 
}

prediction_AddiVortes<-predict(Model_AddiVortes_local,X_test)
rmse_AddiVortes<-sqrt(mean((prediction_AddiVortes-Y_test_mu)^2))
Accuracy_AddiVortes<-sum(round(pnorm(prediction_AddiVortes))==(as.numeric(test_data$Y)))/length(test_data$Y)
print(paste("RMSE for AddiVortes:", rmse_AddiVortes))
print(paste("Accuracy for AddiVortes:", Accuracy_AddiVortes))
sqrt(mean((mean(test_data$Y)-Y_test_mu)^2))
1-sum(as.numeric(test_data$Y))/length(test_data$Y)
plot(sort(pnorm(prediction_AddiVortes)),col = ifelse(test_data$Y[order(prediction_AddiVortes)]==1,'red','black'))
abline(h=0.5)

cat("Fitting local update algorithm...\n")
time_soft <- system.time({
  set.seed(123)
  Model_softBART<-softbart(X_train,training_data$Y,X_test,hypers = Hypers(X_train, training_data$Y,num_tree = 20))#,opts = Opts(num_burn = 200, num_save = 1000))
})
print(time_soft)
plot(Model_softBART$sigma^2, type = "l", col = "blue", lwd = 2)
abline(h = 10, col = "red")
print(Model_softBART$alpha)
plot(colMeans(Model_softBART$s))
print(paste("RMSE for SoftBART:", sqrt(mean((Model_softBART$y_hat_test_mean - Y_test_mu)^2))))

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

      