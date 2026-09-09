#Friedman dataset
{
  
  sim_fried <- function(N,P,sigma) {
    X <- matrix(runif(N * P), nrow = N, ncol = P)
    mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)**2 + 10 * X[,4] + 5 * X[,5]
    Y <- mu + sigma * rnorm(N)
    
    return(list(X = X, Y = Y,mu=mu))
  }
  
  set.seed(453)
  training_data <- sim_fried(200,400 , sqrt(10))
  test_data <- sim_fried(200,400, sqrt(10))
  X_train <- training_data$X
  X_test <-test_data$X
  Y_test_mu<-test_data$mu
}

set.seed(136)

tau<-50
boost_val<-0.5
penalty_val<-0.5
decay_val<-0.7
alpha_val<-1
a_alpha_val <- 0.5
b_alpha_val <- 1

Model_AddiVortes_local <- AddiVortes(training_data$Y, X_train,m=200,thinning = 1,varSelMode = 2,Omega=1,LambdaRate=4,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = TRUE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = 0.6,numChains = 1,IntialSigma = "LASSO",tau=tau,splitMode =1)#,rho_alpha = 1)#,power = p_init_val, p_shape = p_shape_val, p_rate = p_rate_val, p_sd = p_sd_val)

Model_AddiVortes_original <- AddiVortes(training_data$Y, X_train,m=200,thinning = 1,varSelMode = 0,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = FALSE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = 0.6,numChains = 1,IntialSigma = "LASSO",tau=tau,splitMode =1)#,rho_alpha = 1)#,power = p_init_val, p_shape = p_shape_val, p_rate = p_rate_val, p_sd = p_sd_val)

Model_AddiVortes_dir <- AddiVortes(training_data$Y, X_train,m=200,thinning = 1,varSelMode = 1,totalMCMCIter = 5000, mcmcBurnIn = 2500,dirichletWarmup = 1000,nu=6,q=0.9,updateAlpha = TRUE,alpha=alpha_val,a_alpha = a_alpha_val,b_alpha=b_alpha_val,adaptBoost = boost_val, adaptPenalty = penalty_val, momentumDecay = decay_val, kappa = 0.6,numChains = 1,IntialSigma = "LASSO",tau=tau,splitMode =1)#,rho_alpha = 1)#,power = p_init_val, p_shape = p_shape_val, p_rate = p_rate_val, p_sd = p_sd_val)

Model_DART<-wbart(X_train,training_data$Y,X_test,ntree = 200,ndpost = 2500, nskip = 2500,sparse = TRUE)

Model_BART<-wbart(X_train,training_data$Y,X_test,ntree = 200,ndpost = 2500, nskip = 2500,sparse = FALSE)


col_adaptive <- "darkorange"
col_original <- "firebrick"
col_dart <- "forestgreen"
col_dir <- "darkorchid"
col_bart <- "royalblue"



#### Graph 1 Sigma trace ####
{par(mfrow=c(1,1))
  par(mgp = c(2.5, 1, 0))
  plot(Model_AddiVortes_local$posteriorSigma, type = "l", col = col_adaptive, lwd = 2,
       main = expression("Posterior"~ sigma^2~" Trace Plot"),
       xlab = "MCMC Iteration",ylab = expression(sigma^2 ~ "sample"), ylim = c(1,12))
  lines(Model_AddiVortes_original$posteriorSigma, col= col_original)
  lines(Model_AddiVortes_dir$posteriorSigma,  col = col_dir)
  lines(Model_DART$sigma[2500:5000]^2, col = col_dart)
  #lines(Model_BART$sigma[2500:5000]^2, col = col_bart)
  abline(10,0,)
  legend("bottomright", legend = c("Adaptive AddiVortes","Dirichlet AddiVortes", "Original AddiVortes", "DART" ,expression("True"~sigma^2~"value (10)")),
         col = c(col_adaptive, col_dir,col_original, col_dart,col_bart, "black"), 
         lwd = c(2, 2, 2, 2,2), lty = c(1, 1, 1, 1,1), bty = "n", cex = 0.8)
  par(mgp = c(3, 1, 0))
}

### Graph 2 Average #Centres/#Dimensions per tessellation
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
       xlab = "MCMC Iteration", ylab = "Average number of dimensions",ylim=c(0.5,5))
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

### Graph 3 Augmented counts ###
{
  old_par <- par(no.readonly = TRUE)
  par(mfrow=c(1,2))
  
  # --- Model 1: Dirichlet AddiVortes ---
  mean_aug_counts_dir <- rowMeans(Model_AddiVortes_dir$posteriorAugmentedCounts)
  lower_aug_counts_dir <- apply(Model_AddiVortes_dir$posteriorAugmentedCounts, 1, quantile, probs = 0.025)
  upper_aug_counts_dir <- apply(Model_AddiVortes_dir$posteriorAugmentedCounts, 1, quantile, probs = 0.975)
  
  cov_names_dir <- if(!is.null(names(Model_AddiVortes_dir$xCentres))) names(Model_AddiVortes_dir$xCentres) else paste0("X", 1:length(mean_aug_counts_dir))
  ord_dir <- order(mean_aug_counts_dir, decreasing = FALSE)
  
  max_vars_to_plot <- 20
  if (length(mean_aug_counts_dir) > max_vars_to_plot) {
    ord_dir <- tail(ord_dir, max_vars_to_plot)
    main_title_dir <- "Dirichlet AddiVortes"
  } else {
    main_title_dir <- "Mean Augmented Geometric Failures"
  }
  
  max_name_len_dir <- max(nchar(cov_names_dir[ord_dir]))
  par(mar = c(5, max(4.1, max_name_len_dir * 0.6), 4, 2) + 0.1)
  
  # Extend xlim to ensure error bars fit on the plot
  xlim_dir <- c(0, max(upper_aug_counts_dir[ord_dir]) * 1.05)
  
  bp_dir <- barplot(mean_aug_counts_dir[ord_dir], horiz = TRUE, names.arg = cov_names_dir[ord_dir], las = 1,
                    col = col_dir, border = NA, xlab = "Failure Count",
                    main = main_title_dir, cex.names = 0.8, xlim = xlim_dir)
  
  # Add 95% Credible Intervals
  segments(x0 = lower_aug_counts_dir[ord_dir], y0 = bp_dir, 
           x1 = upper_aug_counts_dir[ord_dir], y1 = bp_dir, 
           col = "black", lwd = 1.5)
  
  # Add vertical caps to the intervals
  epsilon <- 0.15
  segments(x0 = lower_aug_counts_dir[ord_dir], y0 = bp_dir - epsilon, 
           x1 = lower_aug_counts_dir[ord_dir], y1 = bp_dir + epsilon, col = "black", lwd = 1.5)
  segments(x0 = upper_aug_counts_dir[ord_dir], y0 = bp_dir - epsilon, 
           x1 = upper_aug_counts_dir[ord_dir], y1 = bp_dir + epsilon, col = "black", lwd = 1.5)
  
  # --- Model 2: Adaptive AddiVortes ---
  mean_aug_counts_local <- rowMeans(Model_AddiVortes_local$posteriorAugmentedCounts)
  lower_aug_counts_local <- apply(Model_AddiVortes_local$posteriorAugmentedCounts, 1, quantile, probs = 0.025)
  upper_aug_counts_local <- apply(Model_AddiVortes_local$posteriorAugmentedCounts, 1, quantile, probs = 0.975)
  
  cov_names_local <- if(!is.null(names(Model_AddiVortes_local$xCentres))) names(Model_AddiVortes_local$xCentres) else paste0("X", 1:length(mean_aug_counts_local))
  ord_local <- order(mean_aug_counts_local, decreasing = FALSE)
  
  if (length(mean_aug_counts_local) > max_vars_to_plot) {
    ord_local <- tail(ord_local, max_vars_to_plot)
    main_title_local <- "Adaptive AddiVortes"
  } else {
    main_title_local <- "Mean Augmented Geometric Failures"
  }
  
  max_name_len_local <- max(nchar(cov_names_local[ord_local]))
  par(mar = c(5, max(4.1, max_name_len_local * 0.6), 4, 2) + 0.1)
  
  xlim_local <- c(0, max(upper_aug_counts_local[ord_local]) * 1.05)
  
  bp_local <- barplot(mean_aug_counts_local[ord_local], horiz = TRUE, names.arg = cov_names_local[ord_local], las = 1,
                      col = col_adaptive, border = NA, xlab = "Failure Count",
                      main = main_title_local, cex.names = 0.8, xlim = xlim_local)
  
  # Add 95% Credible Intervals
  segments(x0 = lower_aug_counts_local[ord_local], y0 = bp_local, 
           x1 = upper_aug_counts_local[ord_local], y1 = bp_local, 
           col = "black", lwd = 1.5)
  
  # Add vertical caps to the intervals
  segments(x0 = lower_aug_counts_local[ord_local], y0 = bp_local - epsilon, 
           x1 = lower_aug_counts_local[ord_local], y1 = bp_local + epsilon, col = "black", lwd = 1.5)
  segments(x0 = upper_aug_counts_local[ord_local], y0 = bp_local - epsilon, 
           x1 = upper_aug_counts_local[ord_local], y1 = bp_local + epsilon, col = "black", lwd = 1.5)
  
  par(old_par)
}

### Graph 4 TRUE vs Predicted
{
  # 1. Compute posterior predictive means for each model
  # Adjust the slot names if your object structures store test predictions differently
  pred_addivortes_local <- predict(Model_AddiVortes_local,as.matrix(X_test))
  pred_addivortes_original <- predict(Model_AddiVortes_original,as.matrix(X_test))
  pred_addivortes_dir   <- predict(Model_AddiVortes_dir,as.matrix(X_test))
  pred_dart            <- colMeans(Model_DART$yhat.test)
  
  # 2. Configure a 2x2 plotting grid with shared graphical parameters
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  
  # Determine common axis limits for consistency across all four plots
  all_values <- c(Y_test_mu, pred_addivortes_local, pred_addivortes_dir, pred_dart)
  plot_limits <- range(all_values, na.rm = TRUE)
  
  # 3. Plot 1: AddiVortes Local
  plot(Y_test_mu, pred_addivortes_local ,
       main = "Adaptive AddiVortes",
       xlab = "True values",
       ylab = "Predicted values",
       xlim = plot_limits, ylim = plot_limits,
       col = col_adaptive, pch = 16)
  abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = 2)
  
  # 4. Plot 2: AddiVortes Dirichlet
  plot(Y_test_mu, pred_addivortes_dir,
       main = "AddiVortes Dirichlet",
       xlab = "True values",
       ylab = "",
       xlim = plot_limits, ylim = plot_limits,
       col = col_dir, pch = 16)
  abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = 2)
  
  # 5. Plot 3: DART
  plot(Y_test_mu, pred_dart,
       main = "DART",
       xlab = "True values",
       ylab = "Predicted values",
       xlim = plot_limits, ylim = plot_limits,
       col = col_dart, pch = 16)
  abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = 2)
  
  plot(Y_test_mu, pred_addivortes_original,
       main = "Original AddiVortes",
       xlab = "True values",
       ylab = "",
       xlim = plot_limits, ylim = plot_limits,
       col = col_original, pch = 16)
  abline(a = 0, b = 1, col = "darkgrey", lty = 2, lwd = 2)
  
  # Reset plotting parameters to default
  par(mfrow = c(1, 1))
  
  print(sqrt(mean((Y_test_mu-pred_addivortes_local)^2)))
  print(sqrt(mean((Y_test_mu-pred_addivortes_original)^2)))
  print(sqrt(mean((Y_test_mu-pred_addivortes_dir)^2)))
  print(sqrt(mean((Y_test_mu-pred_dart)^2)))
}

### Graph 5 Number of Active Dimensions ###
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
       xlab = "MCMC Iteration", ylab = "Number of Active Covariates",
       main = "Trace: Active Dimensions")
  lines(active_orig, col = col_original, lwd = 2)
  lines(active_dart, col = col_dart, lwd = 2)
  lines(active_dir, col = col_dir, lwd = 2)
  lines(active_bart, col = col_bart, lwd = 2)
  
  
  # Add ground truth line (5 active variables)
  abline(h = 5, col = "black", lwd = 2, lty = 2)
  legend("right", legend = c("Adaptive AddiVortes","Dirichlet AddiVortes", "Original AddiVortes", "DART","BART" ,"True Sparsity (5)"),
         col = c(col_adaptive, col_dir,col_original, col_dart,col_bart, "black"), 
         lwd = c(2, 2, 2, 2,2,2), lty = c(1, 1, 1, 1,1,2), bty = "n", cex = 0.8)
  
  plot(active_adapt, type = "l", col = col_adaptive, lwd = 2,
       ylim = c(0, max(active_adapt, active_dart)),
       xlab = "MCMC Iteration", ylab = "Number of Active Covariates",
       main = "Trace: Active Dimensions")
  lines(active_dart, col = col_dart, lwd = 2)
  lines(active_dir, col = col_dir, lwd = 2)
  
  # Add ground truth line (5 active variables)
  abline(h = 5, col = "black", lwd = 2, lty = 2)
  legend("topright", legend = c("Adaptive AddiVortes","Dirichlet AddiVortes", "DART" ,"True Sparsity (5)"),
         col = c(col_adaptive, col_dir, col_dart, "black"), 
         lwd = c(2, 2,2, 2), lty = c(1, 1,1,2), bty = "n", cex = 0.8)
}

### Graph 6 Posterior Dirichlet Weights ###
{
  par(mfrow=c(1,3))
  {# --- 1. Extract the Dirichlet Weight metrics ---
    s_means <- Model_AddiVortes_local$posteriorDirichletWeightsMean
    s_lower <- Model_AddiVortes_local$posteriorDirichletWeightsLower
    s_upper <- Model_AddiVortes_local$posteriorDirichletWeightsUpper
    covariate_index <- 1:length(s_means)
    
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
                  col = col_adaptive,
                  main = "Adaptive AddiVortes",
                  xlab = "Covariate Index", 
                  ylab = expression(s[j]),ylim=c(0,1))
    
    # 3. Add the error bars using the arrows() function
    arrows(x0 = bp, y0 = y_lower, x1 = bp, y1 = y_upper, 
           angle = 90, code = 3, length = 0.05, col = "black")
    
  }
  
  {# --- 1. Extract the Dirichlet Weight metrics ---
    s_means <- Model_AddiVortes_dir$posteriorDirichletWeightsMean
    s_lower <- Model_AddiVortes_dir$posteriorDirichletWeightsLower
    s_upper <- Model_AddiVortes_dir$posteriorDirichletWeightsUpper
    covariate_index <- 1:length(s_means)
    
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
                  col = col_dir,
                  main = "Dirchlet AddiVortes",
                  xlab = "Covariate Index", 
                  ylab = NULL,ylim=c(0,1))
    
    # 3. Add the error bars using the arrows() function
    arrows(x0 = bp, y0 = y_lower, x1 = bp, y1 = y_upper, 
           angle = 90, code = 3, length = 0.05, col = "black")
  }
  
  {
    # 1. Extract and calculate the proportions matrix
    split_counts <- Model_DART$varcount
    proportions_matrix <- split_counts / rowSums(split_counts)
    
    # 2. Extract means and 95% credible intervals
    inc_means <- colMeans(proportions_matrix)
    inc_lower <- apply(proportions_matrix, 2, quantile, probs = 0.025)
    inc_upper <- apply(proportions_matrix, 2, quantile, probs = 0.975)
    
    covariate_index_dart <- 1:length(inc_means)
    number_points <- 10
    
    # 3. Filter for the first 5 and the top 10 of the rest
    ordered_indexes_dart <- order(inc_means[-c(1:5)], decreasing = TRUE)[1:number_points] + 5
    
    x_vals_dart <- c(1:5, covariate_index_dart[ordered_indexes_dart])
    y_vals_dart <- c(inc_means[1:5], inc_means[ordered_indexes_dart])
    y_lower_dart <- c(inc_lower[1:5], inc_lower[ordered_indexes_dart])
    y_upper_dart <- c(inc_upper[1:5], inc_upper[ordered_indexes_dart])
    
    # 4. Draw the bar chart
    y_max <- max(y_upper_dart, na.rm = TRUE) * 1.1 
    
    bp_dart <- barplot(y_vals_dart, 
                       names.arg = x_vals_dart, 
                       col = col_dart,
                       main = "DART",
                       xlab = "Covariate Index", 
                       ylim = c(0, 1))
    
    # 5. Add the error bars
    arrows(x0 = bp_dart, y0 = y_lower_dart, x1 = bp_dart, y1 = y_upper_dart, 
           angle = 90, code = 3, length = 0.05, col = "black")
  }
  par(mfrow=c(1,1))
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

### Graph 7 Ensemble inclusion probabilities ###
{
  par(mfrow=c(2,2))
  {
    inclusion_probs <- Model_AddiVortes_local$ensembleInclusionProbabilities
    bars_displayed <- 10
    
    # 1. Identify the original indices of the top noise variables
    noise_indices <- 6:length(inclusion_probs)
    noise_probs <- inclusion_probs[noise_indices]
    
    # Use order() to get the indices of the highest noise variables, then map back to original index
    top_noise_order <- order(noise_probs, decreasing = TRUE)[1:bars_displayed]
    top_noise_actual_indices <- noise_indices[top_noise_order]
    
    # 2. Combine the values and names for the plot
    plot_values <- c(inclusion_probs[1:5], inclusion_probs[top_noise_actual_indices])
    plot_names <- c(1:5, top_noise_actual_indices)
    
    # 3. Generate the corrected barplot
    barplot(plot_values, 
            names.arg = plot_names,
            main = "Adaptive AddiVortes",
            xlab = "Covariate Index", 
            ylab = "Proportion of inclusion",
            col = col_adaptive,
            border = NA,
            ylim = c(0, 1))
  }
  
  {
    inclusion_probs <- Model_AddiVortes_dir$ensembleInclusionProbabilities
    bars_displayed <- 10
    
    # 1. Identify the original indices of the top noise variables
    noise_indices <- 6:length(inclusion_probs)
    noise_probs <- inclusion_probs[noise_indices]
    
    # Use order() to get the indices of the highest noise variables, then map back to original index
    top_noise_order <- order(noise_probs, decreasing = TRUE)[1:bars_displayed]
    top_noise_actual_indices <- noise_indices[top_noise_order]
    
    # 2. Combine the values and names for the plot
    plot_values <- c(inclusion_probs[1:5], inclusion_probs[top_noise_actual_indices])
    plot_names <- c(1:5, top_noise_actual_indices)
    
    # 3. Generate the corrected barplot
    barplot(plot_values, 
            names.arg = plot_names,
            main = "Dirichlet AddiVortes",
            xlab = "Covariate Index", 
            col = col_dir,
            border = NA,
            ylim = c(0, 1))
  }
  
  {
    # 1. Extract the variable count matrix
    split_counts_wbart <- Model_DART$varcount
    
    # 2. Calculate probability of inclusion in the ensemble (count > 0)
    prob_in_ensemble <- colMeans(split_counts_wbart > 0)
    
    # 3. Define the number of top additional bars to display
    bars_displayed <- 10
    
    # 4. Identify the indices of the top remaining covariates
    ordered_indexes_ensemble <- order(prob_in_ensemble[-c(1:5)], decreasing = TRUE)[1:bars_displayed] + 5
    
    # 5. Construct the vectors for plotting
    y_vals_ensemble <- c(prob_in_ensemble[1:5], prob_in_ensemble[ordered_indexes_ensemble])
    x_names_ensemble <- c(1:5, ordered_indexes_ensemble)
    
    # 6. Generate the bar plot
    barplot(y_vals_ensemble,
            main = "DART",
            xlab = "Covariate Index", 
            ylab = "Proportion of inclusion",
            col = col_dart,
            border = NA,
            ylim = c(0, 1),
            names.arg = x_names_ensemble)
  }
  
  {
    inclusion_probs <- Model_AddiVortes_original$ensembleInclusionProbabilities
    bars_displayed <- 10
    
    # 1. Identify the original indices of the top noise variables
    noise_indices <- 6:length(inclusion_probs)
    noise_probs <- inclusion_probs[noise_indices]
    
    # Use order() to get the indices of the highest noise variables, then map back to original index
    top_noise_order <- order(noise_probs, decreasing = TRUE)[1:bars_displayed]
    top_noise_actual_indices <- noise_indices[top_noise_order]
    
    # 2. Combine the values and names for the plot
    plot_values <- c(inclusion_probs[1:5], inclusion_probs[top_noise_actual_indices])
    plot_names <- c(1:5, top_noise_actual_indices)
    
    # 3. Generate the corrected barplot
    barplot(plot_values, 
            names.arg = plot_names,
            main = "Original AddiVortes",
            xlab = "Covariate Index", 
            col = col_original,
            border = NA,
            ylim = c(0, 1))
  }
  
  par(mfrow=c(1,1))
}

### Graph 7.5 trace of variable count ###
{
  par(mfrow=c(2,5))
  for (j in 1:10) {
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
}

### Graph 8 Alpha trace ###

{
  par(mfrow=c(1,2))
  plot(1:length(Model_AddiVortes_local$posteriorAlpha), Model_AddiVortes_local$posteriorAlpha, type = "l", xlab = "MCMC Iteration",
       ylab = expression(alpha), main = "Adaptive Alpha Trace",
       col = col_adaptive, lwd = 1.5)
  abline(h = mean(Model_AddiVortes_local$posteriorAlpha, na.rm=TRUE), col = "red", lty = 2)
  plot(1:length(Model_AddiVortes_dir$posteriorAlpha), Model_AddiVortes_dir$posteriorAlpha, type = "l", xlab = "MCMC Iteration",
       ylab = expression(alpha), main = "Dirichlet Alpha Trace",
       col = col_dir, lwd = 1.5)
  abline(h = mean(Model_AddiVortes_dir$posteriorAlpha, na.rm=TRUE), col = "red", lty = 2)
  par(mfrow=c(1,1))
}

### Graph 9 Momentum barchart ###
{
  # 1. Calculate the mean, lower (2.5%), and upper (97.5%) bounds
  mean_momentum <- rowMeans(Model_AddiVortes_local$posteriorMomentum)
  lower_momentum <- apply(Model_AddiVortes_local$posteriorMomentum, 1, quantile, probs = 0.025)
  upper_momentum <- apply(Model_AddiVortes_local$posteriorMomentum, 1, quantile, probs = 0.975)
  
  cov_names <- if(!is.null(names(Model_AddiVortes_local$xCentres))) names(Model_AddiVortes_local$xCentres) else paste0("X", 1:length(mean_momentum))
  
  # 2. Sort based on the mean
  ord <- order(mean_momentum, decreasing = FALSE)
  
  max_vars_to_plot <- 20
  if (length(mean_momentum) > max_vars_to_plot) {
    ord <- tail(ord, max_vars_to_plot)
    main_title <- paste("Top", max_vars_to_plot, "Mean Adaptive Momentum (95% CI)")
  } else {
    main_title <- "Mean Adaptive Momentum (95% CI)"
  }
  
  max_name_len <- max(nchar(cov_names[ord]))
  par(mar = c(5, max(4.1, max_name_len * 0.6), 4, 2) + 0.1)
  
  # 3. Calculate dynamic xlim to ensure error bars aren't cut off
  x_min <- min(0, min(lower_momentum[ord])) 
  x_max <- max(upper_momentum[ord])
  
  # 4. Save the bar midpoints into 'bp'
  bp <- barplot(mean_momentum[ord], horiz = TRUE, names.arg = cov_names[ord], las = 1,
                col = col_adaptive, border = NA, xlab = "Mean Momentum Term",
                main = main_title, cex.names = 0.8,
                xlim = c(x_min, x_max * 1.05)) # Added dynamic x-axis limit
  
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
}

### Graph 9.5 Momentum traces ###
{
  par(mfrow=c(2,5))
  plot(Model_AddiVortes_local$posteriorMomentum[1,], type = "l", col = col_adaptive,
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Posterior Momentum for x_1")
  plot(Model_AddiVortes_local$posteriorMomentum[2,], type = "l", col = col_adaptive,
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Posterior Momentum")
  plot(Model_AddiVortes_local$posteriorMomentum[3,], type = "l", col = col_adaptive,
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Posterior Momentum for x_1")
  plot(Model_AddiVortes_local$posteriorMomentum[4,], type = "l", col = col_adaptive,
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Posterior Momentum for x_1")
  plot(Model_AddiVortes_local$posteriorMomentum[5,], type = "l", col = col_adaptive,
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Posterior Momentum for x_1")
  plot(Model_AddiVortes_local$posteriorMomentum[10,], type = "l", col = col_adaptive,
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Posterior Momentum for x_1")
  plot(Model_AddiVortes_local$posteriorMomentum[159,], type = "l", col = col_adaptive,
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Posterior Momentum for x_1")
  plot(Model_AddiVortes_local$posteriorMomentum[475,], type = "l", col = col_adaptive,
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Posterior Momentum for x_1")
  plot(Model_AddiVortes_local$posteriorMomentum[271,], type = "l", col = col_adaptive,
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Posterior Momentum for x_1")
  plot(Model_AddiVortes_local$posteriorMomentum[222,], type = "l", col = col_adaptive,
       main = "Posterior Alpha Trace",
       xlab = "MCMC Iteration", ylab = "Posterior Momentum for x_1")
  
  par(mfrow=c(1,1))
}

### Graph 10 Probability 
{
  par(mfrow=c(1,2))
  sum_active_covariates <- colSums(Model_AddiVortes_local$posteriorDirichletWeights[1:5, ])
  sum_noise_covariates <- colSums(Model_AddiVortes_local$posteriorDirichletWeights[6:500, ])
  plot(sum_active_covariates, type = "l", col = col_adaptive,
       main = "Posteriot Weights",ylim = c(0,1) )
  lines(sum_noise_covariates,col='red')
  abline(a=mean(sum_active_covariates),b=0)
  abline(a=mean(sum_noise_covariates),b=0)
  sum_active_covariates <- colSums(Model_AddiVortes_dir$posteriorDirichletWeights[1:5, ])
  sum_noise_covariates <- colSums(Model_AddiVortes_dir$posteriorDirichletWeights[6:500, ])
  plot(sum_active_covariates, type = "l", col = col_dir,
       main = "Posteriot Weights",ylim = c(0,1) )
  lines(sum_noise_covariates,col='red')
  abline(a=mean(sum_active_covariates),b=0)
  abline(a=mean(sum_noise_covariates),b=0)
  
  par(mfrow=c(1,1))
}












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