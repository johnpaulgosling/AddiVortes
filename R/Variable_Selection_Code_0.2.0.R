if (!requireNamespace("invgamma", quietly = TRUE)) {
  install.packages("invgamma")
}

if (!requireNamespace("FNN", quietly = TRUE)) {
  install.packages("FNN")
}

library('invgamma')
library('FNN')

rdirichlet_custom <- function(alpha_vec) {
  y <- rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)
  return(y / sum(y))
}

Homo_AddiVortes_Algorithm<-function(y,x,m = 200 ,max_iter = 1200,burn_in= 200,nu = 6,q =0.85,k = 3 ,sd = 0.8 ,Omega = 3,lambda_rate = 25,YTest,XTest,IntialSigma = "Linear",thinning=1, alpha=1, a_alpha=0.5, b_alpha=1){
  
  yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
  xScaled=x;
  for (i in 1:length(x[1,])){
    xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  for (i in 1:length(XTest[1,])){
    XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  Pred<-rep(list(matrix(mean(yScaled)/m)),m)
  Dim=vector(length = m)
  Tess=vector(length = m)
  for (i in 1:m){
    Dim[i]<-list(sample(1:length(x[1,]), 1))
    Tess[i]<-(list(matrix(rnorm(1,0,sd))))
  }
  
  p_cov <- length(x[1,])
  rho_alpha <- p_cov
  s_weights <- rep(1 / p_cov, p_cov)
  dirichlet_warmup <- floor(burn_in / 2)
  
  SumOfAllTess=rep(mean(yScaled),length(yScaled))
  SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
  LastTessPred=matrix
  TessIndexes=matrix(1,nrow=length(yScaled),ncol=m)
  
  posterior_samples<-floor((max_iter - burn_in) / thinning)+1
  PredictionMatrix<-array(dim=c(length(y),posterior_samples))
  TestMatrix<-array(dim=c(length(YTest),posterior_samples))
  
  Dirichlet_Weights_Matrix<-array(0, dim=c(p_cov, posterior_samples))
  Variable_Selection_Matrix<-array(0, dim=c(p_cov, posterior_samples)) 
  Alpha_Trace<-vector(length = posterior_samples)
  
  if (IntialSigma=="Naive"){ 
    SigmaSquaredHat=var(yScaled)
  } else {  
    MultiLinear<-lm(yScaled ~ xScaled)
    SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
  }
  Sigma_squared_test<-vector(length = (max_iter-burn_in)/thinning)
  
  lambda=1;
  lambda <- optim(par = 1,
                  fitting_function,
                  method = "Brent",
                  lower = 0.001,
                  upper = 100,
                  q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
  
  for (i in 1:max_iter){
    
    SigmaSquared<-Homo_SigmaSquaredCalculation(yScaled,nu,lambda,SumOfAllTess)
    
    if (i>burn_in ){ 
      Sigma_squared_test[i-burn_in]<-SigmaSquared*((max(y)-min(y))^2)
    }
    
    for (j in 1:m){
      NewTessOutput<-Homo_NewTess(xScaled,j,Tess,Dim,sd,s_weights) 
      TessStar<-NewTessOutput[[1]]  
      DimStar<-NewTessOutput[[2]]
      Modification<-NewTessOutput[[3]]
      
      ResidualsOutput<-Homo_CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred,as.numeric(TessIndexes[,j])) 
      R_ijOld<-ResidualsOutput[[1]]   
      n_ijOld<-ResidualsOutput[[2]]
      R_ijNew<-ResidualsOutput[[3]]
      n_ijNew<-ResidualsOutput[[4]]
      SumOfAllTess<-ResidualsOutput[[5]] 
      IndexesStar<-ResidualsOutput[[6]] 
      
      LOGAcceptenceProb=Homo_AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) 
      
      if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ 
        Tess=TessStar
        Dim=DimStar
        Pred[[j]]=Homo_NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
        LastTessPred=Pred[[j]][IndexesStar]
        TessIndexes[,j]<- as.numeric(IndexesStar)
      } else { 
        Pred[[j]]=Homo_NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared)
        LastTessPred=Pred[[j]][TessIndexes[,j]]
      }
      
      if (j==m){ 
        SumOfAllTess=SumOfAllTess+LastTessPred;
      }
    }
    
    if (i > dirichlet_warmup) {
      m_counts <- rep(0, p_cov)
      for (j in 1:m) {
        current_dims <- Dim[[j]]
        m_counts[current_dims] <- m_counts[current_dims] + 1
        
        d <- length(current_dims)
        if (d > 1) {
          ordered_dims <- sample(current_dims)
          
          for (k in 2:d) {
            ineligible <- ordered_dims[1:(k-1)]
            prob_success <- 1 - sum(s_weights[ineligible])
            
            if (prob_success > 0 && prob_success < 1) {
              failures <- rgeom(1, prob_success)
              
              if (failures > 0) {
                if (length(ineligible) == 1) {
                  m_counts[ineligible] <- m_counts[ineligible] + failures
                } else {
                  failed_draws <- sample(ineligible, failures, replace = TRUE, prob = s_weights[ineligible])
                  tab <- table(factor(failed_draws, levels = 1:p_cov))
                  m_counts <- m_counts + as.numeric(tab)
                }
              }
            }
          }
        }
      }
      s_weights <- rdirichlet_custom(alpha / p_cov + m_counts)
      
      alpha_star <- exp(rnorm(1, mean = log(alpha), sd = 0.5))
      
      ll_alpha <- function(a_val) {
        safe_s_weights <- pmax(s_weights, 1e-16)
        lgamma(a_val) - p_cov * lgamma(a_val / p_cov) + (a_val / p_cov) * sum(log(safe_s_weights)) + (a_alpha - 1) * log(a_val) - (a_alpha + b_alpha) * log(a_val + rho_alpha)
      }
      
      log_accept_ratio <- ll_alpha(alpha_star) - ll_alpha(alpha) + log(alpha_star) - log(alpha)
      
      if (log(runif(1)) < log_accept_ratio) {
        alpha <- alpha_star
      }
      
    } else {
      s_weights <- rep(1 / p_cov, p_cov)
    }
    
    if (i %% 100 == 0){
      cat(sprintf("Iteration %d out of %d", i, max_iter), "\n")
    }
    
    if (i>=burn_in & (i-burn_in) %% thinning ==0 ){ 
      current_sample_index <- 1 + (i - burn_in) / thinning
      PredictionMatrix[,current_sample_index]=SumOfAllTess;
      TestMatrix[,current_sample_index]=Homo_TestPrediction(XTest,m,Tess,Dim,Pred);
      
      Dirichlet_Weights_Matrix[,current_sample_index] <- s_weights
      Alpha_Trace[current_sample_index] <- alpha
      
      current_covariate_counts <- rep(0, p_cov)
      for (j in 1:m) {
        current_covariate_counts[Dim[[j]]] <- current_covariate_counts[Dim[[j]]] + 1
      }
      Variable_Selection_Matrix[,current_sample_index] <- current_covariate_counts
    }
  }
  
  mean_yhat=(rowSums(PredictionMatrix)/(posterior_samples))*(max(y)-min(y))+((max(y)+min(y))/2)
  mean_yhat_Test=(rowSums(TestMatrix)/(posterior_samples))*(max(y)-min(y))+((max(y)+min(y))/2)
  
  inclusion_indicator_matrix <- Variable_Selection_Matrix > 0
  ensemble_inclusion_probabilities <- rowMeans(inclusion_indicator_matrix)
  
  return( 
    list(
      In_sample_RMSE = sqrt(mean((y-mean_yhat)^2)),
      Out_of_sample_RMSE = sqrt(mean((YTest-mean_yhat_Test)^2)),
      SMSE=mean((YTest-mean_yhat_Test)^2)/mean((YTest-mean(y))^2),
      y_yest_hat=mean_yhat_Test,
      Posterior_Dirichlet_Weights_Mean = rowMeans(Dirichlet_Weights_Matrix),
      Posterior_Dirichlet_Weights_Lower = apply(Dirichlet_Weights_Matrix, 1, quantile, probs = 0.025),
      Posterior_Dirichlet_Weights_Upper = apply(Dirichlet_Weights_Matrix, 1, quantile, probs = 0.975),
      Ensemble_Inclusion_Probabilities = ensemble_inclusion_probabilities,
      Posterior_Alpha_Trace = Alpha_Trace
    )
  )
}
Homo_SigmaSquaredCalculation<-function(yScaled,nu,lambda,SumOfAllTess){ 
  n=length(yScaled)
  SigmaSquared<-rinvgamma(1,shape=(nu+n)/2,rate=(nu*lambda+sum((yScaled-SumOfAllTess)^2))/2)
  return(SigmaSquared)
}

Homo_NewTess<-function(x,j,Tess,Dim,sd,s_weights){ 
  p=runif(1,0,1) 
  DimStar=Dim 
  TessStar=Tess 
  
  if (p<0.2 & length(Dim[[j]])!=length(x[1,]) | length(Dim[[j]])==1 & p<0.4){ 
    NumberOfCovariates=1:length(x[1,]) 
    NumberOfCovariates=NumberOfCovariates[-Dim[[j]]] 
    if (length(NumberOfCovariates) == 1) {
      added_dim <- NumberOfCovariates
    } else {
      added_dim <- sample(NumberOfCovariates, 1, prob = s_weights[NumberOfCovariates])
    }
    DimStar[[j]]<-c(Dim[[j]], added_dim) 
    TessStar[[j]]=cbind(Tess[[j]],rnorm(length(Tess[[j]][,1]),0,sd)) 
    Modification="AD"
  } else if (p<0.4){ 
    RemovedDim=sample(1:length(Dim[[j]]),1) 
    DimStar[[j]]=DimStar[[j]][-RemovedDim] 
    TessStar[[j]]=matrix(TessStar[[j]][,-RemovedDim],ncol=length(DimStar[[j]])) 
    Modification="RD"
  } else if (p<0.6 || p<0.8 & length(Tess[[j]][,1])==1){ 
    TessStar[[j]]=rbind(Tess[[j]],rnorm(length(Dim[[j]]),0,sd)) 
    Modification="AC"
  } else if (p<0.8){ 
    CenterRemoved=sample(1:length(TessStar[[j]][,1]),1) 
    TessStar[[j]]=matrix(TessStar[[j]][-CenterRemoved,],ncol=length(Dim[[j]])) 
    Modification="RC"
  } else if (p<0.9 || length(Dim[[j]])==length(x[1,])){ 
    TessStar[[j]][sample(1:length(TessStar[[j]][,1]),1),]=rnorm(length(Dim[[j]]),0,sd) 
    Modification="Change"
  } else { 
    NumberOfCovariates=1:length(x[1,])  
    NumberOfCovariates=NumberOfCovariates[-Dim[[j]]]  
    DimToChange=sample(1:length(Dim[[j]]),1) 
    if (length(NumberOfCovariates) == 1) {
      added_dim <- NumberOfCovariates
    } else {
      added_dim <- sample(NumberOfCovariates, 1, prob = s_weights[NumberOfCovariates])
    }
    DimStar[[j]][DimToChange]=added_dim 
    TessStar[[j]][,DimToChange]=rnorm(length(Tess[[j]][,1]),0,sd) 
    Modification="Swop"
  }
  
  TessStar[[j]]<-matrix(TessStar[[j]],ncol=length(DimStar[[j]])) 
  return(list(TessStar,DimStar,Modification)) 
}

fitting_function<- function(lambda,q,nu,sigmaSquared_hat){ 
  return((sigmaSquared_hat- qinvgamma(q, shape=nu/2, rate=nu*lambda/2))^2)
}

Homo_Indexes<-function(x,Tess,Dim){ 
  if (length(Tess[,1])==1){ 
    CellsForGivenTess=rep(1,length(x[,1]))
  } else { 
    CellsForGivenTess=knnx.index(Tess,matrix(x[,Dim],ncol = length(Dim)),1)
  }
  return(CellsForGivenTess)
}

Homo_AlphaCalculation<-function(x,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate){ 
  
  d=length(DimStar[[j]]);
  NumCovariates=length(x[1,]);
  cStar=length(TessStar[[j]][,1]);
  
  LOGlikelihoodRatio=0.5*(log(prod(n_ijOld*SigmaSquaredMu+SigmaSquared))-log(prod(n_ijNew*SigmaSquaredMu+SigmaSquared)))+((SigmaSquaredMu/(2*SigmaSquared))*(-sum((R_ijOld^2)/(n_ijOld*SigmaSquaredMu+SigmaSquared))+sum((R_ijNew^2)/(n_ijNew*SigmaSquaredMu+SigmaSquared))))
  
  if (Modification == "AD"){ 
    TessStructure=(dbinom(d-1,NumCovariates-1,Omega/NumCovariates))/(dbinom(d-2,NumCovariates-1,Omega/NumCovariates)*(NumCovariates-d+1))
    TransitionRatio=(NumCovariates-d+1)/d;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    if (d==2){
      AcceptenceProb=AcceptenceProb+log(1/2)
    } else if (d==NumCovariates){
      AcceptenceProb=AcceptenceProb+log(2)
    }
  } else if (Modification == "RD"){
    TessStructure=(dbinom(d-1,NumCovariates,Omega/NumCovariates)*(NumCovariates-d))/(dbinom(d,NumCovariates,Omega/NumCovariates))
    TransitionRatio=(d+1)/(NumCovariates-d)
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    if (d==NumCovariates-1){
      AcceptenceProb=AcceptenceProb+log(1/2)
    } else if (d==1){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  } else if(Modification == "AC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar-2,lambda_rate)
    TransitionRatio=1/cStar;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)+0.5*log(SigmaSquared)
    
    if (cStar==2){
      AcceptenceProb=AcceptenceProb+log(1/2);
    }
  } else if (Modification == "RC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar,lambda_rate);
    TransitionRatio=cStar+1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)-0.5*log(SigmaSquared)
    
    if (cStar==1){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  } else if (Modification == "Change"){
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  } else {
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  return(AcceptenceProb)
}

Homo_CalculateResiduals<-function(y,x,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred,TessIndexes){ 
  if (j==1){
    CurrentTessPred<-Pred[[j]][TessIndexes]
    SumOfAllTess=SumOfAllTess-CurrentTessPred
  } else {
    CurrentTessPred<-Pred[[j]][TessIndexes]
    SumOfAllTess=SumOfAllTess+LastTessPred-CurrentTessPred;
  }
  
  IndexesStar=Homo_Indexes(x,TessStar[[j]],DimStar[[j]]);
  R_j<-y-SumOfAllTess
  
  R_ijOld=rep(0,length(Pred[[j]]))
  n_ijOld=rep(0,length(Pred[[j]]))
  
  for (i in 1:length(Pred[[j]])){
    R_ijOld[i]<-sum(R_j[TessIndexes==i])
    n_ijOld[i]<-sum(TessIndexes==i)
  }
  
  R_ijNew=rep(0,length(TessStar[[j]][,1]))
  n_ijNew=rep(0,length(TessStar[[j]][,1]))
  
  for (i in 1:length(TessStar[[j]][,1])){
    R_ijNew[i]<-sum(R_j[IndexesStar==i])
    n_ijNew[i]<-sum(IndexesStar==i)
  }
  
  return(list(R_ijOld,n_ijOld,R_ijNew,n_ijNew,SumOfAllTess,IndexesStar))
}

Homo_NewPredSet<-function(j,Tess,R_ijNew,n_ijNew,sigmaSquaredMu,SigmaSquared){ 
  PredSet=rep(0,length(Tess[[j]][,1]))
  for (i in 1:length(Tess[[j]][,1])){
    PredSet[i]=rnorm(1,sigmaSquaredMu*R_ijNew[i]/(sigmaSquaredMu*n_ijNew[i]+SigmaSquared),((SigmaSquared*sigmaSquaredMu)/(n_ijNew[i]*sigmaSquaredMu+SigmaSquared))^0.5);
  }
  return(PredSet)
}

Homo_TestPrediction<-function(x,m,Tess,Dim,Pred){ 
  Prediction=rep(0,length(x[,1]));
  for (j in 1:m){
    NewTessIndexes=Homo_Indexes(x,Tess[[j]],Dim[[j]]);
    Prediction=Prediction+Pred[[j]][NewTessIndexes]
  }
  return(Prediction)
}