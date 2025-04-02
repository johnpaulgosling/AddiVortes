#' @title AddiVortes
#'
#' @description
#' The AddiVortes function is a Bayesian nonparametric regression model that uses a tessellation to model the relationship between the covariates and the output values. The model uses a backfitting algorithm to sample from the posterior distribution of the output values for each tessellation. The function returns the RMSE value for the test samples.
#'
#'
#' @param y A vector of the output values.
#' @param x A matrix of the covariates.
#' @param m The number of tessellations.
#' @param max_iter The number of iterations.
#' @param burn_in The number of burn in iterations.
#' @param nu The degrees of freedom.
#' @param q The quantile.
#' @param k The number of centers.
#' @param sd The standard deviation.
#' @param Omega The prior probability of adding a dimension.
#' @param lambda_rate The rate of the Poisson distribution for the number of centers.
#' @param YTest A vector of the output values for the test set.
#' @param XTest A matrix of the covariates for the test set.
#' @param IntialSigma The method used to calculate the initial sigma.
#' @param thinning The thinning rate.
#'
#' @return The RMSE value for the test samples.
#'
#' @export
AddiVortes<-function(y,x,m = 200 ,max_iter = 1200,burn_in= 200,nu = 6,q =0.85,k = 3 ,sd = 0.8 ,Omega = 3,lambda_rate = 25,YTest,XTest,IntialSigma = "Linear",thinning=1){

  #Scaling x and y
  yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
  xScaled=x;
  for (i in 1:length(x[1,])){
    xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }


  for (i in 1:length(XTest[1,])){
    XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }

  #Initialize:
  #Prediction Set (A list of vectors with the output values for each tessellation),
  #Dimension set (A list of vectors with the covariates included in the tessellaions);
  #and Tessellation Set (A list of matrices that give the coordinates of the centers in the tessellations)

  Pred<-rep(list(matrix(mean(yScaled)/m)),m)
  Dim=vector(length = m)
  Tess=vector(length = m)
  for (i in 1:m){
    Dim[i]<-list(sample(1:length(x[1,]), 1))
    Tess[i]<-(list(matrix(rnorm(1,0,sd))))
  }

  #Prepare some variables used in the backfitting algorithm
  SumOfAllTess=rep(mean(yScaled),length(yScaled))
  SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
  LastTessPred=matrix

  #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
  posterior_samples<-floor((max_iter - burn_in) / thinning)+1
  PredictionMatrix<-array(dim=c(length(y),posterior_samples))
  TestMatrix<-array(dim=c(length(YTest),posterior_samples))

  #finding lambda
  if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
    SigmaSquaredHat=var(yScaled)
  }
  else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
    MultiLinear<-lm(yScaled ~ xScaled)
    SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
  }

  #Find lambda
  lambda=1;
  lambda <- optim(par = 1,
                  fitting_function,
                  method = "Brent",
                  lower = 0.001,
                  upper = 100,
                  q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par

  for (i in 1:max_iter){

    #Sample Sigma squared using all tessellations to predict the outcome variables
    SigmaSquared=Sample_Sigma_Squared(yScaled,nu,lambda,SumOfAllTess)

    for (j in 1:m){
      NewTessOutput<-Propose_Tessellation(xScaled,j,Tess,Dim,sd) #Propose new Tessellation
      TessStar<-NewTessOutput[[1]]
      DimStar<-NewTessOutput[[2]]
      Modification<-NewTessOutput[[3]]

      ResidualsOutput<-Calculate_Residuals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
      R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
      n_ijOld<-ResidualsOutput[[2]]
      R_ijNew<-ResidualsOutput[[3]]
      n_ijNew<-ResidualsOutput[[4]]
      SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
      IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
      Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.

      if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.

        LOGAcceptanceProb=Acceptance_Probability(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptance probability.

        if (log(runif(n=1, min=0, max=1))<LOGAcceptanceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation.
          Tess=TessStar
          Dim=DimStar
          Pred[[j]]=Sample_mu_values(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
          LastTessPred=Pred[[j]][IndexesStar]
        }
        else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=Sample_mu_values(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
      }
      else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
        Pred[[j]]=Sample_mu_values(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
        LastTessPred=Pred[[j]][Indexes];
      }
      if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
        SumOfAllTess=SumOfAllTess+LastTessPred;
      }
    }

    if (i %% 100 == 0){
      cat(sprintf("Iteration %d out of %d", i, max_iter), "\n")
    }

    if (i>=burn_in & (i-burn_in) %% thinning ==0 ){ #vectors that hold the predictions for each iteration after burn in.
      PredictionMatrix[,1+(i-burn_in)/thinning]=SumOfAllTess;
      TestMatrix[,1+(i-burn_in)/thinning]=Prediction_for_TestSet(XTest,m,Tess,Dim,Pred);
    }
  }

  #finding the mean of the predition over the iterations and then unscaling the predictions.
  mean_yhat=(rowSums(PredictionMatrix)/(posterior_samples))*(max(y)-min(y))+((max(y)+min(y))/2)
  mean_yhat_Test=(rowSums(TestMatrix)/(posterior_samples))*(max(y)-min(y))+((max(y)+min(y))/2)

  return( #Returns the RMSE value for the test samples.
    data.frame(
      In_sample_RMSE = sqrt(mean((y-mean_yhat)^2)),
      Out_of_sample_RMSE = sqrt(mean((YTest-mean_yhat_Test)^2))
    )
  )
}
