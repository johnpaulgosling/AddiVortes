Sample_mu_values<-function(j,Tess,R_ijNew,n_ijNew,sigmaSquaredMu,SigmaSquared){ #Sampling the new output values for the new tessellation.
  PredSet=rep(0,length(Tess[[j]][,1]))
  for (i in 1:length(Tess[[j]][,1])){
    PredSet[i]=rnorm(1,sigmaSquaredMu*R_ijNew[i]/(sigmaSquaredMu*n_ijNew[i]+SigmaSquared),((SigmaSquared*sigmaSquaredMu)/(n_ijNew[i]*sigmaSquaredMu+SigmaSquared))^0.5);
  }
  return(PredSet)
}
