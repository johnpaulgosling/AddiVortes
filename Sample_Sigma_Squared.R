library('invgamma')

Sample_Sigma_Squared<-function(yScaled,nu,lambda,SumOfAllTess){ #Sample sigma squared from inverse gamma distribution
  
  n=length(yScaled)
  SigmaSquared<-rinvgamma(1,shape=(nu+n)/2,rate=(nu*lambda+sum((yScaled-SumOfAllTess)^2))/2)
  
  return(SigmaSquared)
}
