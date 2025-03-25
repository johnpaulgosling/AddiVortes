library('invgamma')

fitting_function<- function(lambda,q,nu,sigmaSquared_hat){ #function that calculates the squared difference between sigma squared hat and the inverse gamma function
  return((sigmaSquared_hat- qinvgamma(q, shape=nu/2, rate=nu*lambda/2))^2)
}
