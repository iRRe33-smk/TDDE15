# Matern32  kernel
k <- function(sigmaf = 1, ell = 1)  
{   
  rval <- function(x, y = NULL) 
  {	r = sqrt(crossprod(x-y))
  return(sigmaf^2*(1+sqrt(3)*r/ell)*exp(-sqrt(3)*r/ell))   
  }   
  class(rval) <- "kernel"   
  return(rval) 
} 

posteriorGP = function(X, y, X_star, sigmaNoise, k, ...){
  n = length(X)
  
  K = k(X,X, ... )
  
  k_star = k(X,X_star, ... )
  
  L = t(chol(K + sigmaNoise^2*diag(n)))
  
  alpha = solve(t(L),solve(L,y))
  
  mean_pred = t(k_star) %*% alpha #predictive mean
  
  v = solve(L,k_star)
  
  var_pred= k(X_star, X_star) - t(v) %*% v # predictive var
  
  logmarglik = -.5 * t(y) %*% alpha - sum(log(diag(L))) - .5*n*log(2*pi)
  
  return(list(mean=mean_pred, var= var_pred, logmarglik = logmarglik))
  
}