ImSamplingWithNormal = function(y,X,nu,alpha,N){
  #--Define some real and array constants
  n = nrow(X)
  p = ncol(X)
  xx = t(X)%*%X
  xx.inv = solve(xx)
  xy = t(X)%*%y
  
  #--Preliminary calculations (OLS, sigma2 hat, SSR, residuals,...)
  betaHat = xx.inv%*%xy
  yHat = X%*%betaHat
  uHat = y-yHat
  SCR = sum(uHat^2)
  sigma2Hat = (1/(n-p))*SCR
  varBetaHat = sigma2Hat*xx.inv
  
  #--We simulate the coefficients beta^{(k)} following N{beta_{ols},Var(beta_{ols})} 
  betaSimulated = rmvnorm(n=N, mean=betaHat, sigma=varBetaHat)
  
  #--Function that calculates the weight (here == prior up to a constant)
  weightBridge = function(beta){
    return(exp(-nu*sum((abs(beta))^alpha)))
  }
  
  #--Vector of weights w(beta^{k}) 
  BigWeight = apply(betaSimulated,1,weightBridge)
  
  #--Calculating estimator components "IS"
  ISestimate =  betaSimulated * (BigWeight)
  ISestimate = apply(ISestimate,2,cumsum)*(1/cumsum(BigWeight))
  
  #--Calculation of criterion ESS
  ESS = (sum(BigWeight,na.rm=T))^2/sum(BigWeight^2,na.rm=T)
  
  #--We return the estimator IS, the vector of weights w(beta^{(k)}), ESS criterion, simulated parameters and hyperparameters......
  return(list(IS=ISestimate,weit=BigWeight,ESS=ESS,OLS=betaHat,betaSimulated=betaSimulated,sigma2=sigma2Hat))
}