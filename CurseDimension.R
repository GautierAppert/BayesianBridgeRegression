rm(list=ls())

######################
#-----PACKAGES-------#
######################

library(mvtnorm)
library(BayesBridge)
library(corrplot)

# import functions for the Gibbs Sampler and Importance Sampling methods
source("BridgeMCMC_TEST.R")

##########################################
#-----IMPORTANCE SAMPLING FUNCTION-------#
##########################################

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


################################
#-----BostonHousing DATA-------#
################################

data("BostonHousing", package="mlbench")

## Setup
y = BostonHousing$medv
n = length(y)

## No interactions or squared terms.
X = model.matrix(medv ~ ., data=BostonHousing);
cnames = colnames(X);
dim(X)

head(X)
data.frame(X)$chas1

## Include interations and squared terms
idc  = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13);
nidc = length(idc)

X.base  = model.matrix(medv ~ . * ., data=BostonHousing);

X = X.base

cnames = colnames(X.base);
for (i in 1:nidc){
  nm = paste(colnames(BostonHousing)[idc[i]], 2, sep="");
  cnames = c(cnames, nm);
  X = cbind(X, BostonHousing[,idc[i]] * BostonHousing[,idc[i]]);
}

## On centre et on reduit les 'X'
y = (y - mean(y))
mX = colMeans(X)
for(i in 1:n){ X[i,] = X[i,] - mX }
X = X[,-1]

for(j in 1:ncol(X)){
  
  X[,j] = X[,j]/sd(X[,j])
  
}

cnames = colnames(X);

head(X)
class(X)
dim(X)


###########################
#-------SIMULATIONs-------#
###########################

# hyperparameters 
alpha = 0.5
nu =  0.467
tau = nu^(-1/alpha)

# Number of simulation
N=1000000

# Go for importance sampling method !! (with hyperparameter nu fixed and sigma2 fixed also )
IS = ImSamplingWithNormal(y,X,nu,alpha,N)
# ESS criterion
(IS$ESS/N)*100
ESS_BOSTON_DATA = (IS$ESS/N)*100

# Go for Gibbs sampler method !! (with hyperparameter nu fixed and sigma2 fixed also )
# mygrid=sort(union(seq(-800,800,length=1001), seq(-10,10,length=201)))
# BayesianBridge = bridge(X, y, 50000, c=1/2, d=1/2, alpha=0.5, mygrid=mygrid, verbose=500, burn=2000, tau.know=tau, sigma2.known = IS$sigma2)

# we keep the weigth of importance sampling and simulation of the beta's
weit = IS$weit
Bayes = IS$IS
betaOLS = IS$OLS
betaSimulated = IS$betaSimulated


#############################
#-----Ozone Data DATA-------#
#############################

data("Ozone", package="mlbench")
summary(Ozone)
class(Ozone)
X = Ozone[,5:13]
dim(X)
## Plot des corrélations des 
# corrplot(cor(X,use = "complete.obs"),order="hclust",addCoef.col="grey")
# head(X)


y = Ozone[,4]
y = y-mean(y,na.rm=T)

## on remplace les Na's
for(j in 1:ncol(X)){
  
  X[,j][which(is.na(X[,j])==T)] = mean(X[,j],na.rm=T)
  
}


## On créer les termes d'interaction croisé et les termes au carrés : 
for(j in 1:9){
  for(c in j:9){
    
    X = cbind(X,X[,j]*X[,c]) 
  }
}


# on centre et on reduit :
for(j in 1:ncol(X)){
  
  X[,j] =  (X[,j]-mean(X[,j]))/sd(X[,j])
  
}


dim(X)

indiceNA = which(is.na(y)==T)
y = y[-indiceNA]
X = X[-indiceNA,]
dim(X)
length(y)


###########################
#-------SIMULATIONs-------#
###########################

# hyperparameters 
alpha = 0.5
nu =  0.467
tau = nu^(-1/alpha)

# Number of simulation
N=500000

# Go for importance sampling method !! (with hyperparameter nu fixed and sigma2 fixed also )
IS = ImSamplingWithNormal(y,X,nu,alpha,N)
# ESS criterion
(IS$ESS/N)*100
ESS_OZONE_DATA = (IS$ESS/N)*100



# Go for Gibbs sampler method !! (with hyperparameter nu fixed and sigma2 fixed also )
# mygrid=sort(union(seq(-800,800,length=1001), seq(-10,10,length=201)))
# BayesianBridge = bridge(X, y, 50000, c=1/2, d=1/2, alpha=0.5, mygrid=mygrid, verbose=500, burn=2000, tau.know=tau, sigma2.known = IS$sigma2)

# we keep the weigth of importance sampling and simulation of the beta's
weit = IS$weit
Bayes = IS$IS
betaOLS = IS$OLS
betaSimulated = IS$betaSimulated


##########################
##-review on stargezer-###
##########################

matESS = matrix(c(ESS_BOSTON_DATA,ESS_OZONE_DATA),nrow=1,dimnames = list(c("ESS"), c("Boston","Ozone")))
stargazer(matESS, summary=FALSE)

