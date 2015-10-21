rm(list=ls())

######################
#-----PACKAGES-------#
######################

library(mvtnorm)
library(BayesBridge)
library(corrplot)
library(stargazer)

##########################################
#-----IMPORTANCE SAMPLING FUNCTION-------#
##########################################

ImSamplingWithHyperPrior = function(y,X,N,c,d){
  #--Define some real and array constants
  n = nrow(X)
  p = ncol(X)
  xx = t(X)%*%X
  xx.inv = solve(xx)
  xy = t(X)%*%y
  yy = t(y)%*%y
  yx = t(y)%*%X
  
  #--Preliminary calculations (OLS, sigma2 hat, SSR, residuals,...)
  betaHat = xx.inv%*%xy
  yHat = X%*%betaHat
  uHat = y-yHat
  SCR = sum(uHat^2)
  sigma2Hat = (1/n)*SCR
  varBetaHat = sigma2Hat*xx.inv
  
  #--Parameters fixed :
  sigma2Hat = (1/n)*SCR #Maximum Likelihood Estimator
  alpha = rep(0.5,N)
  
  #--We simulate the coefficients beta^{(k)} following the importance function N{beta_{ols},Var(beta_{ols})} 
  betaSimulated = rmvnorm(n=N, mean=betaHat, sigma=varBetaHat)
  betaAlpha = cbind(betaSimulated,alpha) #we concatenate beta and alpha in order to simulate hyperparameter nu
  
  #--We simulate the hyperparameter nu^{(k)} with the proposal q(nu|alpha,beta)
  nu = apply(betaAlpha,1,FUN=function(x){rgamma(1,p/x[p+1]+c,rate=d+sum(abs(x[1:p])^x[p+1]))})
  
  #--Function that calculates the importance weights !! (there are simplications behind)
  weightBridge = function(beta,alpha){
    return((d+sum((abs(beta))^alpha))^(-p/alpha-c))
  }
    
  #--Vector of weights w(beta^{k},alpha) 
  BigWeight = apply(betaAlpha,1,FUN=function(x){weightBridge(x[1:p],x[p+1])})
  
  #--Calculating estimator components "IAS"
  ISestimate = betaAlpha[,1:p] * (BigWeight)
  ISestimate = apply(ISestimate,2,cumsum)*(1/cumsum(BigWeight))
  
  #--Calculation of criterion ESS
  ESS = (sum(BigWeight,na.rm=T))^2/sum(BigWeight^2,na.rm=T)
  
  #--We give the estimator "IAS", the vector of weights w(beta^{(k)},alpha), ESS criterion, simulated parameters and hyperparameters...
  return(list(IS=ISestimate, weit=BigWeight, ESS=ESS, OLS=betaHat, betaSimulated=betaSimulated, sigma2=sigma2Hat, nu=nu, alpha=alpha))
  
}

##########################
#-----DIABETE DATA-------#
##########################

## 10 variables 
data(diabetes, package="lars")
X = diabetes$x
y = diabetes$y
y = (y - mean(y))
p = ncol(X)
n = length(y)
for(j in 1:p)
{
  X[,j] = (X[,j] - mean(X[,j]))/sd(X[,j])
}


###########################
#-------SIMULATIONs-------#
###########################

# Number of simulation
N=4000000

# Importance Sampling
IS = ImSamplingWithHyperPrior(y,X,N,2,2)
(IS$ESS/N)*100
ESS_DIABETE_DATA = (IS$ESS/N)*100

# Bayesian Bridge with Gibbs Sampler
BayesianBridge = bridge.reg(y, X, alpha=0.5,nsamp=N, nu.shape=2, nu.rate = 2 ,sig2.true=IS$sigma2 , method="stable",burn=1000)

# we keep the weigth of importance sampling and simulation of the beta's
weit = IS$weit
Bayes = IS$IS
betaOLS = IS$OLS
betaSimulated = IS$betaSimulated


#####################
#-------PLOTs-------#
#####################

pdf(file="poids.pdf")
boxplot(log(weit),main=expression('Logarithme des poids : '*log(w(beta^(k)))))
dev.off()


pdf("diabetesSimulHyperPrior1.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)
#bhat.bayesmode = rep(0,p)
#for(j in 1:6)

plot(density(betaSimulated[,1],weights=weit/sum(weit)),main=colnames(X)[1],col='blue',ylim=c(0,0.35))
lines(density(BayesianBridge$beta[,1]),col="green")
abline(v=Bayes[N,1], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,1]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))



plot(density(betaSimulated[,2],weights=weit/sum(weit)),main=colnames(X)[2],col='blue')
lines(density(BayesianBridge$beta[,2]),col="green")
abline(v=Bayes[N,2], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,2]), col='red', lty='dashed')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,3],weights=weit/sum(weit)),main=colnames(X)[3],col='blue')
lines(density(BayesianBridge$beta[,3]),col="green")
abline(v=Bayes[N,3], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,3]), col='red', lty='dashed')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,4],weights=weit/sum(weit)),main=colnames(X)[4],col='blue')
lines(density(BayesianBridge$beta[,4]),col="green")
abline(v=Bayes[N,4], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,4]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))

dev.off()



pdf("diabetesSimulHyperPrior2.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)

plot(density(BayesianBridge$beta[,5]),col="green",main=colnames(X)[5])
lines(density(betaSimulated[,5],weights=weit/sum(weit)),col='blue',ylim=c(0,0.08))
abline(v=Bayes[N,5], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,5]), col='red', lty='dashed')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,6],weights=weit/sum(weit)),main=colnames(X)[6],col='blue',ylim=c(0,0.15))
lines(density(BayesianBridge$beta[,6]),col="green")
abline(v=Bayes[N,6], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,6]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,7],weights=weit/sum(weit)),main=colnames(X)[7],col='blue',ylim=c(0,0.08))
lines(density(BayesianBridge$beta[,7]),col="green")
abline(v=Bayes[N,7], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,7]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,8],weights=weit/sum(weit)),main=colnames(X)[8],col='blue',ylim=c(0,0.15))
lines(density(BayesianBridge$beta[,8]),col="green")
abline(v=Bayes[N,8], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,8]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))

dev.off()


pdf("diabetesSimulHyperPrior3.pdf",width  = 11)
par(mfrow=c(1,2), mar=c(4,4,2,1), las=1)


plot(density(betaSimulated[,9],weights=weit/sum(weit)),main=colnames(X)[9],col='blue',ylim=c(0,0.1))
lines(density(BayesianBridge$beta[,9]),col="green")
abline(v=Bayes[N,9], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,9]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,10],weights=weit/sum(weit)),main=colnames(X)[10],col='blue',ylim=c(0,0.19))
lines(density(BayesianBridge$beta[,10]),col="green")
abline(v=Bayes[N,10], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,10]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))

dev.off()


##############################
#-----Ozone Data DATA-------##
##############################


data("Ozone", package="mlbench")
summary(Ozone)
class(Ozone)
X = Ozone[,5:13]
colnames(Ozone)
dim(X)
## Plot des corrélations des 
#
head(X)

y = Ozone[,4]
y = y-mean(y,na.rm=T)
## on remplace les Na's
for(j in 1:ncol(X)){
  X[,j][which(is.na(X[,j])==T)] = mean(X[,j],na.rm=T)
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
class(X)
X = as.matrix(X)


###########################
#-------SIMULATIONs-------#
###########################

# Number of simulation
N=4000000

# Importance Sampling
IS = ImSamplingWithHyperPrior(y,X,N,2,2)
(IS$ESS/N)*100
ESS_OZONE_DATA = (IS$ESS/N)*100
# Bayesian Bridge with Gibbs Sampler
BayesianBridge = bridge.reg(y, X, alpha=0.5,nsamp=N, nu.shape=2, nu.rate = 2 ,sig2.true=IS$sigma2 , method="stable",burn=1000)

# we keep the weigth of importance sampling and simulation of the beta's
weit = IS$weit
Bayes = IS$IS
betaOLS = IS$OLS
betaSimulated = IS$betaSimulated


#####################
#-------PLOTs-------#
#####################

pdf(file="poids.pdf")
boxplot(log(weit),main=expression('Logarithme des poids : '*log(w(beta^(k)))))
dev.off()


pdf("ozoneSimulHyperPrior1.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)
#bhat.bayesmode = rep(0,p)
#for(j in 1:6)

plot(density(betaSimulated[,1],weights=weit/sum(weit)),main=colnames(X)[1],col='blue',ylim=c(0,0.35))
lines(density(BayesianBridge$beta[,1]),col="green")
abline(v=Bayes[N,1], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,1]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))



plot(density(betaSimulated[,2],weights=weit/sum(weit)),main=colnames(X)[2],col='blue')
lines(density(BayesianBridge$beta[,2]),col="green")
abline(v=Bayes[N,2], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,2]), col='red', lty='dashed')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,3],weights=weit/sum(weit)),main=colnames(X)[3],col='blue')
lines(density(BayesianBridge$beta[,3]),col="green")
abline(v=Bayes[N,3], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,3]), col='red', lty='dashed')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,4],weights=weit/sum(weit)),main=colnames(X)[4],col='blue')
lines(density(BayesianBridge$beta[,4]),col="green")
abline(v=Bayes[N,4], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,4]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))

dev.off()



pdf("ozoneSimulHyperPrior2.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)

plot(density(BayesianBridge$beta[,5]),col="green",main=colnames(X)[5])
lines(density(betaSimulated[,5],weights=weit/sum(weit)),col='blue',ylim=c(0,0.08))
abline(v=Bayes[N,5], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,5]), col='red', lty='dashed')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,6],weights=weit/sum(weit)),main=colnames(X)[6],col='blue',ylim=c(0,0.15))
lines(density(BayesianBridge$beta[,6]),col="green")
abline(v=Bayes[N,6], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,6]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,7],weights=weit/sum(weit)),main=colnames(X)[7],col='blue',ylim=c(0,0.08))
lines(density(BayesianBridge$beta[,7]),col="green")
abline(v=Bayes[N,7], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,7]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,8],weights=weit/sum(weit)),main=colnames(X)[8],col='blue',ylim=c(0,0.15))
lines(density(BayesianBridge$beta[,8]),col="green")
abline(v=Bayes[N,8], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,8]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))

dev.off()


pdf("ozoneSimulHyperPrior3.pdf",width  = 11)
par(mfrow=c(1,2), mar=c(4,4,2,1), las=1)


plot(density(betaSimulated[,9],weights=weit/sum(weit)),main=colnames(X)[9],col='blue',ylim=c(0,0.1))
lines(density(BayesianBridge$beta[,9]),col="green")
abline(v=Bayes[N,9], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,9]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,10],weights=weit/sum(weit)),main=colnames(X)[10],col='blue',ylim=c(0,0.19))
lines(density(BayesianBridge$beta[,10]),col="green")
abline(v=Bayes[N,10], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,10]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))

dev.off()

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
X = X[,-1]

## On centre et on reduit les 'X'
y = (y - mean(y))
mX = colMeans(X)
for(i in 1:n){ X[i,] = X[i,] - mX }
X = X[,-1]

for(j in 1:ncol(X)){
  
  X[,j] = X[,j]/sd(X[,j])
  
}

###########################
#-------SIMULATIONs-------#
###########################

# Number of simulations
N=2000000

# Importance Sampling
IS = ImSamplingWithHyperPrior(y,X,N,2,2)
(IS$ESS/N)*100
ESS_BOSTON_DATA = (IS$ESS/N)*100


# Bayesian Bridge with Gibbs Sampler
BayesianBridge = bridge.reg(y, X, alpha=0.5,nsamp=5000000, nu.shape=0.5, nu.rate = 2 ,sig2.true=IS$sigma2 , method="stable",burn=1000)

# we keep the weigth of importance sampling and simulation of the beta's
weit = IS$weit
Bayes = IS$IS
betaOLS = IS$OLS
betaSimulated = IS$betaSimulated
par(mar = c(5.1, 5.2, 4, 2))


#####################
#-------PLOTs-------#
#####################


##########################
##-review on stargezer-###
##########################

matESS = matrix(c(ESS_DIABETE_DATA,ESS_OZONE_DATA,ESS_BOSTON_DATA),nrow=1,dimnames = list(c("ESS"), c("Diabete", "Ozone","Boston")))

stargazer(matESS, summary=FALSE)
