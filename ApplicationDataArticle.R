

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

# hyperparameters 
alpha = 0.5
nu =  0.467
tau = nu^(-1/alpha)

# Number of simulation
N=4000000

# Go for importance sampling method !! (with hyperparameter nu fixed and sigma2 fixed also )
IS = ImSamplingWithNormal(y,X,nu,alpha,N)
# ESS criterion
(IS$ESS/N)*100


matESS = matrix((IS$ESS/N)*100,nrow=1,dimnames = list(c("ESS"), c("Diabete data")))
stargazer(matESS, summary=FALSE)


# Go for Gibbs sampler method !! (with hyperparameter nu fixed and sigma2 fixed also )
mygrid=sort(union(seq(-800,800,length=1001), seq(-10,10,length=201)))
BayesianBridge = bridge(X, y, 50000, c=1/2, d=1/2, alpha=0.5, mygrid=mygrid, verbose=500, burn=2000, tau.know=tau, sigma2.known = IS$sigma2)

# we keep the weigth of importance sampling and simulation of the beta's
weit = IS$weit
Bayes = IS$IS
betaOLS = IS$OLS
betaSimulated = IS$betaSimulated

#####################
#-------PLOTs-------#
#####################

corrplot.mixed(cor(X),main="Correlation matrix Diabetes Data", mar = c(0,0,2.5,0))
par(mar = c(5.1, 5.2, 4, 2))

pdf(file="poids.pdf")
boxplot(log(weit),main=expression('Logarithme des poids : '*log(w(beta^(k)))))
dev.off()


pdf("diabetesSimul1.pdf",width  = 11)
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


pdf("diabetesSimul2.pdf",width  = 11)
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


pdf("diabetesSimul3.pdf",width  = 11)
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
N=2000000

# Go for importance sampling method !! (with hyperparameter nu fixed and sigma2 fixed also )
IS = ImSamplingWithNormal(y,X,nu,alpha,N)
# ESS criterion
(IS$ESS/N)*100
ESS_Boston = (IS$ESS/N)*100

# Go for Gibbs sampler method !! (with hyperparameter nu fixed and sigma2 fixed also )
mygrid=sort(union(seq(-800,800,length=1001), seq(-10,10,length=201)))
BayesianBridge = bridge(X, y, 50000, c=1/2, d=1/2, alpha=0.5, mygrid=mygrid, verbose=500, burn=2000, tau.know=tau, sigma2.known = IS$sigma2)

# we keep the weigth of importance sampling and simulation of the beta's
weit = IS$weit
Bayes = IS$IS
betaOLS = IS$OLS
betaSimulated = IS$betaSimulated


#####################
#-------PLOTs-------#
#####################

corrplot.mixed(cor(X),main="Matrice de corrélation Boston Data", mar = c(0,0,2.5,0))
par(mar = c(5.1, 5.2, 4, 2))

pdf(file="poids.pdf")
boxplot(log(weit),main=expression('Logarithme des poids : '*log(w(beta^(k)))))
dev.off()


pdf("BostonSimul1.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)
#bhat.bayesmode = rep(0,p)
#for(j in 1:6)

plot(density(BayesianBridge$beta[,1]),col="green",main=colnames(X)[1])
lines(density(betaSimulated[,1],weights=weit/sum(weit)),col='blue',ylim=c(0,0.35))
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




pdf("BostonSimul2.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)
plot(density(BayesianBridge$beta[,5]),col="green",main=colnames(X)[5])
lines(density(betaSimulated[,5],weights=weit/sum(weit)),col='blue',ylim=c(0,0.08))
abline(v=Bayes[N,5], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,5]), col='red', lty='dashed')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,6],weights=weit/sum(weit)),main=colnames(X)[6],col='blue')
lines(density(BayesianBridge$beta[,6]),col="green")
abline(v=Bayes[N,6], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,6]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,7],weights=weit/sum(weit)),main=colnames(X)[7],col='blue')
lines(density(BayesianBridge$beta[,7]),col="green")
abline(v=Bayes[N,7], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,7]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,8],weights=weit/sum(weit)),main=colnames(X)[8],col='blue')
lines(density(BayesianBridge$beta[,8]),col="green")
abline(v=Bayes[N,8], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,8]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))
dev.off()


pdf("BostonSimul3.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)
plot(density(betaSimulated[,9],weights=weit/sum(weit)),main=colnames(X)[9],col='blue')
lines(density(BayesianBridge$beta[,9]),col="green")
abline(v=Bayes[N,9], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,9]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,10],weights=weit/sum(weit)),main=colnames(X)[10],col='blue')
lines(density(BayesianBridge$beta[,10]),col="green")
abline(v=Bayes[N,10], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,10]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))

plot(density(betaSimulated[,11],weights=weit/sum(weit)),main=colnames(X)[11],col='blue')
lines(density(BayesianBridge$beta[,11]),col="green")
abline(v=Bayes[N,11], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,11]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))

plot(density(betaSimulated[,12],weights=weit/sum(weit)),main=colnames(X)[12],col='blue')
lines(density(BayesianBridge$beta[,12]),col="green")
abline(v=Bayes[N,12], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,12]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))

dev.off()


###############################
##-----Ozone Data DATA-------##
###############################


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
dim(X)

###########################
#-------SIMULATIONs-------#
###########################

# hyperparameters 
alpha = 0.5
nu =  0.467
tau = nu^(-1/alpha)

# Number of simulation
N=2000000

# Go for importance sampling method !! (with hyperparameter nu fixed and sigma2 fixed also )
IS = ImSamplingWithNormal(y,X,nu,alpha,N)
# ESS criterion
(IS$ESS/N)*100
ESS_Ozone = (IS$ESS/N)*100

# Go for Gibbs sampler method !! (with hyperparameter nu fixed and sigma2 fixed also )
mygrid=sort(union(seq(-800,800,length=1001), seq(-10,10,length=201)))
BayesianBridge = bridge(X, y, 50000, c=1/2, d=1/2, alpha=0.5, mygrid=mygrid, verbose=500, burn=2000, tau.know=tau, sigma2.known = IS$sigma2)

# we keep the weigth of importance sampling and simulation of the beta's
weit = IS$weit
Bayes = IS$IS
betaOLS = IS$OLS
betaSimulated = IS$betaSimulated


#####################
#-------PLOTs-------#
#####################

corrplot.mixed(cor(X),main="Matrice de corrélation Ozone Data", mar = c(0,0,2.5,0))
par(mar = c(5.1, 5.2, 4, 2))


pdf(file="Ozonepoids.pdf")
boxplot(log(weit),main=expression('Logarithme des poids : '*log(w(beta^(k)))))
dev.off()


pdf("OzoneSimul1.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)
#bhat.bayesmode = rep(0,p)
#for(j in 1:6)

plot(density(BayesianBridge$beta[,1]),col="green",main=colnames(X)[1])
lines(density(betaSimulated[,1],weights=weit/sum(weit)),col='blue')
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


pdf("OzoneSimul2.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)
plot(density(BayesianBridge$beta[,5]),col="green",main=colnames(X)[5])
lines(density(betaSimulated[,5],weights=weit/sum(weit)),col='blue',ylim=c(0,0.08))
abline(v=Bayes[N,5], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,5]), col='red', lty='dashed')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,6],weights=weit/sum(weit)),main=colnames(X)[6],col='blue')
lines(density(BayesianBridge$beta[,6]),col="green")
abline(v=Bayes[N,6], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,6]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,7],weights=weit/sum(weit)),main=colnames(X)[7],col='blue')
lines(density(BayesianBridge$beta[,7]),col="green")
abline(v=Bayes[N,7], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,7]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(betaSimulated[,8],weights=weit/sum(weit)),main=colnames(X)[8],col='blue')
lines(density(BayesianBridge$beta[,8]),col="green")
abline(v=Bayes[N,8], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,8]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))
dev.off()


pdf("OzoneSimul3.pdf",width  = 11)
par(mfrow=c(1,1), mar=c(4,4,2,1), las=1)
plot(density(betaSimulated[,9],weights=weit/sum(weit)),main=colnames(X)[9],col='blue')
lines(density(BayesianBridge$beta[,9]),col="green")
abline(v=Bayes[N,9], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,9]), col='red', lty='dashed')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))
dev.off()

##########################
##-review on stargezer-###
##########################

matESS = matrix(c(ESS_Boston,ESS_Ozone),nrow=1,dimnames = list(c("ESS"), c("Boston","Ozone")))
stargazer(matESS, summary=FALSE)
