###############################################################################
#-----------------------------------------------------------------------------#
#-------------------------MONTE CARLO SEQUENTIEL------------------------------#
#----------------APPLIQUE A LA REGRESSION BRIDGE BAYESIENNE-------------------#
#-----------------------------------------------------------------------------#
###############################################################################


rm(list=ls())
dev.off()


######################
#-----PACKAGES-------#
######################


library(KernSmooth)
## loi normal multivariee0
library(mvtnorm)
## BayesBridge
library(BayesBridge)
library(corrplot)
#library(plotrix)
library(stargazer)
library(zipfR)
library(monomvn)
library(data.table)

#################################
#---IMPORTATION GIBBs SAMPLER---#
#################################

source("BridgeMCMC_TEST.R")

##########################
#-----SMC FUNCTION-------#
##########################

SMC = function(y,X,nu,alpha,N){
  
  #--Definissons quelques constantes reelles et matricielles
  n = nrow(X)
  p = ncol(X)
  xx = t(X)%*%X
  xx.inv = solve(xx)
  xy = t(X)%*%y
  ## liste des temperatures
  powerSequence = seq(0, 1, length.out = 6)
  
  #--Calculs preliminaires (MCO,sigma2,SCR,residus,...)
  betaHat = xx.inv%*%xy
  yHat = X%*%betaHat
  uHat = y-yHat
  SCR = sum(uHat^2)
  sigma2Hat = (1/(n-p))*SCR
  varBetaHat = sigma2Hat*xx.inv
  
  
  #---Fonction qui calcule les poids---#
  weightBridgeSMC = function(beta,power1,power2){
    exp(-(power1-power2)*nu*sum((abs(beta))^alpha))
  }
  
  
  #---Fonction de reechantillonnage---#
  reSample = function(N,weight){
    sample(1:N, size=N, replace=TRUE, prob = weight)
  }
  
  
  #---Distribution pi_{t}---#
  distrib = function(beta,Power){
    exp(-Power*nu*sum((abs(beta))^alpha))
  }
  
  
  #-----MCMC Hasting by Random Walk-------#
  mutation = function(beta,Power,kappa,beta_ALL){
    xi = rmvnorm(n=1, mean = beta, sigma=var(beta_ALL)*kappa)
    u = runif(1)
    rho  = min(distrib(xi,Power)/distrib(beta,Power),1)
    if(u<=rho){
      return(list(newBeta = xi, rho=rho))
    }
    else{
      return(list(newBeta= beta, rho=rho)) 
    }
      
}
  

#---Etape : initialisation---#
betaSample = rmvnorm(n=N, mean=betaHat, sigma=varBetaHat)
ESS = c()
rhoSample = c()
betaSample.list = list()
betaSample.list[[1]] = betaSample


#------iteration------#
for(l in 2:length(powerSequence))
{
  
  message("iteration : ",l)
  ## Calcul des poids
  Weight = apply(betaSample,1,FUN=function(x){weightBridgeSMC(x,powerSequence[l],powerSequence[l-1])})
  WeigthNorm = Weight/sum(Weight)
  
  ## Calcul du critere ESS a chaque iteration
  ESS = c(ESS,((sum(WeigthNorm,na.rm=T))^2/sum(WeigthNorm^2,na.rm=T))*(100/N))
  
  ## ReSampling
  betaSample2 = betaSample[reSample(N,WeigthNorm),]
  betaSample = betaSample2
  rm(betaSample2)
  
  Weight = apply(betaSample,1,FUN=function(x){weightBridgeSMC(x,powerSequence[l],powerSequence[l-1])})
  WeigthNorm = Weight/sum(Weight)
  ESS = c(ESS,((sum(WeigthNorm,na.rm=T))^2/sum(WeigthNorm^2,na.rm=T))*(100/N))
  
  

  for(i in 1:1){
  ## Mutation (MCMC)
  betaALL = betaSample
  MCMC  = apply(betaSample,1,FUN=function(x){mutation(x,powerSequence[l],0.004,betaALL)})
  betaList = lapply(MCMC,FUN=function(x){return(x[[1]])})
  acceptationList = lapply(MCMC,FUN=function(x){return(x[[2]])})
  betaSample = do.call(rbind,betaList)
  rhoSample = c(rhoSample,mean(unlist(acceptationList)))
  }
  
  betaSample.list[[l]] = betaSample
  
}
  
  
#--On retourne les simulations
return(list(betaSample=betaSample,sigma2=sigma2Hat,ESS=ESS,rhoSample=rhoSample,betaSample.list=betaSample.list))
  
}


##########################
#-----DIABETE DATA-------#
##########################


## DOnnees du Diabete : 10 variables explicatives
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

alpha = 0.5
nu =  0.467
tau = nu^(-1/alpha)
#nu = rgamma(1, shape=2, scale=1/2)
#nu = 0.4674039

### Nombre de simulation
N=50000

smc = SMC(y,X,nu,alpha,N)
BetaSimul = smc$betaSample
dim(BetaSimul)

mygrid=sort(union(seq(-800,800,length=1001), seq(-10,10,length=201)))
BayesianBridge = bridge(X,y,50000,c=1/2,d=1/2,alpha=0.5,mygrid=mygrid,verbose=500,burn=2000,tau.known = tau,sigma2.known = smc$sigma2)


#####################
#-------PLOTs-------#
#####################

pdf(file='ESSevolution.pdf')
plot(unlist(smc$ESS),type='l',col='blue',ylab="ESS%")
points(unlist(smc$ESS),pch=19)
title(main="Evolution de l'ESS avant et après Resampling",cex.main=1.3,font.main=6.7)
dev.off()


pdf(file='rateAcceptation.pdf')
plot(unlist(smc$rhoSample),type='l',col='blue',ylab="")
points(unlist(smc$rhoSample),pch=19)
title(main="Taux d'acceptation du MCMC",cex.main=1.3,font.main=6.7)
dev.off()



pdf("diabetesSMC1.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)

plot(density(BetaSimul[,1]),main=colnames(X)[1],col='blue',ylim=c(0,0.35))
lines(density(BayesianBridge$beta[,1]),col="green")
abline(v=mean(BetaSimul[,1]), col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,1]), col='red', lty='dashed')
legend("topright",c("SMC","Gibbs","Gibbs estimator","SMC estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,2]),col="green",main=colnames(X)[2])
lines(density(BetaSimul[,2]),col='blue')
abline(v=mean(BetaSimul[,2]), col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,2]), col='red', lty='dashed')
legend("topleft",c("SMC","Gibbs","Gibbs estimator","SMC estimator"),lwd=1,col=c("blue","green","red","yellow"))



plot(density(BayesianBridge$beta[,3]),col="green",main=colnames(X)[3])
lines(density(BetaSimul[,3]),col='blue')
abline(v=mean(BetaSimul[,3]), col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,3]), col='red', lty='dashed')
legend("topleft",c("SMC","Gibbs","Gibbs estimator","SMC estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,4]),col="green",main=colnames(X)[4])
lines(density(BetaSimul[,4]),col='blue')
abline(v=mean(BetaSimul[,4]), col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,4]), col='red', lty='dashed')
legend("topright",c("SMC","Gibbs","Gibbs estimator","SMC estimator"),lwd=1,col=c("blue","green","red","yellow"))

dev.off()


pdf("diabetesSMC2.pdf",width  = 11)
par(mfrow=c(2,2), mar=c(4,4,2,1), las=1)
plot(density(BayesianBridge$beta[,5]),col="green",main=colnames(X)[5])
lines(density(BetaSimul[,5]),col='blue')
abline(v=mean(BetaSimul[,5]), col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,5]), col='red', lty='dashed')
legend("topright",c("SMC","Gibbs","Gibbs estimator","SMC estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,6]),col="green",main=colnames(X)[6])
lines(density(BetaSimul[,6]),col='blue')
abline(v=mean(BetaSimul[,6]), col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,6]), col='red', lty='dashed')
legend("topright",c("SMC","Gibbs","Gibbs estimator","SMC estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,7]),col="green",main=colnames(X)[7])
lines(density(BetaSimul[,7]),col='blue')
abline(v=mean(BetaSimul[,7]), col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,7]), col='red', lty='dashed')
legend("topright",c("SMC","Gibbs","Gibbs estimator","SMC estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,8]),col="green",main=colnames(X)[8])
lines(density(BetaSimul[,8]),col='blue')
abline(v=mean(BetaSimul[,8]), col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,8]), col='red', lty='dashed')
legend("topright",c("SMC","Gibbs","Gibbs estimator","SMC estimator"),lwd=1,col=c("blue","green","red","yellow"))
dev.off()


pdf("diabetesSMC3.pdf",width  = 11)
par(mfrow=c(1,2), mar=c(4,4,2,1), las=1)
plot(density(BayesianBridge$beta[,9]),col="green",main=colnames(X)[9])
lines(density(BetaSimul[,9]),col='blue')
abline(v=mean(BetaSimul[,9]), col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,9]), col='red', lty='dashed')
legend("topright",c("SMC","Gibbs","Gibbs estimator","SMC estimator"),lwd=1,col=c("blue","green","red","yellow"))



plot(density(BayesianBridge$beta[,10]),col="green",main=colnames(X)[10])
lines(density(BetaSimul[,10]),col='blue')
abline(v=mean(BetaSimul[,10]), col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,10]), col='red', lty='dashed')
legend("topright",c("SMC","Gibbs","Gibbs estimator","SMC estimator"),lwd=1,col=c("blue","green","red","yellow"))


dev.off()