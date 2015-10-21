####################################################################
#------------------------------------------------------------------#
#---Comparaison des methodes OLS, Bridge(EM), Bayes Bridge---------#
#----------a l'aide de critere de prediction-----------------------#
#------------------------------------------------------------------#
####################################################################


rm(list=ls())
dev.off()

#-------------------------#
#--Importation packages---#
#-------------------------#


# install.packages("KernSmooth")
# install.packages('BayesBrigde')
# install.packages('mlbench')
# install.packages('chemometrics')
# install.packages('corrplot')
# install.packages('pls')
# install.packages('caret')
# install.packages('lars)

library(KernSmooth)
library(BayesBridge)
library(mlbench)
library(chemometrics)
library(corrplot)
library(pls)
## package pour la prediction : echantillon apprentissage et test :
library(caret)
## package pour la regression Ridge :  fonction "lm.ridge()"
library(MASS)
## package pour la regression LASSSO
library(lars)


###########################################################################################
#----------------------------------------------------------------------------------------##
#------------CACUL DES ERREURS DE PREDICTION TYPE RMSE SUR TROIS JEU DE DONNEES----------##
#----------------------------------------------------------------------------------------##
###########################################################################################

################################
#------------------------------#
#-----BostonHousing DATA-------#
#------------------------------#
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

## Plot des corrélations des 
corrplot.mixed(cor(X[,-1]),order="hclust",addCoef.col="grey")

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


# On stocke les MSE pour les différents modeles
List.Linear = list()
List.ClassicalBridge = list()
List.BayesianBridgeFix = list()
List.BayesianBridge = list()
List.Ridge = list()
List.Lasso = list()

for(i in 1:30){
  
  message(i)
  ## on créée tout d'abord un échantillon d'apprentissagge et un échantillon test 
  indiceTrain = createDataPartition(y, p = 0.87, list = FALSE)
  trainX = X[indiceTrain, ]
  trainY = y[indiceTrain]
  
  ## echantillon Test
  testX = X[-indiceTrain, ]
  testY = y[-indiceTrain]
  
  ## Estimation
  # burn = 500 par default
  BayesianBridgeFix = bridge.reg(trainY, trainX, nsamp=30000,alpha=0.5, sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0 , method="stable")
  BayesianBridge = bridge.reg(trainY, trainX, nsamp=30000,sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,method="stable")
  ClassicalBridge = bridge.EM(trainY, trainX, 0.5, 1, 1e8, 1e-9, 30,use.cg=TRUE)
  
  #bridge.EM(y, X, 0.5, 1.0, 1e8, 1e-9, 30, use.cg=TRUE);
  
  
  Ols = lm(trainY~trainX-1)
  Ridge = lm.ridge(trainY~trainX-1,lambda=seq(0,2,0.01))
  Ridge = lm.ridge(trainY~trainX-1,lambda=  as.numeric(names(which.min(Ridge$GCV))))
  #Lasso = lars(trainX,trainY, type = c("lasso"),normalize = FALSE, intercept = FALSE)
  frac.delta = seq(from = 0, to = 1, 0.01)
  mse.cv.Lasso = cv.lars(trainX,trainY,trace=F,plot.it = F,intercept=F,normalize=F, index = frac.delta, type = c("lasso"))
  lambdaOpt = frac.delta[which.min(mse.cv.Lasso$cv)]  
  # plot(frac.delta,  mse.cv.Lasso$cv,xlab="delta",ylab="MSEP",type="l",lwd=2)
  # abline(v=lambdaOpt,col="red")
  Lasso = lars(trainX, trainY, type = c("lasso"), use.Gram=F, intercept=F, normalize=F)


  # Prediction
  predLinear = testX%*%Ols$coefficients
  predBridge = testX%*%ClassicalBridge
  predBayesianBridgeFix = testX%*%apply(BayesianBridgeFix$beta,2,mean)
  predBayesianBridge = testX%*%apply(BayesianBridge$beta,2,mean)
  predRidge = testX%*%Ridge$coef
  predLasso = predict(Lasso, testX, s=lambdaOpt , mode="fraction")
  
  SSE.Linear=sum((testY-predLinear)^2)
  SSE.ClassicalBridge = sum((testY-predBridge)^2)
  SSE.BayesianBridgeFix = sum((testY-predBayesianBridgeFix)^2)
  SSE.BayesianBridge = sum((testY-predBayesianBridge)^2)
  SSE.Ridge = sum((testY-predRidge)^2)
  SSE.Lasso = sum((testY-predLasso$fit)^2)
  
  List.Linear[[i]] =  SSE.Linear
  List.ClassicalBridge[[i]] =  SSE.ClassicalBridge
  List.BayesianBridgeFix[[i]] =   SSE.BayesianBridgeFix
  List.BayesianBridge[[i]] =   SSE.BayesianBridge
  List.Ridge[[i]] = SSE.Ridge
  List.Lasso[[i]] = SSE.Lasso
  
}

#png("/home/gappert/Bridge/boston.png",height=662, width=935)
pdf("bostonPred.pdf",width = 11)
boxplot(unlist(List.Linear),unlist(List.ClassicalBridge),unlist(List.BayesianBridgeFix),unlist(List.BayesianBridge),unlist(List.Ridge),unlist(List.Lasso),
        names=c("OLS","Bridge(EM+GCV)","Bayes(a=0.5)","Bayes","Ridge(GCV)","Lasso(GCV)"),
        col=c("red","yellow","blue","green","violet","orange"),ylab="SSE",ylim = c(0, 3000))
title(main="Boston Housing data",cex.main = 1.5, font.main= 6.7, col.main= "black")
dev.off()

#dev.off()


#
## ACF sur les Beta's du modele Bridge Bayesien
# BayesianBridge = bridge.reg(y, X, nsamp=5000,sig2.shape=100, sig2.scale=100, nu.shape=2.0, nu.rate=2.0,method="stable")
# apply(BayesianBridge$beta,2,acf)
# 
# dev.off()
# for(i in 1:5){
#   
#   acf(BayesianBridge$beta[,i],main=colnames(BayesianBridge$beta)[i])
#   
# }




#############################
#---------------------------#
#-----Ozone Data DATA-------#
#---------------------------#
#############################

data("Ozone", package="mlbench")
summary(Ozone)
class(Ozone)
X = Ozone[,5:13]
dim(X)
## Plot des corrélations des 
corrplot(cor(X,use = "complete.obs"),order="hclust",addCoef.col="grey")
head(X)


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


#  On stocke les MSE pour les différents modeles
List.Linear = list()
List.ClassicalBridge = list()
List.BayesianBridgeFix = list()
List.BayesianBridge = list()
List.Ridge = list()
List.Lasso = list()

for(i in 1:20){
  
  message(i)
  ## on créée tout d'abord un échantillon d'apprentissagge et un échantillon test 
  indiceTrain = createDataPartition(y, p = 0.87, list = FALSE)
  trainX = X[indiceTrain, ]
  trainX = apply(trainX,2,as.double)
  trainY = y[indiceTrain]
  
  testX = X[-indiceTrain, ]
  testX = apply(testX,2,as.double)
  
  testY = y[-indiceTrain]
  # Estimation
  BayesianBridgeFix = bridge.reg(trainY, trainX, nsamp=30000,alpha=0.5,sig2.shape=0, sig2.scale=0, nu.shape=2.0, nu.rate=2.0,method="stable")
  BayesianBridge = bridge.reg(trainY, trainX, nsamp=30000,sig2.shape=0, sig2.scale=0, nu.shape=2.0, nu.rate=2.0,method="stable")
  ClassicalBridge = bridge.EM(trainY, trainX, 0.5, 1, 1e8, 1e-9, 30)
  Ols = lm(trainY~trainX-1)
  Ridge = lm.ridge(trainY~trainX-1,lambda = seq(0,2,0.01))
  Ridge = lm.ridge(trainY~trainX-1,lambda =  as.numeric(names(which.min(Ridge$GCV))))
  #Lasso = lars(trainX,trainY, type = c("lasso"),normalize = FALSE, intercept = FALSE)
  frac.delta = seq(from = 0, to = 1, 0.01)
  mse.cv.Lasso = cv.lars(trainX,trainY,trace=F,plot.it = F,intercept=F,normalize=F, index = frac.delta, type = c("lasso"))
  lambdaOpt = frac.delta[which.min(mse.cv.Lasso$cv)]  
  #   plot(frac.delta,  mse.cv.Lasso$cv,xlab="delta",ylab="MSEP",type="l",lwd=2)
  #   abline(v=lambdaOpt,col="red")
  Lasso = lars(trainX, trainY, type = c("lasso"), use.Gram=F, intercept=F, normalize=F)
  
  
  # Prediction
  predLinear = testX%*%Ols$coefficients
  predBridge = testX%*%ClassicalBridge
  predBayesianBridgeFix = testX%*%apply(BayesianBridgeFix$beta,2,mean)
  predBayesianBridge = testX%*%apply(BayesianBridge$beta,2,mean)
  predRidge = testX%*%Ridge$coef
  predLasso = predict(Lasso, testX, s=lambdaOpt , mode="fraction")
  
  SSE.Linear=sum((testY-predLinear)^2)
  SSE.ClassicalBridge = sum((testY-predBridge)^2)
  SSE.BayesianBridgeFix = sum((testY-predBayesianBridgeFix)^2)
  SSE.BayesianBridge = sum((testY-predBayesianBridge)^2)
  SSE.Ridge = sum((testY-predRidge)^2)
  SSE.Lasso = sum((testY-predLasso$fit)^2)
  
  List.Linear[[i]] =  SSE.Linear
  List.ClassicalBridge[[i]] =  SSE.ClassicalBridge
  List.BayesianBridgeFix[[i]] =   SSE.BayesianBridgeFix
  List.BayesianBridge[[i]] =   SSE.BayesianBridge
  List.Ridge[[i]] = SSE.Ridge
  List.Lasso[[i]] = SSE.Lasso
  
}

#png("/home/gappert/Bridge/ozone.png",height=662, width=935)
pdf("ozonePred.pdf",width = 11)
boxplot(unlist(List.Linear),unlist(List.ClassicalBridge),unlist(List.BayesianBridgeFix),unlist(List.BayesianBridge),unlist(List.Ridge),unlist(List.Lasso),
        names=c("OLS","Bridge(EM+GCV)","Bayes(a=0.5)","Bayes","Ridge(GCV)","Lasso(GCV)"),
        col=c("red","yellow","blue","green","violet","orange"),ylab="SSE",ylim = c(0, 1500))
title(main="Ozone data",cex.main = 1.5, font.main= 6.7, col.main= "black")
dev.off()
#dev.off()


#############################
#---------------------------#
#-----Glucose Data----------#
#---------------------------#
#############################

data(NIR, package="chemometrics")
str(NIR)
head(NIR$x)
dim(NIR$x)
dim(NIR$y)


X = NIR$x
y = NIR$y$Glucose

# ## Plot des corrélations des 
# corrplot(cor(X,use = "complete.obs"),order="hclust",addCoef.col="grey")
# head(X)

# On centre
y = y-mean(y)
length(y)
dim(X)
class(X)


# on centre (et on reduit) :
for(j in 1:ncol(X)){
  
  X[,j] =  (X[,j]-mean(X[,j]))/sd(X[,j])
  
}




#  On stocke les MSE pour les différents modeles
List.Linear = list()
List.ClassicalBridge = list()
List.BayesianBridgeFix = list()
List.BayesianBridge = list()
List.Ridge = list()
List.Lasso = list()
List.PLS = list()

for(i in 1:10){
  
  message(i)
  ## on créée tout d'abord un échantillon d'apprentissagge et un échantillon test 
  indiceTrain = createDataPartition(y, p = 0.87, list = FALSE)
  trainX = X[indiceTrain, ]
  trainX = apply(trainX,2,as.double)
  trainY = y[indiceTrain]
  
  testX = X[-indiceTrain, ]
  testX = apply(testX,2,as.double)
  
  testY = y[-indiceTrain]
  
  # Estimation
  BayesianBridgeFix = bridge.reg(trainY, trainX, nsamp=5000,alpha=0.5,sig2.shape=100, sig2.scale=100, nu.shape=2.0, nu.rate=2.0,method="stable")
  BayesianBridge = bridge.reg(trainY, trainX, nsamp=5000,sig2.shape=100, sig2.scale=100, nu.shape=2.0, nu.rate=2.0,method="stable")
  ClassicalBridge = bridge.EM(trainY, trainX, 0.5, 1, 1e8, 1e-9, 30)
  Ols = lm(trainY~trainX-1)
  Ridge = lm.ridge(trainY~trainX-1,lambda=seq(0,2,0.01))
  Ridge = lm.ridge(trainY~trainX-1,lambda=  as.numeric(names(which.min(Ridge$GCV))))
  frac.delta = seq(from = 0, to = 1, 0.01)
  mse.cv.Lasso = cv.lars(trainX,trainY,trace=F,plot.it = F,intercept=F,normalize=F, index = frac.delta, type = c("lasso"))
  lambdaOpt = frac.delta[which.min(mse.cv.Lasso$cv)]  
  Lasso = lars(trainX, trainY, type = c("lasso"), use.Gram=F, intercept=F, normalize=F)
  PLS = plsr(trainY~trainX-1, ncomp = 50, validation = "LOO")
  PLS = plsr(trainY~trainX-1, ncomp = which.min(PLS$validation$PRESS[1,]), validation = "LOO")
  

  # Prediction
  predLinear = testX%*%Ols$coefficients
  predBridge = testX%*%ClassicalBridge
  predBayesianBridgeFix = testX%*%apply(BayesianBridgeFix$beta,2,mean)
  predBayesianBridge = testX%*%apply(BayesianBridge$beta,2,mean)
  predRidge = testX%*%Ridge$coef
  predLasso = predict(Lasso, testX, s=lambdaOpt , mode="fraction")
  predPLS = predict(PLS,newdata=testX)
  
  
  #Somme des Erreurs de prediction
  SSE.Linear=sum((testY-predLinear)^2)
  SSE.ClassicalBridge = sum((testY-predBridge)^2)
  SSE.BayesianBridgeFix = sum((testY-predBayesianBridgeFix)^2)
  SSE.BayesianBridge = sum((testY-predBayesianBridge)^2)
  SSE.Ridge = sum((testY-predRidge)^2)
  SSE.Lasso = sum((testY-predLasso$fit)^2)
  SSE.PLS = sum((testY-predPLS)^2)
  
  
  # On stocke dans la liste
  List.Linear[[i]] =  SSE.Linear
  List.ClassicalBridge[[i]] =  SSE.ClassicalBridge
  List.BayesianBridgeFix[[i]] =   SSE.BayesianBridgeFix
  List.BayesianBridge[[i]] =   SSE.BayesianBridge
  List.Ridge[[i]] = SSE.Ridge
  List.Lasso[[i]] = SSE.Lasso
  List.PLS[[i]] = SSE.PLS
  
  
}



#png("glucose.png",height=700, width=935)
pdf("glucosePred.pdf",width = 11)
boxplot(unlist(List.Linear),unlist(List.Ridge),unlist(List.Lasso),unlist(List.ClassicalBridge),unlist(List.BayesianBridgeFix),unlist(List.BayesianBridge),
        names=c("OLS","Ridge(GCV)","Lasso(GCV)","Bridge(EM+GCV)","Bayes(a=0.5)","Bayes"),col=c("red","paleturquoise4","yellow","blue","green","violet","orange"),ylab="SSE")
title(main="Glucose data",cex.main = 1.5, font.main= 6.7, col.main= "black")
dev.off()
#dev.off()


