#################################################################
#################################################################
########## IMPORTANCE SAMPLING AVEC UN STUDENT MULTIVARIE #######
#################################################################
#################################################################

rm(list=ls())
dev.off()


## loi normal multivariee0
library(mvtnorm)
library(BayesBridge)

ImSamplingWithStudent = function(y,X,nu,alpha,N){  
  #--Definissons quelques constantes reelles et matricielles
  n = nrow(X)
  p = ncol(X)
  xx = t(X)%*%X
  xx.inv = solve(xx)
  xy = t(X)%*%y
  
  #--Calculs preliminaires (MCO,sigma2,SCR,residus,...)
  betaHat = xx.inv%*%xy
  yHat = X%*%betaHat
  uHat = y-yHat
  SCR = sum(uHat^2)
  sigma2Hat = (1/(n-p))*SCR
  #sigma2Hat = (1/n)*SCR
  varBetaHat = sigma2Hat*xx.inv
  
  #--On simule les beta^{(k)} suivant t{xi,mu,SIGMA} 
  betaStudent = rmvt(N, sigma = varBetaHat*1/30, delta = betaHat, df=3, type = "shifted")
  
  #--Fonction qui calcule les poids pour l'estimateur "IS"
  weitStudent = function(beta){
    N = exp(((-1)/(2*sigma2Hat))*sum((y-X%*%beta)^2))*exp(-nu*sum(abs(beta)^{alpha}))
    D = dmvt(beta, sigma = varBetaHat*1/30, delta = betaHat,df=3, log = F, type = "shifted")
    return(N/D)
  }
  
  #--Vecteur des poids w(beta^{k}) empilés horizontalements 
  BigWeigth = apply(betaStudent,1,weitStudent)
  
  #--Calcule des composantes de l'estimateur "IS"
  ISestimate = betaStudent * (BigWeigth)
  ISestimate = apply(ISestimate,2,cumsum)*(1/cumsum(BigWeigth))
  
  #--Calcule du critere ESS(beta^{k})
  ESS = (sum(BigWeigth,na.rm=T))^2/sum(BigWeigth^2,na.rm=T)
  
  #--On retourne les termes de la somme (a calculer) de l'estimateur IS et le vecteur des poids w(beta^{(k)})
  return(list(IS=ISestimate,weit=BigWeigth,ESS=ESS,OLS=betaHat,betaStudent=betaStudent,sigma2=sigma2Hat ))
  
}

#######################
#---Data SImulation---#
#######################

#nombre dindividu
n = 500

# on considère un modele linaire a 4 variable x1, x2, x3, x4
x1 = rnorm(n, mean = 0, sd = 1)
x2 = rnorm(n, mean = 0, sd = 1)
x3 = rnorm(n, mean = 0, sd = 1)
x4 = rnorm(n, mean = 0, sd = 1)

# terme d'erreur
epsilon = rnorm(n)

# les beta's
beta1 = 2
beta2 = 3
beta3 = -1
beta4 = 0.5

# MODELE : y = Xbeta + epsilon
y = beta1*x1 + beta2*x2 + beta3*x3 + beta4*x4  + epsilon
X = cbind(x1,x2,x3,x4)


############################
###---Hyperparametres----###
############################

alpha = 0.5
nu =  0.4674039
#nu = rgamma(1, shape=2, scale=1/2)
#nu = 0.4674039

### Nombre de simulation
N=10000

IS = ImSamplingWithStudent(y,X,nu,alpha,N)
BayesianBridge = bridge.reg(y, X, nsamp=N,sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,method="stable")

betaOLS = IS$OLS
poid = IS$weit
(sum(poid^2))
(sum(poid))^2
betaStudent=IS$betaStudent

pdf(file='DATA_SIMUL_Student_V2.pdf')

boxplot(log(poid),main=expression('Logarithme des Poids : '*log(w(beta^(k)))))

#boxplot(poid,type='l',col='blue',main=expression('Poids : '*w(beta^(k))),xlab="k",ylab="")
#plot(log(poid),type='l',col='blue',main=expression('Poids : '*log(w(beta^(k)))),xlab="k",ylab="")
plot(cumsum(poid)^2/cumsum(poid^2),type='l',xlab="simulations",ylab="taille effective d'echantillon",main="ESS avec Student",lwd=2)


(IS$ESS/N)*100

#estimateur de Bayes :)
Bayes = IS$IS


par(mar = c(5.1, 5.2, 4, 2))

par(mar = c(5.1, 5.2, 4, 2))
burnIn = 1:N


plot(density(BayesianBridge$beta[,1]),col="green",main=colnames(X)[1],ylim=c(0,20))
lines(density(betaStudent[,1],weights=poid/sum(poid)),main=colnames(X)[1],col='blue')
#ab1line(v=mean(ans$beta[,1]), col='red', lty='dashed')
abline(v=Bayes[N,1], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,1]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,2]),col="green",main=colnames(X)[2],ylim=c(0,20))
lines(density(betaStudent[,2],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,2]), col='red', lty='dashed')
abline(v=Bayes[N,2], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,2]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,3]),col="green",main=colnames(X)[3],ylim=c(0,20))
lines(density(betaStudent[,3],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,3]), col='red', lty='dashed')
abline(v=Bayes[N,3], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,3]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,4]),col="green",main=colnames(X)[4],ylim=c(0,20))
lines(density(betaStudent[,4],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,4], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,4]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


dev.off()


##########################
#------------------------#
#-----DIABETE DATA-------#
#------------------------#
##########################

dev.off()
alpha = 0.5
nu =  0.4674039


### Nombre de simulation
N=100000
## DOnnees du Diabete : 10 variables explicatives
data(diabetes, package="lars")
X = diabetes$x
y = diabetes$y
y = (y - mean(y))
p = ncol(X)
n = length(y)
for(j in 1:p)
{
  X[,j] = (X[,j] - mean(X[,j]))
}


IS = ImSamplingWithStudent(y,X,nu,alpha,N)
BayesianBridge = bridge.reg(y, X, nsamp=N,sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,method="stable")

poid = IS$weit

################################
### ON ENREEGISTRE EN PDF(.) ###
################################
dev.off()

pdf(file='Importance_Sampling_With_Student.pdf')


boxplot(log(poid),main=expression('Logarithme des Poids : '*log(w(beta^(k)))))
plot(cumsum(poid)^2/cumsum(poid^2),type='l',xlab="simulations",ylab="taille effective d'echantillon",main="ESS avec Student",lwd=2)

plot(poid,type='l',col='blue',main=expression('Poids : '*w(beta^(k))),xlab="k",ylab="")
plot(log(poid),type='l',col='blue',main=expression('Poids : '*log(w(beta^(k)))),xlab="k",ylab="")
(IS$ESS/N)*100


#estimateur de Bayes :)
Bayes = IS$IS
par(mar = c(5.1, 5.2, 4, 2))
#par(mfrow=c(1,2))
betaOLS = IS$OLS

plot(Bayes[,1],type='l',col='blue',main=colnames(X)[1],ylab=expression(hat(beta[B])^{(m)}),xlab="k",cex.main = 1.3)
abline(b=0,a=betaOLS[1],col='red')
#abline(b=0,a=mean(ans$beta[,1]), col='green', lty='dashed')
#lines(cumsum(ans$beta[,1])/1:length(ans$beta[,1]),col='green')
lines(cumsum(BayesianBridge$beta[,1])/1:N,type='l',col='green')
legend("bottomright",c("IS","OLS","Gibbs"),lwd=1,col=c("blue","red","green"))


plot(Bayes[,2],type='l',col='blue',main=colnames(X)[2],ylab=expression(hat(beta[B])^{(m)}),xlab="k",cex.main = 1.3,ylim=c(-280,-100))
abline(b=0,a=betaOLS[2],col='red')
#abline(b=0,a=mean(ans$beta[,2]), col='green', lty='dashed')
#lines(cumsum(ans$beta[,2])/1:length(ans$beta[,2]),col='green')
lines(cumsum(BayesianBridge$beta[,2])/1:N,type='l',col='green')
legend("bottomright",c("IS","OLS","Gibbs"),lwd=1,col=c("blue","red","green"))


plot(Bayes[,3],type='l',col='blue',main=colnames(X)[3],ylab=expression(hat(beta[B])^{(m)}),xlab="k",cex.main = 1.3,ylim=c(400,600))
abline(b=0,a=betaOLS[3],col='red')
#abline(b=0,a=mean(ans$beta[,3]), col='green', lty='dashed')
#lines(cumsum(ans$beta[,3])/1:length(ans$beta[,3]),col='green')
lines(cumsum(BayesianBridge$beta[,3])/1:N,type='l',col='green')
legend("bottomright",c("IS","OLS","Gibbs"),lwd=1,col=c("blue","red","green"))


plot(Bayes[,4],type='l',col='blue',main=colnames(X)[4],ylab=expression(hat(beta[B])^{(m)}),xlab="k",cex.main = 1.3)
abline(b=0,a=betaOLS[4],col='red')
#abline(b=0,a=mean(ans$beta[,3]), col='green', lty='dashed')
#lines(cumsum(ans$beta[,4])/1:length(ans$beta[,4]),col='green')
lines(cumsum(BayesianBridge$beta[,4])/1:N,type='l',col='green')
legend("topright",c("IS","OLS","Gibbs"),lwd=1,col=c("blue","red","green"))


plot(Bayes[,5],type='l',col='blue',main=colnames(X)[5],ylab=expression(hat(beta[B])^{(m)}),xlab="k",cex.main = 1.3,ylim=c(-900,mean(BayesianBridge$beta[,5])+100))
abline(b=0,a=betaOLS[5],col='red')
#abline(b=0,a=mean(ans$beta[,5]), col='green', lty='dashed')
#lines(cumsum(ans$beta[,5])/1:length(ans$beta[,5]),col='green')
lines(cumsum(BayesianBridge$beta[,5])/1:N,type='l',col='green')
legend("bottomright",c("IS","OLS","Gibbs"),lwd=1,col=c("blue","red","green"))



plot(Bayes[,6],type='l',col='blue',main=colnames(X)[6],ylab=expression(hat(beta[B])^{(m)}),xlab="k",cex.main = 1.3)
abline(b=0,a=betaOLS[6],col='red')
#abline(b=0,a=mean(ans$beta[,6]), col='green', lty='dashed')
#lines(cumsum(ans$beta[,6])/1:length(ans$beta[,6]),col='green')
lines(cumsum(BayesianBridge$beta[,6])/1:N,type='l',col='green')
legend("topright",c("IS","OLS","Gibbs"),lwd=1,col=c("blue","red","green"))



plot(Bayes[,7],type='l',col='blue',main=colnames(X)[7],ylab=expression(hat(beta[B])^{(m)}),xlab="k",cex.main = 1.3)
abline(b=0,a=betaOLS[7],col='red')
#abline(b=0,a=mean(ans$beta[,7]), col='green', lty='dashed')
#lines(cumsum(ans$beta[,7])/1:length(ans$beta[,7]),col='green')
#lines(cumsum(BayesianBridge$beta[,7])/1:N,type='l',col='green')
lines(cumsum(BayesianBridge$beta[,7])/1:N,type='l',col='green')
legend("topright",c("IS","OLS","Gibbs"),lwd=1,col=c("blue","red","green"))



plot(Bayes[,8],type='l',col='blue',main=colnames(X)[8],ylab=expression(hat(beta[B])^{(m)}),xlab="k",cex.main = 1.3,ylim=c(0,200))
abline(b=0,a=betaOLS[8],col='red')
#abline(b=0,a=mean(ans$beta[,8]), col='green', lty='dashed')
#lines(cumsum(ans$beta[,8])/1:length(ans$beta[,8]),col='green')
lines(cumsum(BayesianBridge$beta[,8])/1:N,type='l',col='green')
legend("topright",c("IS","OLS","Gibbs"),lwd=1,col=c("blue","red","green"))



plot(Bayes[,9],type='l',col='blue',main=colnames(X)[9],ylab=expression(hat(beta[B])^{(m)}),xlab="k",cex.main = 1.3,ylim=c(250,800))
abline(b=0,a=betaOLS[9],col='red')
#abline(b=0,a=mean(ans$beta[,8]), col='green', lty='dashed')
#lines(cumsum(ans$beta[,9])/1:length(ans$beta[,9]),col='green')
lines(cumsum(BayesianBridge$beta[,9])/1:N,type='l',col='green')
legend("topright",c("IS","OLS","Gibbs"),lwd=1,col=c("blue","red","green"))



plot(Bayes[,10],type='l',col='blue',main=colnames(X)[10],ylab=expression(hat(beta[B])^{(m)}),xlab="k",cex.main = 1.3)
abline(b=0,a=betaOLS[10],col='red')
#abline(b=0,a=mean(ans$beta[,8]), col='green', lty='dashed')
#lines(cumsum(ans$beta[,10])/1:length(ans$beta[,10]),col='green')
lines(cumsum(BayesianBridge$beta[,10])/1:N,type='l',col='green')
legend("topright",c("IS","OLS","Gibbs"),lwd=1,col=c("blue","red","green"))


##############################################################
### On recommence mais avec les densités non parametrique!!--#
##############################################################

dev.off()

alpha = 0.5
nu =  0.4674039
tau = nu^{-1/alpha}
#nu = rgamma(1, shape=2, scale=1/2)
#nu = 0.4674039

### Nombre de simulation
N=1000000
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



IS = ImSamplingWithStudent(y,X,nu,alpha,N)
BayesianBridge = bridge.reg(y,X,alpha=alpha ,nsamp=N,sig2.true=IS$sigma2,tau.true=tau,method="stable")
poid = IS$weit
(IS$ESS/N)*100
  

#ESS = 0.02%




#estimateur de Bayes :)
Bayes = IS$IS
par(mar = c(5.1, 5.2, 4, 2))
#par(mfrow=c(1,2))
betaOLS = IS$OLS
#les betas simulées 
betaStudent = IS$betaStudent
dim(betaStudent)


pdf(file='Importance_Sampling_With_Student_Diabete_Data.pdf')

boxplot(log(poid),main=expression('Logarithme des Poids : '*log(w(beta^(k)))))
plot(cumsum(poid)^2/cumsum(poid^2),type='l',xlab="simulations",ylab="taille effective d'echantillon",lwd=2,main="ESS fontion d'importance Student")

#mylabs = c("Age", "Female", "Body Mass Index", "Blood Pressure", "TC", "LDL", "HDL", "TCH", "LTG", "Glucose")
colnames(X)
par(mar = c(5.1, 5.2, 4, 2))

burnIn = 1:N

plot(density(BayesianBridge$beta[,1]),col="green",main=colnames(X)[1],ylim=c(0,0.25))
lines(density(betaStudent[,1],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,1]), col='red', lty='dashed')
abline(v=Bayes[N,1], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,1]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


#plot(density(ans$beta[,2]),col="green",main=colnames(X)[2],xlim=c(-250,-100),ylim=c(0,0.025))
plot(density(BayesianBridge$beta[,2]),col="green",main=colnames(X)[2],ylim=c(0,0.145))
lines(density(betaStudent[,2],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,2]), col='red', lty='dashed')
abline(v=Bayes[N,2], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,2]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,3]),col="green",ylim=c(0,0.14),main=colnames(X)[3])
lines(density(betaStudent[,3],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,3]), col='red', lty='dashed')
abline(v=Bayes[N,3], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,3]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,4]),col="green",main=colnames(X)[4],ylim=c(0,0.14))
lines(density(betaStudent[,4],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,4], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,4]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))



plot(density(BayesianBridge$beta[,5]),col="green",main=colnames(X)[5],ylim=c(0,0.08))
lines(density(betaStudent[,5],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,5], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,5]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,6]),col="green",main=colnames(X)[6],ylim=c(0,0.18))
lines(density(betaStudent[,6],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,6], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,6]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))



plot(density(BayesianBridge$beta[,7]),col="green",main=colnames(X)[7],ylim=c(0,0.112))
lines(density(betaStudent[,7],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,7], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,7]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,8]),col="green",main=colnames(X)[8],ylim=c(0,0.112))
lines(density(betaStudent[,8],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,8], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,8]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,9]),col="green",main=colnames(X)[9],ylim=c(0,0.115))
lines(density(betaStudent[,9],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,9], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,9]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))



plot(density(BayesianBridge$beta[,10]),col="green",ylim=c(0,0.155),main=colnames(X)[10])
lines(density(betaStudent[,10],weights=poid/sum(poid)),main=colnames(X)[10],col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,10], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,10]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


dev.off()


#############################
#---------------------------#
#-----Ozone Data DATA-------#
#---------------------------#
#############################


alpha = 0.5
nu =  0.4674039
#nu = rgamma(1, shape=2, scale=1/2)
#nu = 0.4674039

### Nombre de simulation
N=30000

data("Ozone", package="mlbench")
summary(Ozone)
class(Ozone)
X = Ozone[,5:13]
colnames(Ozone)
dim(X)
## Plot des corrélations des 
#corrplot(cor(X,use = "complete.obs"),order="hclust",addCoef.col="grey")
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


IS = ImSamplingWithStudent(y,X,nu,alpha,N)
BayesianBridge = bridge.reg(y, X, nsamp=N,sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,method="stable")
poid = IS$weit

#estimateur de Bayes :)
Bayes = IS$IS
par(mar = c(5.1, 5.2, 4, 2))
#par(mfrow=c(1,2))
betaOLS = IS$OLS
#les betas simulées 
betaStudent = IS$betaStudent
dim(betaStudent)
(IS$ESS/N)*100


pdf(file='Importance_Sampling_With_Ozone_Data_Student.pdf')


boxplot(log(poid),main=expression('Logarithme des Poids : '*log(w(beta^(k)))))
plot(cumsum(poid)^2/cumsum(poid^2),type='l',xlab="simulations",ylab="taille effective d'echantillon",lwd=2,main="ESS fontion d'importance Normal")

#mylabs = c("Age", "Female", "Body Mass Index", "Blood Pressure", "TC", "LDL", "HDL", "TCH", "LTG", "Glucose")
colnames(X)
par(mar = c(5.1, 5.2, 4, 2))

burnIn = 1:N

plot(density(BayesianBridge$beta[,1]),col="green",main=colnames(X)[1])
lines(density(betaStudent[,1],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,1]), col='red', lty='dashed')
abline(v=Bayes[N,1], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,1]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,2]),col="green",main=colnames(X)[2])
lines(density(betaStudent[,2],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,2]), col='red', lty='dashed')
abline(v=Bayes[N,2], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,2]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,3]),col="green",main=colnames(X)[3])
lines(density(betaStudent[,3],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,3]), col='red', lty='dashed')
abline(v=Bayes[N,3], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,3]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,4]),col="green",main=colnames(X)[4])
lines(density(betaStudent[,4],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,4], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,4]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,5]),col="green",main=colnames(X)[5])
lines(density(betaStudent[,5],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,5], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,5]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topleft",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,6]),col="green",main=colnames(X)[6])
lines(density(betaStudent[,6],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,6], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,6]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,7]),col="green",main=colnames(X)[7])
lines(density(betaStudent[,7],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,7], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,7]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,8]),col="green",main=colnames(X)[8])
lines(density(betaStudent[,8],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,8], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,8]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))


plot(density(BayesianBridge$beta[,9]),col="green",main=colnames(X)[9])
lines(density(betaStudent[,9],weights=poid/sum(poid)),col='blue')
#abline(v=mean(ans$beta[,4]), col='red', lty='dashed')
abline(v=Bayes[N,9], col='yellow',lwd=3)
abline(v=mean(BayesianBridge$beta[,9]), col='red', lty='dashed')
#abline(v=betaOLS[1], col='red')
legend("topright",c("IS","Gibbs","Gibbs estimator","IS estimator"),lwd=1,col=c("blue","green","red","yellow"))



dev.off()

