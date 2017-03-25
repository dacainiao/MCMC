#1.a
library(mvtnorm)
logic <- function(nn=100, X=cbind(rep(1,nn), rbinom(nn,1,0.5),rbinom(nn,1,0.5)), 
beta=NULL, priormean=rep(0,3), priorvar=diag(1000,3)) {
if (is.null(beta)) beta <- t(rmvnorm(1, priormean, priorvar)) 
yy <- rbinom(nn, 1, exp(X%*%beta)/(1+exp(X%*%beta)))
return(list(yy=yy, X=X))
}


#1.b
loglike <- function(beta, X, yy) 
sum(dbinom(yy, 1, exp(X%*%beta)/(1+exp(X%*%beta)), log=TRUE))
prior <- function(beta, priormean, priorvar)
dmvnorm(beta, priormean, priorvar, log=TRUE)
prop <-function(index,beta.old,width,X,yy,priormean=rep(0,dim(X)[2]),priorvar=diag(1000,dim(X)[2])) {
beta.prop <- beta.old
beta.prop[index] <- rnorm(1,beta.old[index],width)
ratio <-loglike(beta.prop, X, yy)+prior(beta.prop,priormean,priorvar)-loglike(beta.old,X,yy)-prior(beta.old, priormean,priorvar)
accepted <- rbinom(1,1,min(1,exp(ratio)))
return(beta.prop*accepted + beta.old*(1-accepted))
}


#1.c
test<- logic(nn=50,beta=c(-1,0,1))
testwidth <- function(index, width) {
beta.old <- rnorm(3,0,1)
num <- 0
for (kk in 1:1000) {
beta.prop <- prop(index, beta.old, width,test$X, test$yy)
if (beta.prop[index]!=beta.old[index]) num=num+1
beta.old <- beta.prop
}
acptrate <- num/1000
return(acptrate)
}
testwidth(1,1)
testwidth(2,1)
testwidth(3,1)


#1.d
gibbs <- function(yy, X, iterations=1000, thin=1,
burnin=100,priormean=rep(0,dim(X)[2]),priorvar=diag(1000, dim(X)[2]),
beta.start=rnorm(dim(X)[2], 0, 1)) {
width <- 1
beta.out <- array(NA, c(dim(X)[2], iterations))
beta <- beta.start
for (ii in 1:(burnin+thin*iterations)) {
for (jj in 1:length(beta)) {
beta.prop <- prop(index=jj, beta, width, X, yy)
beta <- beta.prop
}
iter <- (ii-burnin)/thin
if (iter>0 & trunc(iter)==iter) beta.out[,iter] <- beta
}
return (beta.out)
}

quantiles <- array(NA, c(3, 30))
for (kk in 1:30) {
trues <- rmvnorm(1, cbind(rep(0,3)), diag(1, 3))
fake <- logic(beta=t(trues))
runner <- gibbs(fake$yy, fake$X, priormean=rep(0,3),priorvar=diag(1, 3))
quantiles[,kk] <- c(mean(runner[1,]>trues[1]), 
                    mean(runner[2,]>trues[2]),
                    mean(runner[3,]>trues[3]))
}
par(mfrow=c(3,1))
hist(quantiles[1,])
hist(quantiles[2,])
hist(quantiles[3,])

gibbs2 <- function(yy, X, iterations=1000, thin=1, width=1,
burnin=100,priormean=rep(0,dim(X)[2]),priorvar=diag(1000, dim(X)[2]),
beta.start=rnorm(dim(X)[2], 0, 1)) {
beta.out <- array(NA, c(dim(X)[2], iterations))
beta <- beta.start
for (ii in 1:(burnin+thin*iterations)) {
for (jj in 1:length(beta)) {
beta.prop <- prop(index=jj, beta, width, X, yy)
beta <- beta.prop
}
iter <- (ii-burnin)/thin
if (iter>0 & trunc(iter)==iter) beta.out[,iter] <- beta
}
rate1<-mean(beta.out[1,][-1]!=beta.out[1,][-1000])
rate2<-mean(beta.out[2,][-1]!=beta.out[2,][-1000])
rate3<-mean(beta.out[3,][-1]!=beta.out[3,][-1000])
return (c(rate1,rate2,rate3))
}
gibbs2(fake$yy, fake$X, width=1, priormean=rep(0,3),priorvar=diag(1, 3))
gibbs2(fake$yy, fake$X, width=0.5, priormean=rep(0,3),priorvar=diag(1, 3))
gibbs2(fake$yy, fake$X, width=3, priormean=rep(0,3),priorvar=diag(1, 3))


#1.e
data1 <-read.table("lupus.txt")
outcome <-data1[,1]
pred <-cbind(rep(1,55),data1[,2],data1[,3])
chain1<-gibbs(outcome,pred,priormean=rep(0,3),priorvar=diag(1,3)¡ê?beta.start=rnorm(3,0,1))
report1<-array(NA,c(2,3))
for (pp in 1:3){
report1[,pp]<-c(mean(chain1[pp,]),sd(chain1[pp,]))
}
rownames(report1)<-c("Mean","SD")
colnames(report1)<-c("intercept","IGD","IA")

chain2<-gibbs(outcome,pred,priormean=rep(0,3),priorvar=diag(1,3)¡ê?beta.start=c(0,10,-10))
report2<-array(NA,c(2,3))
for (pp in 1:3){
report2[,pp]<-c(mean(chain2[pp,]),sd(chain2[pp,]))
}
rownames(report2)<-c("Mean","SD")
colnames(report2)<-c("intercept","IGD","IA")


chain3<-gibbs(outcome,pred,priormean=rep(0,3),priorvar=diag(1,3)¡ê?beta.start=c(0,50,-50))
report3<-array(NA,c(2,3))
for (pp in 1:3){
report3[,pp]<-c(mean(chain3[pp,]),sd(chain3[pp,]))
}
rownames(report3)<-c("Mean","SD")
colnames(report3)<-c("intercept","IGD","IA")

chain4<-gibbs(outcome,pred,priormean=rep(0,3),priorvar=diag(1,3)¡ê?beta.start=c(0,100,-100))
report4<-array(NA,c(2,3))
for (pp in 1:3){
report4[,pp]<-c(mean(chain4[pp,]),sd(chain4[pp,]))
}
rownames(report4)<-c("Mean","SD")
colnames(report4)<-c("intercept","IGD","IA")

     intercept      IGD       IA
Mean -5.227082 12.27425 6.960843
SD    2.071431  3.79081 2.742859
> report2
     intercept       IGD        IA
Mean -8.047809 17.374547 10.540647
SD    2.549493  4.961341  3.116832
> report3
     intercept      IGD       IA
Mean -7.191967 17.93732 9.623171
SD    4.010897  8.98511 5.542036
> report4
     intercept IGD  IA
Mean       NaN NaN NaN
SD          NA  NA  NA






#2.a
data2 <- read.csv("efron-morris-batting.csv")
yy<-data2[,3]
nn<-data2[,1]
gibbs3 <- function (yy, nn, cc.prior=c(60,2), phi.prior=c(2,20),pp.start=NULL, cc.start=NULL, phi.start=NULL,
  iterations=100, thin=10, burnin=100) {
  count <- length(yy)
  cc.tune <- 1
  phi.tune <- 10
  pp.out <- array(NA, c(count, iterations))
  cc.out <- phi.out <- rep(NA, iterations)
  if (is.null(cc.start)) cc <- rgamma(1, cc.prior[1], cc.prior[2]) else cc <- cc.start
  if (is.null(phi.start)) phi <- rbeta(1, phi.prior[1], phi.prior[2]) else phi <- phi.start
  if (is.null(pp.start)) pp <- rbeta(count, cc*phi, cc*(1-phi)) else pp <- pp.start
  for (iter in 1:(burnin+thin*iterations)) {
# sample p    
    pp <- rbeta(count, cc*phi + yy, cc*(1-phi) + nn - yy)
# sample c   
    cc.prop <- rgamma(1, cc*cc.tune, cc.tune)
    loglik.cc.orig <- sum(dbeta(pp, cc*phi, cc*(1-phi), log=TRUE))
    loglik.cc.prop <- sum(dbeta(pp, cc.prop*phi, cc.prop*(1-phi), log=TRUE))
    log.acc.rat <- loglik.cc.prop - loglik.cc.orig +
    dgamma(cc.prop,cc.prior[1],cc.prior[2], log=TRUE)-dgamma(cc,cc.prior[1],cc.prior[2], log=TRUE)+      
    dgamma(cc, cc.prop*cc.tune, cc.tune, log=TRUE) - dgamma(cc.prop, cc*cc.tune, cc.tune, log=TRUE) 
    trial <- 1*(log(runif(1)) < log.acc.rat)
    cc <- trial*cc.prop + (1-trial)*cc
# sample phi
    phi.prop <- rbeta(1, phi.tune*phi, phi.tune*(1-phi))
    loglik.phi.orig <- sum(dbeta(pp, cc*phi, cc*(1-phi), log=TRUE))
    loglik.phi.prop <- sum(dbeta(pp, cc*phi.prop, cc*(1-phi.prop), log=TRUE))
    log.acc.rat <- loglik.phi.prop - loglik.phi.orig +
    dbeta(phi.prop, phi.prior[1], phi.prior[2], log=TRUE)-dbeta(phi, phi.prior[1], phi.prior[2], log=TRUE)+
    dbeta(phi, phi.tune*phi.prop, phi.tune*(1-phi.prop), log=TRUE) - dbeta(phi.prop, phi.tune*phi, phi.tune*(1-phi), log=TRUE)
    trial <- 1*(log(runif(1)) < log.acc.rat)
    phi <- trial*phi.prop + (1-trial)*phi
    index.point <- (iter-burnin)/thin
    if (index.point > 0 & index.point == trunc(index.point)) {
    pp.out[,index.point] <- pp
    cc.out[index.point] <- cc
    phi.out[index.point] <- phi
    }
  } 
  aa=cc.out*phi.out
  bb=cc.out*(1-phi.out)
  out <- list(pp=pp.out,aa=aa,bb=bb)
  return(out)
}
gibbs3(yy,nn)


#2.b
p=data2[,2]/1000
estab<-function(p) {
n=18
p1=mean(p)
p2=var(p)*(n-1)/n
a<-p1*(p1*(1-p1)/p2-1)
b<-(1-p1)*a/p1
return(c(a,b))
}
jsestmu <- mean(yy)+(1-(n-2)*var(yy)/sum((yy-mean(yy))^2))*(yy-mean(yy))
jsestp <-jsestmu/45

#2.c
tt<-array(NA,c(18,1000))
quant<-array(NA,c(18,5))
a<-11
b<-30.5
for (i in 1:18) {
tt[i,]<-rbeta(1000,a+yy[i],b+45+yy[i])
quant[i,]<-quantile(tt[i,],c(0.025,0.25,0.5,0.75,0.975),name=FALSE)
}
colnames(quant)<-c(0.025,0.25,0.5,0.75,0.975)

#2.d
new<-data2[,5]/1000
# (0.25,0.75)
frac1=0
for (i in 1:18){
if (new[i]> quant[i,2] & new[i]<quant[i,4])  frac1=frac1+1
}
frac1=frac1/18
# (0.025,0.975)
frac2=0
for (i in 1:18){
if (new[i]> quant[i,1] & new[i]<quant[i,5])  frac2=frac2+1
}
frac2=frac2/18

