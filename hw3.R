1.
start=runif(1)
target <- function(kk) dbeta(kk,6,4)
prop <- function(c,xold) rbeta(1,c*xold,c*(1-xold))
propden <- function(x,c,xold) dbeta(x,c*xold,c*(1-xold))
mh <- function(nn=1000, c=1) {
chain <- rep(NA,nn)
beta <- start
for (kk in 1:nn) {
beta.try <- 1; reps <- 0; while (beta.try == 0 | beta.try == 1) {
beta.try <- prop(c, beta)
reps <- reps+1
if (reps > 1000) stop ("Sampler is jamming.")
}
acpt <- target(beta.try)/target(beta)*propden(beta,c,beta.try)/propden(beta.try,c,beta)
if (runif(1)<acpt) beta <- beta.try
chain[kk] <- beta
}
return(chain)
}
mhdraw <- mh(c=1)
par(mfrow=c(1,3))
plot(mhdraw)
acf(mhdraw)
hist(mhdraw)
mhdraw <- mh(c=0.1)
par(mfrow=c(1,3))
plot(mhdraw)
acf(mhdraw)
hist(mhdraw)
mhdraw <- mh(c=2.5)
par(mfrow=c(1,3))
plot(mhdraw)
acf(mhdraw)
hist(mhdraw)
mhdraw <- mh(c=10)
par(mfrow=c(1,3))
plot(mhdraw)
acf(mhdraw)
hist(mhdraw)


2.b
chain <- function (nn=1000, start=c(0,0,0), thin=2, burnin=10) {
  yy.mean <- c(1,2,3)
  sig1 <- array(4*c(1,0.3,0.3,1), c(2,2))
  var <- 4 - 4*rbind(c(0.3,0.3)) %*% solve(sig1) %*% cbind(c(0.3, 0.3))*4
  yy.out <- array(NA, c(3, nn))
  yy <- start
  for (kk in 1:(burnin+thin*nn)) {
    for (ll in 1:3)
      {
      yy[ll] <- rnorm(1, yy.mean[ll] + rbind(c(1.2,1.2)) %*% solve(sig1)%*%cbind(yy[-ll]-yy.mean[-ll]),
       sqrt(var))
      }
  
  if (kk>burnin & (kk-burnin)/thin == floor((kk-burnin)/thin)) yy.out[,(kk-burnin)/thin]<-yy
  }

  rownames(yy.out) <- c("yy1","yy2","yy3")
  return(yy.out)
}
chain(1000,c(0,0,0),2,10)->t
par(mfrow=c(1,3))
for (i in 1:3 ) acf(t[i,],100)

2.c
library(coda)
tenchain <- list(NA)
for (kk in 1:10) tenchain[[kk]] <- mcmc(t(chain(start=rnorm(3,0,10))))
tenchain <- as.mcmc.list(tenchain)
gelman.diag(tenchain)


2.d
library(mvtnorm)
mean <- c(1,2,3)
var <- array(4*c(1,0.3,0.3,0.3,1,0.3,0.3,0.3,1), c(3,3))
truth.draw <- rmvnorm(1000, mean, var)
par(mfrow=c(1,3))
for (p in 1:3) {
qqplot(tenchain[[6]][,p], truth.draw[,p],xlab="Gibbs Sampler", ylab="Multivariate Draw")
abline(a=0,b=1,col=2,lwd=3)}


3.a
gammadraw <- rgamma(100, 5, 5)
mean(gammadraw)-> meangamma
var(gammadraw)-> vargamma

3.b
joint <- function(xx,yy) {
dgamma(xx,1,0.01)*dgamma(yy,1,0.01)*prod(dgamma(gammadraw, xx, rate=yy))}
xx <- yy <- seq(0.025,14.975,by=0.05)
nn <- length(xx)
densities <- array(NA, c(nn,nn))
for (kk in 1:nn) for (jj in 1:nn) {
densities[kk,jj] <- joint(xx[kk],yy[jj])
}
nn2 <- 10000
multidraw <- sample(1:length(densities), nn2, replace=TRUE, prob=as.vector(densities))
draws <- cbind(xx[(multidraw-1)%%nn+1], yy[trunc((multidraw-1)/nn)+1])
draws <- draws+runif(2*nn2, -0.025, 0.025)
plot(draws)
mean(draws[,1]<5)
mean(draws[,2]<5)

3.c
density.sum <- sum(densities)*0.05*0.05

3.e
max(densities/density.sum)
prop1=runif(100)
prop2=runif(100)
M<-1.2
target.3e<- function(xx,yy) {
dgamma(xx,1,0.01)*dgamma(yy,1,0.01)*prod(dgamma(gammadraw, xx, rate=yy))/density.sum
}
acceptance.rates=target.3e(prop1,prop2)/M
mean(target.3e(prop1,prop2)/M)


3.f
alpha <- 5
nn <- length(gammadraw)
beta.all <- rep(NA, 1000)
for (kk in 1:1000) {
beta <- rgamma(1, nn*alpha+1, 0.01+sum(gammadraw))
beta.all[kk] <- beta
}
mean(beta.all)
var(beta.all)
mean(beta.all<5))

3.g
beta <-5
alpha<-1
width <- 1
thin <- 5
burnin <- 100
alpha.prop <- rnorm(1, alpha, width)
length.out <- 1000*thin+burnin
alpha.out <- rep(NA, length.out)
for (kk in 1:length.out) {
alpha.prop <- rnorm(1, alpha, width)
if (alpha.prop>0) {
old <- sum(dgamma(gammadraw, alpha, beta, log=TRUE))+
dgamma(alpha, 1, 0.01, log=TRUE)
prop <- sum(dgamma(gammadraw, alpha.prop, beta, log=TRUE))+
dgamma(alpha, 1, 0.01, log=TRUE)
accept <- rbinom (1, 1, min(1, exp(prop-old)))
alpha <- accept*alpha.prop + (1-accept)*alpha
}
alpha.out[kk] <- alpha
}
alpha.out <- alpha.out[seq(burnin+thin, length.out, by=thin)]
c(mean(alpha.out), var(alpha.out), mean(alpha.out< 5))


3.h
beta <- 5
width <- 1
alpha <- 1
thin <- 5
burnin <- 100
nn=100
length.out <- 1000*thin+burnin
alpha.out <- beta.out <- rep(NA, length.out)
for (kk in 1:length.out) {
beta <- rgamma(1, nn*alpha+1, 0.01+sum(gammadraw))
alpha.prop <- rnorm(1, alpha, prop.width)
if (alpha.prop>0) {
old2 <- sum(dgamma(gammadraw, alpha, beta, log=TRUE))+
dgamma(alpha, 1, 0.01, log=TRUE)
prop2 <- sum(dgamma(gammadraw, alpha.prop, beta, log=TRUE))+
dgamma(alpha, 1, 0.01, log=TRUE)
accept <- rbinom (1, 1, min(1, exp(prop2-old2)))
alpha <- accept*alpha.prop + (1-accept)*alpha
}
beta.out[kk] <- beta
alpha.out[kk] <- alpha
}
beta.out <- beta.out[seq(burnin+thin, length.out, by=thin)]
alpha.out <- alpha.out[seq(burnin+thin, length.out, by=thin)]
c(mean(alpha.out), var(alpha.out), mean(alpha.out < 5))
c(mean(beta.out), var(beta.out), mean(beta.out < 5))

par(mfrow=c(1,2))
qqplot(draws[,1],alpha.out)
abline(a=0,b=1,col=2,lwd=3)
qqplot(draws[,2],beta.out)
abline(a=0,b=1,col=2,lwd=3)
cov(alpha.out,beta.out)
cov(draws[,1],draws[,2])