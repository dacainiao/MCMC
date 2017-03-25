1.c
PP<-array(c(0.75,0.5,0.25,0.5),c(2,2))
eigen(t(PP))

1.d
twostate.markov <- function (nn, initial.draw=1, p1 = 0.75, p2 = 0.5) {

  current.draw <- initial.draw
  draw.output <- rep(NA, nn)
  for (kk in 1:nn) {
    current.draw <- rbinom(1,1,p1)*(current.draw==1)+rbinom(1,1,p2)*(current.draw==0)
    draw.output[kk] <- current.draw
  }

  return (draw.output)
}
jiminy <- twostate.markov (110)
jim <- jiminy[11:110]
sum(jim)

1.e
r<-function(n){  
i<-1
t<-PP
while(i<n){
t<-t%*%PP
i<-i+1}
t}
r(50)
r(51)



2.c
PP<- array(c(0,1/2,1/3,1/2,1/3,0,1/3,0,1/3,1/2,0,1/2,1/3,0,1/3,0),c(4,4))
eigen(t(PP))

2.d
fourstate.markov <- function (nn, initial.draw=1) {

  current.draw <- initial.draw
  draw.output <- rep(NA, nn)
  PP <- array(c(0,1/2,1/3,1/2,1/3,0,1/3,0,1/3,1/2,0,1/2,1/3,0,1/3,0),c(4,4))
  for (kk in 1:nn) {
    current.draw <- sample(1:4,1,prob=PP[current.draw,])
    draw.output[kk] <- current.draw
  }
  return (draw.output)
}

cricket <- fourstate.markov(550)
cri<-cricket[51:550]
length(cri[cri==1])
length(cri[cri==2])
length(cri[cri==3])
length(cri[cri==4])

2.e
r<-function(n){  
i<-1
t<-PP
while(i<n){
t<-t%*%PP
i<-i+1}
t}
r(50)
r(51)

3.b
rbeta(1,16.5,14.5)
rbeta(1,54.5,16.5)
PP<-array(c(0.491,0.190,0.509,0.810),c(2,2))
lamda<-eigen(t(PP))$vectors[,1]
lamda<-lamda/sum(lamda)

3.c
obs <- array(c(14,16,16,54),c(2,2))
simu <- function(nn) {
outcome <- replicate(nn,{
p0 <- rbeta(1,0.5+obs[1,2],0.5+obs[1,1])
p1 <- rbeta(1,0.5+obs[2,2],0.5+obs[2,1])
tt <- array(c(1-p0,p0, 1-p1,p1), c(2,2))
ei <- eigen(tt)$vectors[,1]
ei <- ei/sum(ei)
return(ei[2])
})
return(outcome)
}
simu(100)

3.d
hist(simu(100))
abline(v=2/3,col=2)

4.a
obs2 <- table(walker[-1],walker[-length(walker)])

4.b
rdirichlet(1,obs2[1,]+0.5)
rdirichlet(1,obs2[2,]+0.5)
rdirichlet(1,obs2[3,]+0.5)
rdirichlet(1,obs2[4,]+0.5)
eigen(t(PP2))

4.c
simu <- function(nn) {
outcome <- replicate(nn,{
tt <- rbind(rdirichlet(1,obs2[1,]+0.5),
rdirichlet(1,obs2[2,]+0.5),
rdirichlet(1,obs2[3,]+0.5),
rdirichlet(1,obs2[4,]+0.5))
eigie <- eigen(t(tt))$vectors[,1]
eigie <- eigie/sum(eigie)
return(Re(eigie))
})
return(outcome)
}

4.d
z<-simu(100)
state1=z[1,]
hist(state1)
abline(v=0.3,col=2)

5.a
target<- function(t) dt(t,6)
metro<- function(nn,starter=1,propvar) {
  chain <- rep(NA,nn)
  curr.val <- starter
  for (kk in 1:nn) {
    prop <- curr.val+rnorm(1,0,sqrt(propvar))
    accepted <- rbinom(1,1,min(1,target(prop)/target(curr.val)))
    curr.val <- accepted*prop + (1-accepted)*curr.val
    chain[kk] <- curr.val
  }
  return(chain)
}
chain1 <- metro(10000,1,1)
meanrate<- mean(chain1[-1]!=chain1[-10000])
meanabs<-mean(abs(chain1[-1]-chain1[-10000]))
acf(chain1)
plot(chain1)
hist(chain1)



5.b
target<- function(t) dt(t,6)
metro<- function(nn,starter=1,propvar) {
  chain <- rep(NA,nn)
  curr.val <- starter
  for (kk in 1:nn) {
    prop <- curr.val+rnorm(1,0,propvar)
    accepted <- rbinom(1,1,min(1,target(prop)/target(curr.val)))
    curr.val <- accepted*prop + (1-accepted)*curr.val
    chain[kk] <- curr.val
  }
 
  meanac<- mean(chain[-1]!=chain[-10000])
  meanlong<- 1/meanac
  k=ts(chain)
  ll=cor(k[-1],k[-10000])
  return(c(meanac,meanlong,ll))
}

sapply(seq(0.1, 6, by=0.1),FUN=function(propvar) metro(10000,1,propvar))->c
t(c)


6.a
dpareto.simple <- function(xx) {
outcome <- xx
outcome[xx>=1] <- (2/xx[xx>=1]^3 + 3/xx[xx>=1]^4)/2
outcome[xx<1] <- 0
return(outcome)
}
metro.2 <- function(nn=10000,starter=1,propvar) {
chain <- rep(NA,nn)
curr.val <- starter
for (kk in 1:nn) {
prop <- curr.val+rnorm(1,0,sqrt(propvar))
accepted <- rbinom(1,1,min(1,dpareto.simple(prop)/dpareto.simple(curr.val)))
curr.val <- accepted*prop + (1-accepted)*curr.val
chain[kk] <- curr.val
}
return(chain)
}
chain2 <- metro.2(10000,1,1)
meanrate<- mean(chain2[-1]!=chain2[-10000])
meanabs<-mean(abs(chain2[-1]-chain2[-10000]))
acf(chain2)
plot(chain2)
hist(chain2)

6.b
metro2 <- function(nn=10000,starter=1,propvar) {
  chain <- rep(NA,nn)
  curr.val <- starter
  for (kk in 1:nn) {
  prop <- curr.val+rnorm(1,0,sqrt(propvar))
  accepted <- rbinom(1,1,min(1,dpareto.simple(prop)/dpareto.simple(curr.val)))
  curr.val <- accepted*prop + (1-accepted)*curr.val
  chain[kk] <- curr.val
  }
  meanac<- mean(chain[-1]!=chain[-10000])
  meanlong<- 1/meanac
  k=ts(chain)
  ll=cor(k[-1],k[-10000])
  return(c(meanac,meanlong,ll))
}

sapply(c(0.01,0.09,0.25,4,16), FUN=function(propvar) metro2(10000,1,propvar))->c2
t(c2)


 




