2.b
u <- runif(10000)
y = u^(1/5)
par(mfrow=c(2,1))
hist(y)
beta.y <-rbeta(10000,5,1)
hist(beta.y)

3.b
gams <- rgamma(10000, shape=5/2, rate=5/2)
draws <-rnorm(10000, 0, sqrt(1/gams))
par(mfrow=c(2,1))
hist(draws)
t.draw <- rt(10000, 5)
hist(t.draw)

3.d
gams <- rgamma(10000, 6, 2)
draws <- rpois(10000, gams)
par(mfrow=c(2,1))
hist(draws)
nb.draw <- rnbinom(10000,6,2/3)
hist(nb.draw)

5.4
cbind(seq(0.52,0.54,length=100), dbeta(seq(0.52,0.54,length=100),321,281))

6.a
x<- seq(0,1,0.00001)
plot(x,(1/3)*dbeta(x,6,2)+(2/3)*dbeta(x,2,6))

6.b
pdf <- function(x) 1/3*dbeta(x,6,2)+2/3*dbeta(x,2,6)
a <- seq(0,1,length=100001)
a.pdf <- pdf(a)
argmax <- a[which(a.pdf==max(a.pdf))]
print(c(argmax,max(a.pdf)))
n <- 100000
M <- 1.877
proposal <- runif(n)
trials <- 1*(runif(n)<(pdf(proposal)/M/dunif(proposal)))
draws <- proposal[which(trials==1)]
hist(draws)

6.c
efficiency <- length(draws)/n
print(c(efficiency,1/M))

6.d
k <- length(draws)
beta1 <- rbinom(1,k,1/3)
draws2 <- c(rbeta(beta1,6,2), rbeta(k-beta1,2,6))
x11()
qqplot(draws,draws2)


