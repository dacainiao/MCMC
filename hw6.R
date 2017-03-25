1.a
prob<-rep(NA,100)
for (i in 1:100) {
rt(1000,3)->a
prob[i]<-mean(1*(a>2))
}
prob[1]
var(prob)

1.b
var.es<- function(mu) {
prob<-rep(NA,100)
for (i in 1:100) {
rt(1000,3)+mu->a
prob[i]<-mean((1*(a>2))*dt(a,3)/dt(a-mu,3))
}
out<-var(prob)
return(c(prob[1],out))
}
var.es(0.5)
var.es(1)
var.es(1.5)
var.es(2)
var.es(2.5)
var.es(3)

2.a
ea<-rep(NA,100)
for (i in 1:100) {
rt(1000,3)->candi
draw<-candi[candi>2]
ea[i]<-mean(draw)
}
ea[1]
var(ea)

2.b
var.es2<- function(lamda) {
ea<-rep(NA,100)
for (i in 1:100) {
rexp(1000,lamda)+2->w
ea[i]<-mean(w/pt(-2,3)*dt(w,3)/dexp(w-2,lamda))
}
out<-var(ea)
return(c(ea[1],out))
}
var.es2(0.5)
var.es2(1)
var.es2(1.5)
var.es2(2)
var.es2(2.5)
var.es2(3)

2.c
var.es3<- function(lamda) {
ea<-rep(NA,100)
for (i in 1:100) {
rexp(1000,lamda)+2->w
ea[i]<-sum(w*dt(w,3)/dexp(w-2,lamda))/sum(dt(w,3)/dexp(w-2,lamda))
}
out=var(ea)
return(c(ea[1],out))
}
var.es3(0.5)
var.es3(1)
var.es3(1.5)
var.es3(2)
var.es3(2.5)
var.es3(3)


3.a
VAR=25*(XX%*%t(XX))+diag(1,100)
ml1<-dmvnorm(YY1,rep(0,100),VAR)
ml2<-dmvnorm(YY2,rep(0,100),VAR)

3.b
ev1<-rep(NA,100)
ev2<-rep(NA,100)
for (k in 1:100){
beta<-rmvnorm(1000,rep(0,3),diag(25,3))
like1<-rep(NA,1000)
like2<-rep(NA,1000)
for (i in 1:1000){
like1[i]<-dmvnorm(YY1,XX%*%beta[i,],diag(1,100))
like2[i]<-dmvnorm(YY2,XX%*%beta[i,],diag(1,100))
}
ev1[k]<-mean(like1)
ev2[k]<-mean(like2)
}
ev1[1]
ev2[1]
var(ev1)
var(ev2)

3.c
ev1<-rep(NA,100)
ev2<-rep(NA,100)
betamean<-rep(0,3)
betavar<-diag(25,3)
A<-t(XX)%*%XX+solve(betavar)
B1<-t(XX)%*%YY1+solve(betavar)%*%betamean
B2<-t(XX)%*%YY2+solve(betavar)%*%betamean
for (k in 1:100){
beta1<-rmvnorm(1000,solve(A)%*%B1,solve(A))
beta2<-rmvnorm(1000,solve(A)%*%B2,solve(A))
like1<-rep(NA,1000)
like2<-rep(NA,1000)
for (i in 1:1000){
like1[i]<-dmvnorm(YY1,XX%*%beta1[i,],diag(1,100))
like2[i]<-dmvnorm(YY2,XX%*%beta2[i,],diag(1,100))
}
ev1[k]<-1/mean(1/like1)
ev2[k]<-1/mean(1/like2)
}
ev1[1]
ev2[1]
var(ev1)
var(ev2)






