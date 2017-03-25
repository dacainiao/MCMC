1.a,b
mar<- function(Y,X,betamean=rep(1,3),betavar=diag(1,3),beta=c(1,1,-1)){
k=dim(X)[2]
n=length(Y)
A<-t(X)%*%X+solve(betavar)
B<-t(X)%*%Y+solve(betavar)%*%betamean
K<-t(Y)%*%Y+t(betamean)%*%solve(betavar)%*%betamean
part<-(det(betavar%*%t(X)%*%X+diag(k)))^(-0.5)
out<-((2*pi)^(-n/2))*part*exp(-0.5*(K-t(B)%*%solve(A)%*%B))

like<-dmvnorm(YY1,XX%*%beta,diag(1,100))
prior<-dmvnorm(beta,betamean,betavar)
post<-dmvnorm(beta,solve(A)%*%B,solve(A))
out2<-like*prior/post
return(c(out,out2))
}
mar(YY1,XX)
mar(YY1,XX,betamean=rep(2,3),betavar=diag(2,3),beta=c(2,2,1))
mar(YY1,XX,betamean=c(3,1,4),betavar=diag(3,3),beta=c(1,1,1))


2.a
mar2<- function(Y,X,betamean=rep(0,3),betavar){
k=dim(X)[2]
n=length(Y)
A<-t(X)%*%X+solve(betavar)
B<-t(X)%*%Y+solve(betavar)%*%betamean
K<-t(Y)%*%Y+t(betamean)%*%solve(betavar)%*%betamean
part<-(det(betavar%*%t(X)%*%X+diag(k)))^(-0.5)
out<-((2*pi)^(-n/2))*part*exp(-0.5*(K-t(B)%*%solve(A)%*%B))
return(out)
}
data1ml<-rep(NA,5)
data2ml<-rep(NA,5)
for (i in 1:5){
data1ml[i]<-mar2(YY1,XX,betavar=diag(i,3))
data2ml[i]<-mar2(YY2,XX,betavar=diag(i,3))
}
ml<-rbind(data1ml,data2ml)
colnames(ml)<-c("model1","model2","model3","model4","model5")
rownames(ml)<-c("YY1","YY2")
print(ml)

2.b
data1post<-rep(NA,5)
data2post<-rep(NA,5)
for (i in 1:5){
data1post[i]<-data1ml[i]*0.2/sum(data1ml*0.2)
data2post[i]<-data2ml[i]*0.2/sum(data2ml*0.2)
}
postprob<-rbind(data1post,data2post)
colnames(postprob)<-c("model1","model2","model3","model4","model5")
rownames(postprob)<-c("YY1","YY2")
print(postprob)

2.c
inter<-function(Y,X){
betaout<-rep(NA,5)
for (i in 1:5){
betamean=rep(0,3)
betavar=diag(i,3)
A<-t(X)%*%X+solve(betavar)
B<-t(X)%*%Y+solve(betavar)%*%betamean
M<-solve(A)%*%B
V<-solve(A)
betaout[i]<-pnorm(0,M[1],V[1,1])
}
return(betaout)
}
sum(inter(YY1,XX)%*%data1post)
sum(inter(YY2,XX)%*%data2post)


3.a
mar2<- function(Y,X,betamean,betavar){
k=dim(X)[2]
n=length(Y)
A<-t(X)%*%X+solve(betavar)
B<-t(X)%*%Y+solve(betavar)%*%betamean
K<-t(Y)%*%Y+t(betamean)%*%solve(betavar)%*%betamean
part<-(det(betavar%*%t(X)%*%X+diag(k)))^(-0.5)
out<-((2*pi)^(-n/2))*part*exp(-0.5*(K-t(B)%*%solve(A)%*%B))
return(out)
}
ml1<-rep(NA,4)
ml2<-rep(NA,4)
ml1[1]<-mar2(Y=YY1,X=cbind(XX[,1]),betamean=0,betavar=1000)
ml2[1]<-mar2(Y=YY2,X=cbind(XX[,1]),betamean=0,betavar=1000)
ml1[2]<-mar2(Y=YY1,X=cbind(XX[,-3]),betamean=rep(0,2),betavar=diag(1000,2))
ml2[2]<-mar2(Y=YY2,X=cbind(XX[,-3]),betamean=rep(0,2),betavar=diag(1000,2))
ml1[3]<-mar2(Y=YY1,X=cbind(XX[,-2]),betamean=rep(0,2),betavar=diag(1000,2))
ml2[3]<-mar2(Y=YY2,X=cbind(XX[,-2]),betamean=rep(0,2),betavar=diag(1000,2))
ml1[4]<-mar2(Y=YY1,X=cbind(XX),betamean=rep(0,3),betavar=diag(1000,3))
ml2[4]<-mar2(Y=YY2,X=cbind(XX),betamean=rep(0,3),betavar=diag(1000,3))
ml<-rbind(ml1,ml2)
colnames(ml)<-c("model1","model2","model3","model4")
rownames(ml)<-c("YY1","YY2")
print(ml)


3.b
data1post<-rep(NA,4)
data2post<-rep(NA,4)
for (i in 1:4){
data1post[i]<-ml1[i]*0.25/sum(ml1*0.25)
data2post[i]<-ml2[i]*0.25/sum(ml2*0.25)
}
postprob<-rbind(data1post,data2post)
colnames(postprob)<-c("model1","model2","model3","model4")
rownames(postprob)<-c("YY1","YY2")
print(postprob)

3.c
data1post[2]+data1post[4]
data2post[2]+data2post[4]

