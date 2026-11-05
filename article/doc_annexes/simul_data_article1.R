library(mvtnorm)
setwd("C:/Users/j.galharret/Documents/GitHub/OMCOA/Rcodes/")
source("algoOMCOA.R")
source("algoMCOA.R")
source("algoMCOAv.R")
d<-9
sig<-1
SIG<-diag(rep(sig,d))

n<-50
x1<-rmvnorm(n,mean = rep(3,d),sigma = SIG)

x2<-rmvnorm(n,mean = rep(-3,d),sigma = SIG)

x3<-rmvnorm(n,mean = c(-3,rep(3,d-1)),sigma = SIG)



X<-rbind(x1,x2,x3)

N<-dim(X)[1]

U<-svd(X)$u

###  First scenario: same structure for the three blocks 
R=3
U1=list(U[,1:R],U[,1:R],U[,1:R])
D<-list()
D[[1]]<-diag(c(0.8,0.22,0))
D[[2]]<-diag(c(0,0.2,0.8))
D[[3]]<-diag(c(0.2,0,0.8))

V<-list()
p<-list()
p[[1]]<-10
V[[1]]<-svd(matrix(rnorm(n*p[[1]]),ncol=p[[1]]),nu=R,nv=R)$v
p[[2]]<-15
V[[2]]<-svd(matrix(rnorm(n*p[[2]]),ncol=p[[2]]),nu=R,nv=R)$v
p[[3]]<-20
V[[3]]<-svd(matrix(rnorm(n*p[[3]]),ncol=p[[3]]),nu=R,nv=R)$v


Y<-Yh<-EPS<-list()
alpha<-c(0.1,0.05,0.15)

B=5000
r=2
K=3
MCOA<-NULL
OMCOA<-NULL
Block<-NULL
Dim<-NULL

for(b in 1:B){
  for(k in 1:3){
    Yh[[k]]<-U1[[k]]%*%D[[k]]%*%t(V[[k]])
    EPS[[k]]<-matrix(rnorm(3*n*p[[k]],sd=sqrt(alpha[k]*sum(D[[k]]^2))),ncol=p[[k]])/(3*n)
    Y[[k]]<-Yh[[k]]+EPS[[k]]
  }
  
  res1<-algoMCOAv(Y,R=r)
  res2<-algoOMCOA(Y,R=r)
  for(k in 1:K){
    Block<-c(Block,rep(k,r))
    Dim<-c(Dim,1:r)}
  MCOA<-c(MCOA,unlist(res1$SIG))
  OMCOA<-c(OMCOA,unlist(res2$SIG))
}
res=as.data.frame(rbind(
  cbind(Block,Dim,OMCOA),cbind(Block,Dim,MCOA)))

res$method<-c(rep("OMCOA",B),rep("MCOA",B))

save(res,file = "res1.RData")
colnames(res)[3]="value"

res |>
  group_by(Block,Dim,method) |>
  summarise(M=mean(value),SD=sd(value))



# SECOND SCENARIO

r=3
K=3
MCOA<-NULL
OMCOA<-NULL
Block<-NULL
Dim<-NULL


U1=list(U[,1:R],U[,(R+1):(R+3)],U[,(R+4):(R+6)])
D<-list()
D[[1]]<-diag(c(0.8,0.2,0))
D[[2]]<-diag(c(0,0.2,0.8))
D[[3]]<-diag(c(0.2,0,0.8))

V<-list()
p<-list()
p[[1]]<-10
V[[1]]<-svd(matrix(rnorm(n*p[[1]]),ncol=p[[1]]),nu=R,nv=R)$v
p[[2]]<-15
V[[2]]<-svd(matrix(rnorm(n*p[[2]]),ncol=p[[2]]),nu=R,nv=R)$v
p[[3]]<-20
V[[3]]<-svd(matrix(rnorm(n*p[[3]]),ncol=p[[3]]),nu=R,nv=R)$v



for(b in 1:B){
  for(k in 1:3){
    Yh[[k]]<-U1[[k]]%*%D[[k]]%*%t(V[[k]])
    EPS[[k]]<-matrix(rnorm(3*n*p[[k]],sd=sqrt(alpha[k]*sum(D[[k]]^2))),ncol=p[[k]])/(3*n)
    Y[[k]]<-Yh[[k]]+EPS[[k]]
  }
  
  res1<-algoMCOAv(Y,R=r)
  res2<-algoOMCOA(Y,R=r)
  for(k in 1:K){
    Block<-c(Block,rep(k,r))
    Dim<-c(Dim,1:r)}
  MCOA<-c(MCOA,unlist(res1$SIG))
  OMCOA<-c(OMCOA,unlist(res2$SIG))
}
res=as.data.frame(rbind(
  cbind(Block,Dim,OMCOA),cbind(Block,Dim,MCOA)))

res$method<-c(rep("OMCOA",B),rep("MCOA",B))

save(res,file = "res2.RData")
colnames(res)[3]="value"

res |>
  group_by(Block,Dim,method) |>
  summarise(M=mean(value),SD=sd(value))



D<-list()
D[[1]]<-diag(c(0.5,0.22,0,0,0.28))
D[[2]]<-diag(c(0,0.22,0.5,0.14,0.14))
D[[3]]<-diag(c(0,0.22,0.14,0.5,0.14))

V<-list()
p<-list()
p[[1]]<-10
V[[1]]<-svd(matrix(rnorm(n*p[[1]]),ncol=p[[1]]),nu=d,nv=d)$v
p[[2]]<-15
V[[2]]<-svd(matrix(rnorm(n*p[[2]]),ncol=p[[2]]),nu=d,nv=d)$v
p[[3]]<-20
V[[3]]<-svd(matrix(rnorm(n*p[[3]]),ncol=p[[3]]),nu=d,nv=d)$v

Y<-Yh<-EPS<-list()

for(k in 1:3){
  Yh[[k]]<-U%*%D[[k]]%*%t(V[[k]])
  EPS[[k]]<-matrix(rnorm(3*n*p[[k]]),ncol=p[[k]])/(3*n)
  Y[[k]]<-Yh[[k]]+EPS[[k]]
}

res1<-algoOMCOA(Y,R=1)
res2<-algoMCOA(Y,R=1)
