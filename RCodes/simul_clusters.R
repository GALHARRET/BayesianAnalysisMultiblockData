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

orth <- function(n, h, center = FALSE) {
  M <- matrix(rnorm(n * h), n, h)
  if (center) M <- scale(M, scale = FALSE)
  qr.Q(qr(M))
}

U<-orth(3*n,R,center=T)


###  First scenario: same structure for the three blocks 
R=3
D<-list()


Sigma_t<-function(t){
  Sigma0<-list(c(0.7,0.5,0.3),c(0.7,0.5,0.3),c(0.7,0.5,0.3))
  Sigma1<-list(c(0.7,0.5,0.3),c(0.3,0.7,0.5),c(0.5,0.3,0.7))
  Sigma<-list()
  for(k in 1:length(Sigma0)) Sigma[[k]]<-(1-t)*Sigma0[[k]]+t*Sigma1[[k]]
  return(Sigma)
}

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
    Yh[[k]]<-U%*%diag(Sigma_t(1)[[k]])%*%t(V[[k]])
    EPS[[k]]<-matrix(rnorm(3*n*p[[k]],sd=sqrt(alpha[k]*sum(Sigma_t(1)[[k]]^2))),ncol=p[[k]])/(3*n)
    Y[[k]]<-Yh[[k]]+EPS[[k]]
  }
  
  res<-gibbs_omcoa(Y,H =3,n_iter = 6000,burn = 1000) 
  apply(res$S[,,1],2,mean)
  apply(res$S[,,2],2,mean)
  apply(res$S[,,3],2,mean)
  
  
  
  for(k in 1:K){for(k inmean 1:K){
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

