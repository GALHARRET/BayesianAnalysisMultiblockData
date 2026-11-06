library(BayesCanoFacto)

N <- 100; P <- c(10, 10, 10); H <- 3; K <- 3
sig_true <- rbind(c(.7, .5,.3), c(.5, .3,.7), c(.3,.7,.5))
orth <- function(n, h) qr.Q(qr(matrix(rnorm(n * h), n, h)))
Sigma=list(sig_true[1,],sig_true[2,],sig_true[3,])
orth <- function(n, h) qr.Q(qr(scale(matrix(rnorm(n * h), n, h), scale = FALSE)))

U0 <- orth(N,H)
V <- vector("list", K)
X <- vector("list", K)
Xs <- vector("list", K)
normX<-rep(NA,K)
sigma<-0
for(k in seq_len(K)){
      E <- matrix(rnorm(N * P[k], 0, sigma), nrow = N, ncol = P[k])
      V[[k]]<-orth(P[k],H)
      X[[k]] <- U0 %*% diag(Sigma[[k]]) %*% t(V[[k]]) + E
      X[[k]] <- scale(X[[k]], center = TRUE, scale = FALSE)   
      normX[k] <- norm(X[[k]], "F")                           
      Xs[[k]] <- X[[k]] / normX[k]
    }

    
### mcoa with your function ###
m0 <- BayesCanoFacto::mcoa(Xs,H)
apply(m0$V[[1]],2,function(x) norm(as.matrix(x),"F"))
align<-align_solution(U=m0$U,V = m0$V,S = m0$S,U_ref=U0)
SIG<-align$S
align$cosines
for(k in 1:K) SIG[k,]<-SIG[k,]*normX[k]
summary_sigma(SIG)


#### omcoa ####
m2<-omcoa_als(Xs,H)  
align<-align_solution(U=m2$U,V = m2$V,S = m2$S,U_ref=U0)
SIG<-align$S
align$cosines
for(k in 1:K) SIG[k,]<-SIG[k,]*normX[k]
summary_sigma(SIG)


### mcoa with ade4 ####
Xt=list()
K<-length(X)
for(k in 1:K) Xt[[k]]=as.data.frame(X[[k]])
XX=ade4::ktab.list.df(Xt)
res=ade4::mcoa(XX, scannf = FALSE, nf = H,option = "uniform")
U=as.matrix(res$SynVar)
V<-list()
SIG<-list()
for(k in 1:K){
  V[[k]]=t(X[[k]])%*%U
  SIG[[k]]=sqrt(diag(t(U)%*%X[[k]]%*%V[[k]]))
}
m0<-list(U=U,S=SIG,V=V)
align<-align_solution(U=m0$U,V = m0$V,S = m0$S,U_ref=U0)
SIG<-align$S
align$cosines

