#remotes::install_github("GALHARRET/BayesCanoFacto")
library(BayesCanoFacto)
B<-1
B<-500

resGSOMCOA1<-list()
resGSOMCOA2<-list()
resGSOMCOA3<-list()

N <- 100; P <- c(10, 5, 15); R <- 4; K <- 3
sig_true <- rbind(c(.7, .5,0 ,0), c(.5, .7, 0,0), c(.3,0, .5,0))
Sigma=list(sig_true[1,],sig_true[2,],sig_true[3,])
orth <- function(n, h) qr.Q(qr(scale(matrix(rnorm(n * h), n, h), scale = FALSE)))
sigma<-c(0.02,0.04,0.06,0.08,0.12,0.16,0.18)

t0<-Sys.time()
for(j in 1:length(sigma)){
for(b in 1:B){
  U0 <- orth(N,R)
  V <- vector("list", K)
  X <- vector("list", K)
  Xs <- vector("list", K)
  normX<-rep(NA,K)
  for(k in seq_len(K)) {
    E <- matrix(rnorm(N * P[k], 0, sigma[j]), nrow = N, ncol = P[k])
    V[[k]]<-orth(P[k],R)
    X[[k]] <- U0 %*% diag(Sigma[[k]]) %*% t(V[[k]]) + E
    X[[k]] <- scale(X[[k]], center = TRUE, scale = FALSE)
    normX[k] <- norm(X[[k]], "F")
    Xs[[k]] <- X[[k]] #/ normX[k]
  }

  m0 <- mcoa(Xs,R)

  SIG<-align_solution(U=m0$U,U_ref=U0,S=m0$S)$S
  #for(k in 1:K) SIG[k,]<-SIG[k,]#*normX[k]
  resMCOA[[b]]<-summary_sigma(SIG)


  m2<-omcoa_als(Xs,R)
  SIG<-align_solution(U=m2$U,U_ref=U0,S=m2$S)$S
  #for(k in 1:K) SIG[k,]<-SIG[k,]#*normX[k]
  resOMCOA[[b]]<-summary_sigma(SIG)

  dr <- gibbs_omcoa(X, R, n_iter = 2000, burn = 1000, init = m2)
  dr <- align_draws(dr, U0)
  resGOMCOA[[b]]<-summary_sigma(dr)

  dr <- gibbs_sparse_omcoa(X, R, n_iter = 2000, burn = 1000, init = m2,lambda=0.1)
  dr <- align_draws(dr,U0)
  resGSOMCOA1[[b]]<-list(summary=summary_sigma(dr,sparse=T),rank=rank_posterior(dr))

  dr <- gibbs_sparse_omcoa(X, R, n_iter = 2000, burn = 1000, init = m2,lambda=1)
  dr <- align_draws(dr, U0)
  resGSOMCOA2[[b]]<-list(summary=summary_sigma(dr,sparse=T),rank=rank_posterior(dr))

  dr <- gibbs_sparse_omcoa(X, R, n_iter = 2000, burn = 1000, init = m2,lambda = 10)
  dr <- align_draws(dr, U0)
  resGSOMCOA3[[b]]<-list(summary=summary_sigma(dr,sparse=T),rank=rank_posterior(dr))
}
t1<-Sys.time()

res<-list(mcoa=resMCOA,omcoa=resOMCOA,
           gomcoa=resGOMCOA,
           gsomcoa1=resGSOMCOA1,
           gsomcoa2=resGSOMCOA2,
           gsomcoa3=resGSOMCOA3)

save(res,file=paste("res",j,".RData",sep=""))
}


