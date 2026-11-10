#remotes::install_github("GALHARRET/BayesCanoFacto")
library(BayesCanoFacto)
B<-250

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

 

  m2<-omcoa_als(Xs,R)
  
  dr <- gibbs_sparse_omcoa(X, R, n_iter = 2000, burn = 1000, init = m2,lambda=1)
  dr <- align_draws(dr,U0)
  resGSOMCOA1[[b]]<-list(summary=summary_sigma(dr,sparse=T),rank=rank_posterior(dr))

  dr <- gibbs_sparse_omcoa(X, R, n_iter = 2000, burn = 1000, init = m2,
                           slab = "noise")
  dr <- align_draws(dr, U0)
  resGSOMCOA2[[b]]<-list(summary=summary_sigma(dr,sparse=T),rank=rank_posterior(dr))

  dr <- gibbs_sparse_omcoa(X, R, n_iter = 2000, burn = 1000, init = m2,
                           slab = "noise", g = 10 * 99 * c(10, 5, 15))
  dr <- align_draws(dr, U0)
  resGSOMCOA3[[b]]<-list(summary=summary_sigma(dr,sparse=T),rank=rank_posterior(dr))
}
t1<-Sys.time()

res<-list(gsomcoa1=resGSOMCOA1,
           gsomcoa2=resGSOMCOA2,
           gsomcoa3=resGSOMCOA3)

save(res,file=paste("~/Dropbox/res",j,".l_vs_g.RData",sep=""))
}
t1<-Sys.time()

