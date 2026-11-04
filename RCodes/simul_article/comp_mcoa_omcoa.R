path<-getwd()

setwd(path)

summary_sigma_freq<-function(out){
  K<-dim(out)[1]
  H<-dim(out)[2]
  block<-NULL
  for(k in 1:K) block<-c(block,rep(k,H))
  SIG<-as.vector(t(out))
  res<-as.data.frame(matrix(NA,ncol=3,nrow=length(SIG)))
  colnames(res)<-c("block","dim","estimate")
  res$block<-block
  res$dim<-rep(1:H,K)
  res$estimate<-SIG
  return(res)
}
source("~/Documents/GitHub/BayesianAnalysisMultiblockData/RCodes/align_solution.R")
 
library(BayesCanoFacto)

N <- 100; P <- c(10, 5, 15); H <- 3; K <- 3
sig_true <- rbind(c(.7, .5,.3), c(.5, .3,.7), c(.3,.7,.5))
orth <- function(n, h) qr.Q(qr(matrix(rnorm(n * h), n, h)))
sigma<-c(0.02,0.04,0.06,0.08,0.12,0.16,0.18,0.20,0.22)
Sigma=list(sig_true[1,],sig_true[2,],sig_true[3,])
orth <- function(n, h) qr.Q(qr(scale(matrix(rnorm(n * h), n, h), scale = FALSE)))

U0 <- orth(N,H)
V <- vector("list", K)
X <- vector("list", K)
Xs <- vector("list", K)
normX<-rep(NA,K)
j=1
for(k in seq_len(K)){
      E <- matrix(rnorm(N * P[k], 0, sigma[j]), nrow = N, ncol = P[k])
      V[[k]]<-orth(P[k],H)
      X[[k]] <- U0 %*% diag(Sigma[[k]]) %*% t(V[[k]]) + E
      X[[k]] <- scale(X[[k]], center = TRUE, scale = FALSE)   
      normX[k] <- norm(X[[k]], "F")                           
      Xs[[k]] <- X[[k]] / normX[k]
    }
    
m0 <- mcoa(Xs,H)
apply(m0$U,2,function(x) norm(as.matrix(x),"F"))
SIG<-align_solution(U=m0$U,V = m0$V,sigma = m0$S,U_ref=U0)$sigma
for(k in 1:K) SIG[k,]<-SIG[k,]*normX[k]
summary_sigma_freq(SIG)
  

m2<-omcoa_als(Xs,R = 5)
SIG<-align_solution(U=m2$U,V = m2$V,sigma = m2$S,U_ref=U0)$sigma
for(k in 1:K) SIG[k,]<-SIG[k,]*normX[k]
summary_sigma_freq(SIG)
  
dr <- gibbs_omcoa(Xs, H, n_iter = 10000, burn = 5000,init=m2)
dr <- align_draws(dr, U0)
for(i in 1:1000){
  for (k in 1:K) dr$S[i,k,]<-dr$S[i,k,]*normX[k]
} 
summary_sigma(dr)
  
dr <- gibbs_sparse_omcoa(Xs,n_iter = 10000, burn = 5000, init = m2,R = 5)
dr <- align_sparse(dr, U0)

Reff<-apply(dr$G,1,function(g) sum(colSums(g)>0))
hist(Reff)


post_Reff <- table(factor(Reff, levels = 0:5)) / length(Reff)   # P(Reff = j | X)
mode_Reff <- as.integer(names(which.max(post_Reff)))
plot(Reff, type = "s")  
rank_posterior(dr)

colSums(dr$G[1,,])

list(summary=summary_sigma(dr,sparse=T),rank=rank_posterior(dr))
  
  dr <- gibbs_sparse_omcoa(Xs, H, n_iter = 2000, burn = 1000, init = m2,lambda = 0.1)
  dr <- align_sparse(dr, U0)
  resGSOMCOA2[[b]]<-list(summary=summary_sigma(dr,sparse=T),rank=rank_posterior(dr))
  
  dr <- gibbs_sparse_omcoa(Xs, H, n_iter = 2000, burn = 1000, init = m2,lambda = 10)
  dr <- align_sparse(dr, U0)
  resGSOMCOA3[[b]]<-list(summary=summary_sigma(dr,sparse=T),rank=rank_posterior(dr))
}


res<-list(mcoa=resMCOA,mcoa2=resMCOA2,omcoa=resOMCOA,
           gomcoa=resGOMCOA,
           gsomcoa1=resGSOMCOA1,
           gsomcoa2=resGSOMCOA2,
           gsomcoa3=resGSOMCOA3)

save(res,file=paste("res",k,".RData",sep=""))}
