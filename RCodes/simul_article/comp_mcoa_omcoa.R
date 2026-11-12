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
      Xs[[k]] <- X[[k]] 
    }

omcoa_als <- function(Xs, R = default_rank(Xs), maxit = 2000, tol = 1e-10) {
  K <- length(Xs)
  X <- do.call(cbind, Xs)
  U <- svd(X, nu = R, nv = 0)$u
  V <- lapply(Xs, function(Xk) procrustes(crossprod(Xk, U)))
  S <- lapply(seq_len(K), function(k) diag(t(U) %*% Xs[[k]] %*% V[[k]]))
  crit <- function() sum(vapply(seq_len(K), function(k)
    sum((Xs[[k]] - U %*% diag(S[[k]], R) %*% t(V[[k]]))^2), numeric(1)))
  old <- Inf
  iter <- 0
  for (it in seq_len(maxit)) {
    iter <- it
    Hm <- do.call(cbind, lapply(seq_len(K), function(k) diag(S[[k]], R) %*% t(V[[k]])))
    U <- procrustes(X %*% t(Hm))
    V <- lapply(seq_len(K), function(k) procrustes(crossprod(Xs[[k]], U) %*% diag(S[[k]], R)))
    S <- lapply(seq_len(K), function(k) diag(t(U) %*% Xs[[k]] %*% V[[k]]))
    f <- crit()
    if (old - f < tol * f) break
    old <- f
  }
  for (k in seq_len(K)) {
    sg <- sign(S[[k]])
    sg[sg == 0] <- 1
    V[[k]] <- sweep(V[[k]], 2, sg, "*")
    S[[k]] <- S[[k]] * sg
  }
  names(V) <- names(S) <- names(Xs)
  list(U = U, V = V, S = S, crit = f, iter = iter)
}    




#### omcoa ####
m2<-omcoa_als(Xs,H)  
m2$S
t(as.matrix(m2$U))%*%as.matrix(U0)


### mcoa with ade4 ####
Xt=list()
K<-length(X)
for(k in 1:K) Xt[[k]]=as.data.frame(X[[k]])
XX=ade4::ktab.list.df(Xt)
res=ade4::mcoa(XX, scannf = FALSE, nf = H,option = "uniform")

U=as.matrix(res$SynVar)

for(j in 1:ncol(U)) U[,j]/norm(as.matrix(U[,j]),"F")

## cos de l'angle  
t(as.matrix(U))%*%as.matrix(U0)


V<-list()
SIG<-list()
for(k in 1:K){
  V[[k]]=t(X[[k]])%*%U
  SIG[[k]]=sqrt(diag(t(U)%*%X[[k]]%*%V[[k]]))
}

SIG
