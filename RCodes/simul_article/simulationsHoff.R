library(BayesCanoFacto)
N <- 100; P <- c(10, 5, 15); H <- 5; K <- 3
sig_true <- rbind(c(.7, .5, .3,0,0), c(.5, .3, .7,0,0), c(.3, .7, .5,0,0))
orth <- function(n, h) qr.Q(qr(matrix(rnorm(n * h), n, h)))
sigma<-c(0.02,0.04,0.06,0.08,0.12,0.16,0.18)
U0 <- orth(N,H)
V <- vector("list", K)
X <- vector("list", K)
Xs <- vector("list", K)
normX<-rep(NA,K)
Sigma=list(sig_true[1,],sig_true[2,],sig_true[3,])
resHoff<-NULL
B<-500

t0<-Sys.time()
for(j in 1:length(sigma)){
for(b in 1:B){
  U0 <- orth(N,H)
  V <- vector("list", K)
  X <- vector("list", K)
  Xs <- vector("list", K)
  normX<-rep(NA,K)
  for(k in seq_len(K)) {
    E <- matrix(rnorm(N * P[k], 0, sigma[j]), nrow = N, ncol = P[k])
    V[[k]]<-orth(P[k],H)
    X[[k]] <- U0 %*% diag(Sigma[[k]]) %*% t(V[[k]]) + E
    X[[k]] <- scale(X[[k]], center = TRUE, scale = FALSE)   
    normX[k] <- norm(X[[k]], "F")                           
    Xs[[k]] <- X[[k]] 
  }
  Y<-NULL
 for(k in seq_len(K)) Y<-cbind(Y,Xs[[k]]) 
fit <- hoff_svd(Y, R = 5, n_iter = 1000, burn = 500, center = TRUE)

resHoff<-rbind(resHoff,c(fit$post,sigma[j]))
}
}
t1<-Sys.time()

save(resHoff,file="resHoff.RData")

