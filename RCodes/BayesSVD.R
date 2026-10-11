rang_max(X){
  max_r<-dim()
  for(k in 1:length(X)) max_r
}


sparse_bayes_canonical_factorisation <- function(X, R = rang_max(X), niter = 5000, burnin = 1000, thin = 4,
                                   a0 = 1e-3, b0 = 1e-3,
                                   a_pi = 1, b_pi = 1,
                                   spike_prec = 1e6,
                                   c0 = 1e-3, d0 = 1e-3,
                                   verbose = TRUE, seed = 123) {
  concat=function(X){
    K=length(X)
    concat=NULL
    for(k in 1:K) concat=cbind(concat,X[[k]])
    return(concat)
  }
  
  K=length(X)
  Y=concat(X)
  n<-dim(X[[1]])[1]
  set.seed(seed)
  n <- nrow(Y); p <- ncol(Y)
  # initialize
  W <- matrix(rnorm(n * R, sd = 0.1), n, R)
  H <- matrix(rnorm(p * R, sd = 0.1), p, R)
  alpha <- rep(1, R)
  pi_k <- rep(0.5, R)
  Z <- matrix(rbinom(p * R, 1, 0.5), p, R)  # inclusion indicators
  tau <- 1.0

  
  ## Results
  U=array(NA,dim=c(niter,n,R))
  SIGk=list()
  Vk=list()
  for(k in 1:K){
    Vk[[k]]=array(NA,dim=c(niter,dim(X[[k]])[2],R))
    SIGk[[k]]=matrix(NA,niter,R)
  } 
  Zs=array(NA,dim=c(niter,p,R))
  pi=matrix(NA,nrow=niter,ncol=R)
  alpha_s=matrix(NA,nrow=niter,ncol=R)
  tau_s=rep(NA,niter)
  
  for (iter in 1:niter) {
    # --- sample W (rows) given H, alpha, tau
    S_H <- crossprod(H)                     
    Precision_w <- tau * S_H + diag(alpha)  
    Rw <- chol(Precision_w)                 
    Sigma_w <- chol2inv(Rw)
    # compute means for all rows: T = tau * Sigma_w %*% t(Y %*% H)
    Tmat <- tau * Sigma_w %*% t(Y %*% H)
    Zmat <- matrix(rnorm(R * n), nrow = R)
    RinvZ <- backsolve(Rw, Zmat)            
    Wt <- Tmat + RinvZ
    W <- t(Wt)                              
    
    # --- sample H (rows) given W, alpha, tau, Z (sparse prior)
    S_W <- crossprod(W)                      # K x K
    # sample each row h_j (length K) separately because prior precision diagonal depends on Z[j,]
    for (j in 1:p) {
      prior_prec_vec <- ifelse(Z[j, ] == 1, alpha, spike_prec) # length K
      Precision_h <- tau * S_W + diag(prior_prec_vec)          # K x K
      Rh <- chol(Precision_h)
      Sigma_h <- chol2inv(Rh)
      Tmat_h <- tau * Sigma_h %*% (t(W) %*% Y[, j])            # K x 1
      z_h <- rnorm(R)
      Rinv_z <- backsolve(Rh, z_h)
      hj <- as.vector(Tmat_h + Rinv_z)
      H[j, ] <- hj
    }
    
    # --- sample Z (inclusion indicators) given H and pi_k
    # P(z_jk = 1 | h_jk) ∝ pi_k * N(h_jk | 0, var_slab=1/alpha_k)
    #                       (1-pi_k) * N(h_jk | 0, var_spike = 1/spike_prec)
    for (k in 1:R) {
      # compute log-probabilities vectorized over j
      log_p1 <- log(pi_k[k]) + dnorm(H[, k], mean = 0, sd = sqrt(1 / alpha[k]), log = TRUE)
      log_p0 <- log(1 - pi_k[k]) + dnorm(H[, k], mean = 0, sd = sqrt(1 / spike_prec), log = TRUE)
      # stabilize
      mlog <- pmax(log_p1, log_p0)
      p1 <- exp(log_p1 - mlog) / (exp(log_p1 - mlog) + exp(log_p0 - mlog))
      Z[, k] <- rbinom(p, size = 1, prob = p1)
      # if Z=0 force H small? we already sample H under spike precision; optionally set to exactly 0:
      # H[Z[,k] == 0, k] <- 0  # optional hard-thresholding; commented out to keep conjugacy
    }
    Zs[iter,,]=Z
    # --- sample pi_k (Beta prior)
    for (k in 1:R) {
      sumz <- sum(Z[, k])
      pi_k[k] <- rbeta(1, a_pi + sumz, b_pi + p - sumz)
    }
    pi[iter,]=pi_k
    # --- sample alpha_k (ARD) using only W and slab entries of H
    for (k in 1:R) {
      # slab entries are those with Z==1
      Jk <- sum(Z[, k])
      shape <- a0 + n / 2 + Jk / 2
      # sum squares: W column + only H entries where Z==1
      ss_w <- sum(W[, k]^2)
      ss_h_slab <- sum(H[Z[, k] == 1, k]^2)
      rate <- b0 + 0.5 * (ss_w + ss_h_slab)
      alpha[k] <- rgamma(1, shape = shape, rate = rate)
    }
    alpha_s[iter,]=alpha
    
    # --- sample tau (noise precision)
    E <- Y - W %*% t(H)
    ressq <- sum(E^2)
    shape_t <- c0 + (n * p) / 2
    rate_t <- d0 + 0.5 * ressq
    tau <- rgamma(1, shape = shape_t, rate = rate_t)
    tau_s[iter]<-tau
    
    # --- save posterior summaries
    Xhat=W %*% t(H)
    S=svd(Xhat,nu=R,nv=R)
    U[iter,,]=S$u
    for(k in 1:K){
      S1=svd(t(X[[k]])%*%U[iter,,],nu=R,nv=R)
      Vk[[k]][iter,,]=(S1$u)%*%t(S1$v)
      SIGk[[k]][iter,]=diag(t(U[iter,,])%*%X[[k]]%*%Vk[[k]][iter,,])
    }
    
    if (verbose && (iter %% max(1, floor(niter/10)) == 0))
      cat(sprintf("Iter %d / %d  (alpha[1]=%.3g, pi[1]=%.3g, tau=%.3g)\n",
                  iter, niter, alpha[1], pi_k[1], tau))
  }
  
  
  list(U = U,
       SIGk = SIGk,
       Vk = Vk,
       tau = tau_s,
       Z=Zs,
       alpha=alpha_s,
       pi=pi
       )
}




# Installation si nécessaire
if (!requireNamespace("rstiefel", quietly = TRUE)) {
  install.packages("rstiefel")
}
library(rstiefel)

bsvd_vMF <- function(Y, r,
                     niter = 5000, burn = 1000, thin = 1,
                     a0 = 1, b0 = 1, tau2 = 1) {
  # Y: matrice n x p
  n <- nrow(Y); p <- ncol(Y)
  iters_kept <- floor((niter - burn) / thin)
  if (iters_kept <= 0) stop("Augmente niter ou réduis burn/thin")
  
  # initialisation par SVD
  s0 <- svd(Y)
  U <- s0$u[, 1:r, drop = FALSE]
  V <- s0$v[, 1:r, drop = FALSE]
  d <- s0$d[1:r]
  sigma2 <- sum((Y - U %*% diag(d) %*% t(V))^2) / (n * p)
  
  # stockage des échantillons
  S_U <- array(0, dim = c(n, r, iters_kept))
  S_V <- array(0, dim = c(p, r, iters_kept))
  S_D <- matrix(0, nrow = iters_kept, ncol = r)
  S_sigma2 <- numeric(iters_kept)
  
  keep_idx <- 0
  for (iter in 1:niter) {
    # --- échantillonnage de U ---
    # paramètre matriciel A_U = (Y %*% V) %*% diag(d) / sigma2
    A_U <- (Y %*% V) %*% diag(d) / sigma2
    # fonction d'échantillonnage matrix-vMF dans rstiefel :
    # selon la version la fonction peut être rmf.matrix(A, r) ou rmf.matrix(A)
    # On tente l'appel courant : rstiefel::rmf.matrix(A_U)
    U <- rstiefel::rmf.matrix(A_U)  # adapter si la signature diffère
    
    # --- échantillonnage de V ---
    A_V <- (t(Y) %*% U) %*% diag(d) / sigma2
    V <- rstiefel::rmf.matrix(A_V)
    
    # --- échantillonnage des singular values d_j ---
    a_j <- numeric(r)
    for (j in 1:r) a_j[j] <- crossprod(U[, j], Y %*% V[, j])
    var_d <- 1 / (1 / sigma2 + 1 / tau2)
    mean_d <- var_d * (a_j / sigma2)
    d <- rnorm(r, mean = mean_d, sd = sqrt(var_d))
    
    # --- échantillonnage de sigma2 ---
    res <- Y - U %*% diag(d) %*% t(V)
    shape_post <- a0 + (n * p) / 2
    scale_post <- b0 + sum(res^2) / 2
    # sigma2 ~ Inv-Gamma(shape_post, scale_post) => 1/rgamma(shape_post, rate = scale_post)
    sigma2 <- 1 / rgamma(1, shape = shape_post, rate = scale_post)
    
    # stockage si beyond burn et selon thin
    if (iter > burn && ((iter - burn) %% thin == 0)) {
      keep_idx <- keep_idx + 1
      S_U[, , keep_idx] <- U
      S_V[, , keep_idx] <- V
      S_D[keep_idx, ] <- d
      S_sigma2[keep_idx] <- sigma2
    }
    
    # (optionnel) message de progression
    if (iter %% 500 == 0) cat("iter:", iter, " sigma2:", round(sigma2,4), "\n")
  }
  
  list(U = S_U, V = S_V, D = S_D, sigma2 = S_sigma2,
       posterior_mean = list(U = apply(S_U, c(1,2), mean),
                             V = apply(S_V, c(1,2), mean),
                             D = colMeans(S_D),
                             sigma2 = mean(S_sigma2)))
}
