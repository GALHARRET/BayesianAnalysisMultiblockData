# Bayesian SVD / Bayesian low-rank matrix factorization with ARD (Gibbs)
bayes_svd_gibbs <- function(Y, K = 10, niter = 2000, burnin = 1000, thin = 5,
                            a0 = 1e-3, b0 = 1e-3, c0 = 1e-3, d0 = 1e-3,
                            verbose = TRUE) {
  # Y: n x p data matrix
  # K: max rank (components)
  # returns list with posterior mean Xhat, mean alpha, mean tau, svals_mat (if saved)
  n <- nrow(Y); p <- ncol(Y)
  # initialize
  set.seed(123)
  W <- matrix(rnorm(n * K, sd = 0.1), n, K)
  H <- matrix(rnorm(p * K, sd = 0.1), p, K)
  alpha <- rep(1, K)
  tau <- 1.0
  
  # storage
  nsave <- floor((niter - burnin) / thin)
  if (nsave < 1) nsave <- 0
  Xsum <- matrix(0, n, p)
  alpha_sum <- rep(0, K)
  tau_sum <- 0
  svals_mat <- matrix(NA, nrow = nsave, ncol = K)
  save_i <- 0
  
  # pre-allocate small matrices
  for (iter in 1:niter) {
    # --- sample W (rows) given H, alpha, tau
    S_H <- crossprod(H)                      # K x K
    Precision_w <- tau * S_H + diag(alpha)  # K x K
    Rw <- chol(Precision_w)                  # upper tri such that Precision = t(Rw) %*% Rw
    Sigma_w <- chol2inv(Rw)                  # K x K
    # compute means for all rows: T = tau * Sigma_w %*% t(Y %*% H)  (K x n)
    Tmat <- tau * Sigma_w %*% t(Y %*% H)
    Z <- matrix(rnorm(K * n), nrow = K)
    RinvZ <- backsolve(Rw, Z)                # K x n
    Wt <- Tmat + RinvZ
    W <- t(Wt)                               # n x K
    
    # --- sample H (rows) given W, alpha, tau
    S_W <- crossprod(W)                      # K x K
    Precision_h <- tau * S_W + diag(alpha)
    Rh <- chol(Precision_h)
    Sigma_h <- chol2inv(Rh)
    Tmat_h <- tau * Sigma_h %*% (t(W) %*% Y) # K x p
    Z_h <- matrix(rnorm(K * p), nrow = K)
    RinvZ_h <- backsolve(Rh, Z_h)            # K x p
    Ht <- Tmat_h + RinvZ_h
    H <- t(Ht)                               # p x K
    
    # --- sample alpha_k (ARD)
    for (k in 1:K) {
      ss <- sum(W[, k]^2) + sum(H[, k]^2)
      shape <- a0 + (n + p) / 2
      rate <- b0 + 0.5 * ss
      alpha[k] <- rgamma(1, shape = shape, rate = rate)
    }
    
    # --- sample tau (noise precision)
    E <- Y - W %*% t(H)
    ressq <- sum(E^2)
    shape_t <- c0 + (n * p) / 2
    rate_t <- d0 + 0.5 * ressq
    tau <- rgamma(1, shape = shape_t, rate = rate_t)
    
    # --- save posterior summaries
    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      save_i <- save_i + 1
      Xcur <- W %*% t(H)
      Xsum <- Xsum + Xcur
      alpha_sum <- alpha_sum + alpha
      tau_sum <- tau_sum + tau
      # save singular values of current X (thin optional)
      svals <- svd(Xcur, nu = 0, nv = 0)$d
      m <- length(svals)
      if (m < K) svals <- c(svals, rep(0, K - m))
      svals_mat[save_i, ] <- svals[1:K]
    }
    
    if (verbose && (iter %% max(1, floor(niter/10)) == 0))
      cat(sprintf("Iter %d / %d  (alpha[1]=%.3g, tau=%.3g)\n", iter, niter, alpha[1], tau))
  }
  
  if (nsave > 0) {
    Xhat <- Xsum / nsave
    alpha_mean <- alpha_sum / nsave
    tau_mean <- tau_sum / nsave
  } else {
    Xhat <- W %*% t(H)
    alpha_mean <- alpha
    tau_mean <- tau
  }
  
  res <- list(Xhat = Xhat,
              alpha_mean = alpha_mean,
              tau_mean = tau_mean,
              svals_mat = if (nsave > 0) svals_mat else NULL,
              last_W = W,
              last_H = H)
  return(res)
}

# ---------------------------
# Small synthetic example
set.seed(2026)
n <- 60; p <- 40; r_true <- 3
U_true <- matrix(rnorm(n * r_true), n, r_true)
V_true <- matrix(rnorm(p * r_true), p, r_true)
D_true <- diag(c(5, 2.5, 1.2))
X_true <- U_true %*% D_true %*% t(V_true)    # true low-rank matrix
sigma <- 0.5

Y <- X_true #+ matrix(rnorm(n * p, sd = sigma), n, p)

# Run sampler (short run for demo)
res <- bayes_svd_gibbs(Y, K = 3, niter = 1200, burnin = 600, thin = 6,
                       a0 = 1e-3, b0 = 1e-3, c0 = 1e-3, d0 = 1e-3,
                       verbose = TRUE)

# Posterior mean estimate
Xhat <- res$Xhat
cat("Posterior mean tau:", res$tau_mean, "\n")
cat("Posterior mean alpha:", round(res$alpha_mean, 3), "\n")

# SVD of posterior mean
svd_hat <- svd(Xhat)
cat("Top singular values (posterior mean X):\n")
print(round(svd_hat$d[1:3], 3))

# Compare to true top singular values
svd_true <- svd(X_true)
cat("True top singular values:\n")
print(round(svd_true$d[1:6], 3))

# Optional: inspect singular-value traces
if (!is.null(res$svals_mat)) {
  svals_mean <- colMeans(res$svals_mat, na.rm = TRUE)
  cat("Mean of saved singular values (across saved iterations):\n")
  print(round(svals_mean, 3))
}


# Bayesian sparse matrix factorization (Gibbs) — spike-and-slab approx on H
bayes_svd_sparse_gibbs <- function(Y, K = 10, niter = 2000, burnin = 1000, thin = 5,
                                   a0 = 1e-3, b0 = 1e-3,
                                   a_pi = 1, b_pi = 1,
                                   spike_prec = 1e6,
                                   c0 = 1e-3, d0 = 1e-3,
                                   verbose = TRUE, seed = 123) {
  set.seed(seed)
  n <- nrow(Y); p <- ncol(Y)
  # initialize
  W <- matrix(rnorm(n * K, sd = 0.1), n, K)
  H <- matrix(rnorm(p * K, sd = 0.1), p, K)
  alpha <- rep(1, K)
  pi_k <- rep(0.5, K)
  Z <- matrix(rbinom(p * K, 1, 0.5), p, K)  # inclusion indicators
  tau <- 1.0
  
  # storage
  nsave <- floor((niter - burnin) / thin)
  if (nsave < 1) nsave <- 0
  Xsum <- matrix(0, n, p)
  alpha_sum <- rep(0, K)
  pi_sum <- rep(0, K)
  tau_sum <- 0
  svals_mat <- matrix(NA, nrow = nsave, ncol = K)
  save_i <- 0
  
  for (iter in 1:niter) {
    # --- sample W (rows) given H, alpha, tau
    S_H <- crossprod(H)                      # K x K
    Precision_w <- tau * S_H + diag(alpha)  # K x K
    Rw <- chol(Precision_w)                  # upper triangular
    Sigma_w <- chol2inv(Rw)
    # compute means for all rows: T = tau * Sigma_w %*% t(Y %*% H)
    Tmat <- tau * Sigma_w %*% t(Y %*% H)
    Zmat <- matrix(rnorm(K * n), nrow = K)
    RinvZ <- backsolve(Rw, Zmat)             # K x n
    Wt <- Tmat + RinvZ
    W <- t(Wt)                               # n x K
    
    # --- sample H (rows) given W, alpha, tau, Z (sparse prior)
    S_W <- crossprod(W)                      # K x K
    # sample each row h_j (length K) separately because prior precision diagonal depends on Z[j,]
    for (j in 1:p) {
      prior_prec_vec <- ifelse(Z[j, ] == 1, alpha, spike_prec) # length K
      Precision_h <- tau * S_W + diag(prior_prec_vec)          # K x K
      Rh <- chol(Precision_h)
      Sigma_h <- chol2inv(Rh)
      Tmat_h <- tau * Sigma_h %*% (t(W) %*% Y[, j])            # K x 1
      z_h <- rnorm(K)
      Rinv_z <- backsolve(Rh, z_h)
      hj <- as.vector(Tmat_h + Rinv_z)
      H[j, ] <- hj
    }
    
    # --- sample Z (inclusion indicators) given H and pi_k
    # P(z_jk = 1 | h_jk) ∝ pi_k * N(h_jk | 0, var_slab=1/alpha_k)
    #                       (1-pi_k) * N(h_jk | 0, var_spike = 1/spike_prec)
    for (k in 1:K) {
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
    
    # --- sample pi_k (Beta prior)
    for (k in 1:K) {
      sumz <- sum(Z[, k])
      pi_k[k] <- rbeta(1, a_pi + sumz, b_pi + p - sumz)
    }
    
    # --- sample alpha_k (ARD) using only W and slab entries of H
    for (k in 1:K) {
      # slab entries are those with Z==1
      Jk <- sum(Z[, k])
      shape <- a0 + n / 2 + Jk / 2
      # sum squares: W column + only H entries where Z==1
      ss_w <- sum(W[, k]^2)
      ss_h_slab <- sum(H[Z[, k] == 1, k]^2)
      rate <- b0 + 0.5 * (ss_w + ss_h_slab)
      alpha[k] <- rgamma(1, shape = shape, rate = rate)
    }
    
    # --- sample tau (noise precision)
    E <- Y - W %*% t(H)
    ressq <- sum(E^2)
    shape_t <- c0 + (n * p) / 2
    rate_t <- d0 + 0.5 * ressq
    tau <- rgamma(1, shape = shape_t, rate = rate_t)
    
    # --- save posterior summaries
    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      save_i <- save_i + 1
      Xcur <- W %*% t(H)
      Xsum <- Xsum + Xcur
      alpha_sum <- alpha_sum + alpha
      pi_sum <- pi_sum + pi_k
      tau_sum <- tau_sum + tau
      svals <- svd(Xcur, nu = 0, nv = 0)$d
      m <- length(svals)
      if (m < K) svals <- c(svals, rep(0, K - m))
      svals_mat[save_i, ] <- svals[1:K]
    }
    
    if (verbose && (iter %% max(1, floor(niter/10)) == 0))
      cat(sprintf("Iter %d / %d  (alpha[1]=%.3g, pi[1]=%.3g, tau=%.3g)\n",
                  iter, niter, alpha[1], pi_k[1], tau))
  }
  
  if (nsave > 0) {
    Xhat <- Xsum / nsave
    alpha_mean <- alpha_sum / nsave
    pi_mean <- pi_sum / nsave
    tau_mean <- tau_sum / nsave
  } else {
    Xhat <- W %*% t(H)
    alpha_mean <- alpha
    pi_mean <- pi_k
    tau_mean <- tau
  }
  
  list(Xhat = Xhat,
       alpha_mean = alpha_mean,
       pi_mean = pi_mean,
       tau_mean = tau_mean,
       svals_mat = if (nsave > 0) svals_mat else NULL,
       last_W = W,
       last_H = H,
       last_Z = Z)
}

# ---------------------------
# Exemple synthétique (H sparse)
set.seed(2026)
n <- 60; p <- 80; r_true <- 3
# true factors
U_true <- matrix(rnorm(n * r_true), n, r_true)
V_true <- matrix(0, p, r_true)
# build V_true with sparse columns
for (k in 1:r_true) {
  idx <- sample(1:p, size = floor(p * 0.15)) # 15% nonzeros per component
  V_true[idx, k] <- rnorm(length(idx))
}
D_true <- diag(c(5, 2.5, 1.2))
X_true <- U_true %*% D_true %*% t(V_true)
sigma <- 0.5
Y <- X_true + matrix(rnorm(n * p, sd = sigma), n, p)

# Run sampler (short demo)
res <- bayes_svd_sparse_gibbs(Y, K = 6, niter = 10000, burnin = 1000, thin = 5,
                              a0 = 1e-3, b0 = 1e-3,
                              a_pi = 1, b_pi = 1,
                              spike_prec = 1e6,
                              c0 = 1e-3, d0 = 1e-3,
                              verbose = TRUE, seed = 42)


# Résultats sommaires
cat("Posterior mean tau:", res$tau_mean, "\n")
cat("Posterior mean alpha:", round(res$alpha_mean, 3), "\n")
cat("Posterior mean pi (prob incl.):", round(res$pi_mean, 3), "\n")

svd_hat <- svd(res$Xhat)
cat("Top singular values (posterior mean X):\n"); print(round(svd_hat$d[1:6], 3))
svd(X_true)$d[1:3]

D_true



library(rstan)
model=stan_model(file = "OneDrive/JDS2026/SVD_sparse.stan")
Y


sigma <- 0.5

Y <- X_true + matrix(rnorm(n * p, sd = sigma), n, p)

res=sampling(model,data=list(r=3,n=dim(Y)[1],p=dim(Y)[2],Y=Y))

fit=extract(res)

mean(fit$sigma)

svd(X_true)$d[1:3]
for(i in 1)
