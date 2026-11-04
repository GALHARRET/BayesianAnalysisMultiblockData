## =============================================================================
## Rank selection for a (concatenated) data matrix X:
##   1. Bayesian SVD with model averaging over the rank  (Hoff, 2007, JASA)
##   2. Gavish & Donoho (2014) optimal hard threshold     (IEEE Trans. Inf. Th.)
##
## Both assume  X = signal of low rank + i.i.d. Gaussian noise  on ALL entries.
## For the multiblock simulation: apply them to the concatenated, column-centred
## X = [X_1 | ... | X_K] WITHOUT block standardisation (a common noise level s
## for all blocks satisfies the i.i.d. assumption; block standardisation breaks it).
## =============================================================================


## -----------------------------------------------------------------------------
## Utilities
## -----------------------------------------------------------------------------

## log(sum(exp(x))) computed stably
.lse <- function(x) { m <- max(x); if (!is.finite(m)) return(m); m + log(sum(exp(x - m))) }

## Orthonormal basis of the orthogonal complement of the columns of A (m x k)
.null_basis <- function(A, m) {
  if (is.null(A) || ncol(A) == 0) return(diag(m))
  qr.Q(qr(A), complete = TRUE)[, -(seq_len(ncol(A))), drop = FALSE]
}

## Orthonormal basis of the orthogonal complement of 1_N (used to remove centring)
centring_basis <- function(N) .null_basis(matrix(1 / sqrt(N), N, 1), N)

## One draw from the vector von Mises-Fisher distribution vMF(b) on the unit
## sphere of R^p (density proportional to exp(b'x)); Wood (1994).
rvmf <- function(b) {
  p <- length(b); kappa <- sqrt(sum(b^2))
  if (p == 1) return(if (runif(1) < plogis(2 * b)) 1 else -1)
  if (kappa < 1e-12) { z <- rnorm(p); return(z / sqrt(sum(z^2))) }
  mu <- b / kappa
  b0 <- (p - 1) / (2 * kappa + sqrt(4 * kappa^2 + (p - 1)^2))
  x0 <- (1 - b0) / (1 + b0)
  cc <- kappa * x0 + (p - 1) * log(1 - x0^2)
  repeat {
    z <- rbeta(1, (p - 1) / 2, (p - 1) / 2)
    w <- (1 - (1 + b0) * z) / (1 - (1 - b0) * z)
    if (kappa * w + (p - 1) * log(1 - x0 * w) - cc >= log(runif(1))) break
  }
  v <- rnorm(p); v <- v - sum(v * mu) * mu; v <- v / sqrt(sum(v^2))
  w * mu + sqrt(1 - w^2) * v
}

## log E[X^(2l)], l = 0..L, for X ~ N(a, s2)   (recursion on raw moments)
.log_even_moments <- function(a, s2, L) {
  a <- abs(a)                                     # even moments are symmetric in a
  if (a < 1e-12 * sqrt(s2)) {                     # centred case: s^(2l) (2l-1)!!
    l <- 0:L
    return(l * log(s2) + lgamma(2 * l + 1) - l * log(2) - lgamma(l + 1))
  }
  n <- 2 * L
  lm <- numeric(n + 1)                            # lm[k+1] = log E[X^k]
  r <- a                                          # r_k = E[X^k] / E[X^(k-1)]
  lm[2] <- log(a)
  for (k in 2:n) {
    r <- a + (k - 1) * s2 / r
    lm[k + 1] <- lm[k] + log(r)
  }
  lm[seq(1, n + 1, by = 2)]
}

## M_l = E[(lambda' q)^l], l = 0..L, q ~ Dirichlet(1/2, ..., 1/2) of dimension
## mdim, lambda = non-zero eigenvalues (normalised to sum to 1 for stability);
## recursion from prod_i (1 - t lambda_i)^(-1/2) = sum_l (A)_l t^l M_l / l!,
## A = mdim / 2  (Hoff, 2007, appendix).
.dirichlet_moments <- function(lambda, mdim, L) {
  A <- mdim / 2
  pk <- sapply(1:(L + 1), function(k) 0.5 * sum(lambda^k))  # p_k, k = 1..L+1
  M <- numeric(L + 1); M[1] <- 1
  lpoch <- function(x, l) lgamma(x + l) - lgamma(x)
  for (l in 0:(L - 1)) {
    k <- 0:l
    lr <- lgamma(l + 1) + lpoch(A, l - k) - lpoch(A, l + 1) - lgamma(l - k + 1)
    M[l + 2] <- sum(exp(lr) * M[l - k + 1] * pk[k + 1])
  }
  M
}

## -----------------------------------------------------------------------------
## Bayes factor of d_j != 0 against d_j = 0, with (u_j, v_j) uniform on the
## spheres of the null spaces and d_j ~ N(mu, 1/psi) integrated out (Hoff, 2007):
##   BF = sqrt(psi/omega) exp(-mu^2 psi phi / (2 omega))
##        * sum_l (phi^2/4)^l  E_u[(u'EE'u)^l]  E[D^(2l)] / ((n~/2)_l l!),
##   omega = phi + psi,  D ~ N(mu psi / omega, 1 / omega).
## E_u[(u'EE'u)^l] = ||E||^(2l) E[(lambda~'q)^l], q ~ Dirichlet(1/2) ("exact");
## it lies between sigma1^(2l) (1/2)_l / (m~/2)_l ("lower", rank-one part of E)
## and sigma1^(2l) ("upper"). All terms are positive, so the bounds bracket BF.
## -----------------------------------------------------------------------------

## common part of the log terms (everything except E_u[(u'EE'u)^l])
.log_common <- function(nt, phi, mu, psi, L, lmom) {
  omega <- phi + psi; l <- 0:L
  0.5 * log(psi / omega) - 0.5 * mu^2 * psi * phi / omega +
    l * log(phi^2 / 4) + lmom - (lgamma(nt / 2 + l) - lgamma(nt / 2)) - lgamma(l + 1)
}
.log_Eu <- function(sv, mt, L, type) {
  l <- 0:L; lam <- sv^2
  switch(type,
    upper = l * log(lam[1]),
    lower = l * log(lam[1]) + lgamma(0.5 + l) - lgamma(0.5) -
            (lgamma(mt / 2 + l) - lgamma(mt / 2)),
    ## normalising by the largest eigenvalue keeps M_l >= (1/2)_l / (m~/2)_l,
    ## which avoids underflow for large l (normalising by the sum does not)
    exact = { M <- .dirichlet_moments(lam / lam[1], mt, L)
              if (any(M <= 0)) stop("underflow in the Dirichlet moments")
              l * log(lam[1]) + log(M) })
}

## cheap lower bound on log BF (no recursion): E[D^(2l)] >= max(a^(2l), s^(2l)(2l-1)!!)
.cheap_lower_logBF <- function(sv, mt, nt, phi, mu, psi, L) {
  omega <- phi + psi; a <- abs(mu * psi / omega); l <- 0:L
  lmom <- pmax(2 * l * log(max(a, 1e-300)),
               -l * log(omega) + lgamma(2 * l + 1) - l * log(2) - lgamma(l + 1))
  max(.log_common(nt, phi, mu, psi, L, lmom) + .log_Eu(sv, mt, L, "lower"))
}

## number of terms of the series (peak near l ~ phi sigma1^2 + phi sigma1 |a|)
.guess_L <- function(sv, phi, mu, psi) {
  a <- abs(mu * psi / (phi + psi))
  ceiling(3 * (phi * sv[1]^2 + phi * sv[1] * a) + 100)
}

## draw d from the density proportional to |d|^(2l) exp(-omega (d - a)^2 / 2)
.draw_d_given_l <- function(l, a, omega) {
  if (l == 0) return(rnorm(1, a, 1 / sqrt(omega)))
  r <- sqrt(a^2 + 8 * l / omega)
  lo <- min((a - r) / 2, 0) - 10 / sqrt(omega); hi <- max((a + r) / 2, 0) + 10 / sqrt(omega)
  g <- seq(lo, hi, length.out = 4001); h <- g[2] - g[1]
  ld <- 2 * l * log(pmax(abs(g), 1e-300)) - omega * (g - a)^2 / 2
  i <- sample.int(length(g), 1, prob = exp(ld - max(ld)))
  g[i] + runif(1, -h / 2, h / 2)
}

## Hoff's "empirical Bayes" hyperparameters: averages over the least-squares
## fits of ranks 0..n (n = min(dim(Y))), with nu0 = eta0 = 2 (Hoff 2007, Sect. 5)
hoff_hyper <- function(Y) {
  m <- nrow(Y); n <- ncol(Y); sv <- svd(Y, nu = 0, nv = 0)$d; r <- length(sv)
  ss <- sum(Y^2)
  s2k  <- sapply(0:r, function(k) (ss - sum(sv[seq_len(k)]^2)) / (m * n))
  muk  <- sapply(1:r, function(k) mean(sv[1:k]))
  tauk <- sapply(1:r, function(k) mean((sv[1:k] - mean(sv[1:k]))^2))
  list(nu0 = 2, s20 = mean(s2k), mu0 = mean(muk), v20 = var(muk),
       eta0 = 2, tau20 = mean(tauk))
}

## =============================================================================
## Hoff (2007) Bayesian SVD with unknown rank
## =============================================================================
#' @param Y       m x n data matrix (for column-centred data, pass t(Q) %*% X
#'                with Q = centring_basis(N): this removes the lost degree of freedom)
#' @param Rmax    maximal rank; prior on the rank uniform on {0, ..., Rmax}
#'                (Hoff uses Rmax = min(dim(Y)); a smaller value is much faster)
#' @param n_iter, burn  MCMC length and burn-in
#' @param hyper   list(nu0, s20, mu0, v20, eta0, tau20); default = hoff_hyper(Y)
#' @return list with K (sampled ranks), post_K (posterior of the rank), d (draws
#'         of the singular values), phi (noise precision), hyper
hoff_svd <- function(Y, Rmax = min(dim(Y)), n_iter = 3000, burn = 500,
                     hyper = NULL, n_uv_steps = 10, verbose = FALSE) {
  Y <- as.matrix(Y); m <- nrow(Y); n <- ncol(Y)
  if (Rmax > min(m, n)) stop("Rmax must be <= min(dim(Y))")
  h <- modifyList(hoff_hyper(Y), as.list(hyper))

  ## initialisation as in Hoff (2007): K = 0, hyperparameters at their prior means
  U <- matrix(0, m, Rmax); V <- matrix(0, n, Rmax); d <- numeric(Rmax)
  phi <- 1 / h$s20; mu <- h$mu0; psi <- 1 / h$tau20

  K_tr <- integer(n_iter); d_tr <- matrix(0, n_iter, Rmax); phi_tr <- numeric(n_iter)

  residual_wo <- function(j) {
    o <- setdiff(which(d != 0), j)
    E <- if (length(o)) Y - U[, o, drop = FALSE] %*% (d[o] * t(V[, o, drop = FALSE])) else Y
    list(E = E, o = o)
  }

  for (it in 1:n_iter) {
    ## ---- (A) variable-rank step: is d_j = 0 ?  (u_j, v_j, d_j) integrated out --
    for (j in 1:Rmax) {
      rw <- residual_wo(j); o <- rw$o
      Nu <- .null_basis(U[, o, drop = FALSE], m); Nv <- .null_basis(V[, o, drop = FALSE], n)
      Et <- crossprod(Nu, rw$E %*% Nv); mt <- nrow(Et); nt <- ncol(Et)
      sv <- svd(Et, nu = 0, nv = 0)$d; sv <- sv[sv > 1e-12 * max(sv, 1e-300)]
      K_o <- length(o)
      lprior <- log((K_o + 1) / (Rmax - K_o))          # uniform prior on the rank
      omega <- phi + psi; a <- mu * psi / omega
      L <- .guess_L(sv, phi, mu, psi)

      lt <- NULL
      if (lprior + .cheap_lower_logBF(sv, mt, nt, phi, mu, psi, L) > 35) {
        active <- TRUE                                  # strong signal: no doubt
      } else {
        lmom <- .log_even_moments(a, 1 / omega, L)
        lc <- .log_common(nt, phi, mu, psi, L, lmom)
        up <- .lse(lc + .log_Eu(sv, mt, L, "upper"))
        if (lprior + up < -35) {
          active <- FALSE                               # pure noise: no doubt
        } else {
          lt <- lc + .log_Eu(sv, mt, L, "exact")        # exact Bayes factor
          active <- runif(1) < plogis(lprior + .lse(lt))
        }
      }
      if (!active) { d[j] <- 0; U[, j] <- 0; V[, j] <- 0; next }

      if (!is.null(lt)) {
        ## d_j from its exact marginal: mixture over l of |d|^(2l) N(a, 1/omega)
        l_draw <- sample.int(L + 1, 1, prob = exp(lt - max(lt))) - 1
        d[j] <- .draw_d_given_l(l_draw, a, omega)
      } else {
        ## strong signal: the marginal of d_j is concentrated near sigma1(Et);
        ## start there, (u, v, d) are then refreshed by exact Gibbs steps
        d[j] <- sign(a + 1e-300) * sv[1]
      }

      ## (u, v) | d_j : density exp(phi d u' Et v); start at the leading singular
      ## pair and run a few alternating vMF steps (approximation used by Hoff)
      A <- phi * d[j] * Et
      s1 <- svd(A, nu = 1, nv = 1); u <- s1$u[, 1]; v <- s1$v[, 1]
      for (s in seq_len(n_uv_steps)) { u <- rvmf(A %*% v); v <- rvmf(crossprod(A, u)) }
      if (is.null(lt))
        d[j] <- rnorm(1, (phi * sum(u * (Et %*% v)) + mu * psi) / omega, 1 / sqrt(omega))
      U[, j] <- Nu %*% u; V[, j] <- Nv %*% v
    }

    ## ---- (B) fixed-rank Gibbs updates of the active components --------------
    for (j in which(d != 0)) {
      rw <- residual_wo(j); o <- rw$o; E <- rw$E
      Nu <- .null_basis(U[, o, drop = FALSE], m)
      U[, j] <- Nu %*% rvmf(phi * d[j] * crossprod(Nu, E %*% V[, j]))
      Nv <- .null_basis(V[, o, drop = FALSE], n)
      V[, j] <- Nv %*% rvmf(phi * d[j] * crossprod(Nv, crossprod(E, U[, j])))
      d[j] <- rnorm(1, (phi * sum(U[, j] * (E %*% V[, j])) + mu * psi) / (phi + psi),
                    1 / sqrt(phi + psi))
    }

    ## ---- (C) hyperparameters -------------------------------------------------
    act <- which(d != 0); K <- length(act)
    R <- if (K) Y - U[, act, drop = FALSE] %*% (d[act] * t(V[, act, drop = FALSE])) else Y
    phi <- rgamma(1, (h$nu0 + m * n) / 2, (h$nu0 * h$s20 + sum(R^2)) / 2)
    prec <- psi * K + 1 / h$v20
    mu  <- rnorm(1, (psi * sum(d[act]) + h$mu0 / h$v20) / prec, 1 / sqrt(prec))
    psi <- rgamma(1, (h$eta0 + K) / 2, (h$eta0 * h$tau20 + sum((d[act] - mu)^2)) / 2)

    K_tr[it] <- K; d_tr[it, ] <- d; phi_tr[it] <- phi
    if (verbose && it %% 500 == 0) cat("iter", it, " K =", K, "\n")
  }

  keep <- (burn + 1):n_iter
  list(K = K_tr[keep],
       post_K = table(factor(K_tr[keep], levels = 0:Rmax)) / length(keep),
       d = d_tr[keep, , drop = FALSE], phi = phi_tr[keep], hyper = h)
}


## =============================================================================
## Gavish & Donoho (2014): optimal hard threshold for singular values
## =============================================================================
## Y is m x n; beta = min/max aspect ratio. Singular values above the threshold
## are kept; the estimated rank is their number.
##   known noise sd s :  tau = lambda*(beta) * sqrt(max(m, n)) * s
##   unknown noise    :  tau = omega(beta)   * median singular value of Y
## lambda*(beta) = sqrt(2(beta+1) + 8 beta / (beta + 1 + sqrt(beta^2 + 14 beta + 1)))
## omega(beta)   = lambda*(beta) / sqrt(mu_beta), mu_beta = median of the
##                 Marchenko-Pastur law of ratio beta (computed numerically);
##                 approximation omega ~ 0.56 b^3 - 0.95 b^2 + 1.82 b + 1.43.



gd_threshold <- function(Y, sigma = NULL) {
  gd_lambda_star <- function(beta)
    sqrt(2 * (beta + 1) + 8 * beta / (beta + 1 + sqrt(beta^2 + 14 * beta + 1)))

  mp_median <- function(beta) {
    a <- (1 - sqrt(beta))^2; b <- (1 + sqrt(beta))^2
    ## substitution t = a + (b - a) sin^2(theta) removes the endpoint singularities
    dens_theta <- function(th) {
      t <- a + (b - a) * sin(th)^2
      ((b - a) * sin(th) * cos(th))^2 * 2 / (2 * pi * beta * pmax(t, 1e-300))
    }
    cdf <- function(x) {
      thx <- asin(sqrt(min(max((x - a) / (b - a), 0), 1)))
      integrate(dens_theta, 0, thx, rel.tol = 1e-10)$value
    }
    uniroot(function(x) cdf(x) - 0.5, c(a, b), tol = 1e-12)$root
  }

  gd_omega <- function(beta) gd_lambda_star(beta) / sqrt(mp_median(beta))
  m <- nrow(Y); n <- ncol(Y); beta <- min(m, n) / max(m, n)
  sv <- svd(Y, nu = 0, nv = 0)$d
  tau <- if (is.null(sigma)) gd_omega(beta) * median(sv)
         else gd_lambda_star(beta) * sqrt(max(m, n)) * sigma
  list(threshold = tau, rank = sum(sv > tau), singular_values = sv, beta = beta)
}
