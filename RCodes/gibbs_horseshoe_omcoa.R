## =====================================================================
##  OMCOA bayésienne parcimonieuse : prior horseshoe sur les Sigma_k
## ---------------------------------------------------------------------
##  Modèle :  X_k = U Sigma_k V_k^T + E_k,   E_k ~ N(0, 1/tau_k)
##            sigma_k^(r) | eta, zeta ~ N+(0, eta_kr^2 zeta^2)
##            eta_kr ~ C+(0,1)  (local),   zeta ~ C+(0, zeta0)  (global)
##  Représentation par variables auxiliaires (Makalic & Schmidt, 2016) :
##            eta_kr^2 | nu_kr ~ IG(1/2, 1/nu_kr),  nu_kr ~ IG(1/2, 1)
##            zeta^2   | xi    ~ IG(1/2, 1/xi),     xi    ~ IG(1/2, 1/zeta0^2)
##  -> toutes les lois conditionnelles sont standard (Gibbs exact).
##
##  Pas de zéros exacts : on utilise le facteur de shrinkage
##      kappa_kr = 1 / (1 + tau_k eta_kr^2 zeta^2)  in (0,1)
##  (kappa ~ 1 : poids écrasé vers 0 ; kappa ~ 0 : poids laissé libre)
##  et on déclare sigma_k^(r) actif si E[1 - kappa_kr | X] > 1/2
##  (Carvalho, Polson & Scott, 2010).
##
##  Réutilise les outils de gibbs_sparse_omcoa.R (même dossier).
## =====================================================================

source("gibbs_sparse_omcoa.R")   # rvmf_vector, rmf_matrix_sweep, rtnorm_pos, omcoa_als, perms

rinvgamma <- function(n, shape, rate) 1 / rgamma(n, shape, rate)

gibbs_horseshoe_omcoa <- function(Xs, R, n_iter = 6000, burn = 2000, init = NULL,
                                  a0 = 1e-3, b0 = 1e-3, zeta0 = 1,
                                  scale_by_noise = TRUE, verbose = TRUE) {
  ## scale_by_noise = TRUE : sigma_k^(r) ~ N+(0, eta^2 zeta^2 / tau_k)  (Wang et al., 2025)
  ## scale_by_noise = FALSE: sigma_k^(r) ~ N+(0, eta^2 zeta^2)
  K <- length(Xs); N <- nrow(Xs[[1]]); P <- sapply(Xs, ncol)
  stopifnot(R <= min(P))
  Q <- qr.Q(qr(scale(diag(N), scale = FALSE)))[, 1:(N - 1)]

  if (is.null(init)) init <- omcoa_als(Xs, R)
  s  <- svd(crossprod(Q, init$U)); Ut <- s$u %*% t(s$v); U <- Q %*% Ut
  V  <- lapply(init$V, procrustes)
  Sm <- do.call(rbind, lapply(init$S, function(x) pmax(abs(x), 1e-3)))   # K x R
  eta2 <- matrix(1, K, R); nu <- matrix(1, K, R); zeta2 <- 1; xi <- 1
  rss <- function(k) sum((Xs[[k]] - U %*% diag(Sm[k, ], R) %*% t(V[[k]]))^2)
  tau <- sapply(seq_len(K), function(k) ((N - 1) * P[k]) / rss(k))

  n_keep <- n_iter - burn
  out <- list(S = array(NA, c(n_keep, K, R)), kappa = array(NA, c(n_keep, K, R)),
              U = array(NA, c(n_keep, N, R)), zeta2 = numeric(n_keep),
              tau = matrix(NA, n_keep, K), trace_S = array(NA, c(n_iter, K, R)))

  for (it in seq_len(n_iter)) {
    ## 1) U | .
    FU <- Reduce(`+`, lapply(seq_len(K), function(k) tau[k] * Xs[[k]] %*% V[[k]] %*% diag(Sm[k, ], R)))
    Ut <- rmf_matrix_sweep(crossprod(Q, FU), Ut); U <- Q %*% Ut

    for (k in seq_len(K)) {
      ## 2) V_k | .
      V[[k]] <- rmf_matrix_sweep(tau[k] * crossprod(Xs[[k]], U) %*% diag(Sm[k, ], R), V[[k]])
      ## 3) sigma_k | . ~ N+(m, 1/prec),  prec = tau_k + 1/(eta^2 zeta^2)
      d    <- diag(t(U) %*% Xs[[k]] %*% V[[k]])
      w    <- if (scale_by_noise) tau[k] else 1          # facteur d'échelle du prior
      prec <- tau[k] + w / (eta2[k, ] * zeta2)
      Sm[k, ] <- rtnorm_pos(tau[k] * d / prec, 1 / sqrt(prec))
      ## 4) tau_k | .  (le prior de sigma_k contribue si scale_by_noise)
      if (scale_by_noise)
        tau[k] <- rgamma(1, a0 + ((N - 1) * P[k] + R) / 2,
                         b0 + rss(k) / 2 + sum(Sm[k, ]^2 / (eta2[k, ] * zeta2)) / 2)
      else
        tau[k] <- rgamma(1, a0 + (N - 1) * P[k] / 2, b0 + rss(k) / 2)
    }
    ## 5) échelles locales et globale (horseshoe)
    W     <- if (scale_by_noise) matrix(tau, K, R) else matrix(1, K, R)
    eta2  <- matrix(rinvgamma(K * R, 1, 1 / nu + W * Sm^2 / (2 * zeta2)), K, R)
    nu    <- matrix(rinvgamma(K * R, 1, 1 + 1 / eta2), K, R)
    zeta2 <- rinvgamma(1, (K * R + 1) / 2, 1 / xi + sum(W * Sm^2 / eta2) / 2)
    xi    <- rinvgamma(1, 1, 1 / zeta0^2 + 1 / zeta2)

    out$trace_S[it, , ] <- Sm
    if (it > burn) {
      i <- it - burn
      out$S[i, , ] <- Sm
      out$kappa[i, , ] <- if (scale_by_noise) 1 / (1 + eta2 * zeta2) else 1 / (1 + tau * eta2 * zeta2)
      out$U[i, , ] <- U; out$zeta2[i] <- zeta2; out$tau[i, ] <- tau
    }
    if (verbose && it %% 1000 == 0) cat("itération", it, "/", n_iter, "\n")
  }
  out
}

align_hs <- function(out, U_ref = out$U[1, , ]) {
  R <- ncol(U_ref); PP <- perms(R)
  for (i in seq_len(dim(out$S)[1])) {
    C <- crossprod(U_ref, out$U[i, , ])
    p <- PP[which.max(apply(PP, 1, function(q) sum(abs(C[cbind(1:R, q)])))), ]
    out$S[i, , ] <- out$S[i, , p]; out$kappa[i, , ] <- out$kappa[i, , p]
    out$U[i, , ] <- sweep(out$U[i, , p], 2, sign(C[cbind(1:R, p)]), "*")
  }
  out
}

## Résumé : pour chaque (k, r)
##  - 1 - E[kappa]  : règle classique du horseshoe (Carvalho et al., 2010) ;
##  - P(share > c)  : probabilité a posteriori que la dimension r explique plus
##                    d'une part c de la variance du bloc k (share = sigma^2 / ||X_k||^2).
hs_summary <- function(out, Xs, c = 0.01) {
  K <- dim(out$S)[2]; R <- dim(out$S)[3]
  do.call(rbind, lapply(1:K, function(k) do.call(rbind, lapply(1:R, function(r) {
    x <- out$S[, k, r]; share <- x^2 / sum(Xs[[k]]^2)
    data.frame(block = k, dim = r,
               one_minus_kappa = 1 - mean(out$kappa[, k, r]),
               p_share = mean(share > c),
               mean = mean(x), lower = unname(quantile(x, .025)), upper = unname(quantile(x, .975)),
               share_mean = mean(share))
  }))))
}

## =====================================================================
##  Démonstration : même simulation que gibbs_sparse_omcoa.R
## =====================================================================
if (sys.nframe() == 0) {
  set.seed(2026)
  N <- 100; P <- c(10, 5, 15); K <- 3; R_true <- 3; R <- 4
  sig_true <- rbind(c(.7, .5, 0), c(.6, .4, 0), c(.5, 0, .6))
  orth <- function(n, h) qr.Q(qr(matrix(rnorm(n * h), n, h)))
  U0 <- qr.Q(qr(scale(matrix(rnorm(N * R_true), N, R_true), scale = FALSE)))
  Xs <- lapply(1:K, function(k) {
    S <- U0 %*% diag(sig_true[k, ]) %*% t(orth(P[k], R_true))
    E <- scale(matrix(rnorm(N * P[k]), N, P[k]), scale = FALSE)
    S + E * sqrt(0.10 / 0.90 * sum(S^2) / sum(E^2))
  })

  t0 <- Sys.time()
  out <- gibbs_horseshoe_omcoa(Xs, R)
  cat("Durée :", format(Sys.time() - t0), "\n")
  U_ref <- cbind(U0, qr.Q(qr(cbind(U0, scale(rnorm(N), scale = FALSE))))[, R_true + 1])
  out <- align_hs(out, U_ref)

  res <- hs_summary(out, Xs, c = 0.01)
  res$truth <- cbind(sig_true, 0)[cbind(res$block, res$dim)]
  print(res, digits = 3)
  cat("\nRègle 1 - E[kappa] > 1/2 :\n")
  print(round(matrix(res$one_minus_kappa, K, R, byrow = TRUE), 3))
  cat("Structure retrouvée :", all((res$truth > 0) == (res$one_minus_kappa > 0.5)), "\n")
  cat("\nRègle P(share > 1 %) > 1/2 :\n")
  print(round(matrix(res$p_share, K, R, byrow = TRUE), 3))
  cat("Structure retrouvée :", all((res$truth > 0) == (res$p_share > 0.5)), "\n")
}
