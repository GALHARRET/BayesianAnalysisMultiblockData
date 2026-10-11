## =====================================================================
##  OMCOA bayésienne : échantillonneur de Gibbs
## ---------------------------------------------------------------------
##  Modèle :  X_k = U Sigma_k V_k^T + E_k,   E_k ~ N(0, 1/tau_k) i.i.d.
##            U  (N x H), U^T U = I, colonnes centrées      (Stiefel)
##            V_k (p_k x H), V_k^T V_k = I                   (Stiefel)
##            Sigma_k = diag(sigma_k^(1..H)),  sigma_k^(h) > 0
##
##  Priors :  U, V_k uniformes sur la variété de Stiefel
##            sigma_k^(h) ~ N+(0, 1/lambda)
##            tau_k       ~ Gamma(a0, b0)
##
##  Lois conditionnelles complètes (cf. Hoff, 2007, 2009) :
##   U       | . ~ von Mises-Fisher matriciel,  F_U  = sum_k tau_k X_k V_k Sigma_k
##   V_k     | . ~ von Mises-Fisher matriciel,  F_Vk = tau_k X_k^T U Sigma_k
##   sigma_k^(h) | . ~ N+( tau_k d_kh / (tau_k + lambda), 1 / (tau_k + lambda) ),
##                      d_kh = u_h^T X_k v_k^(h)
##   tau_k   | . ~ Gamma( a0 + (N-1) p_k / 2,  b0 + ||X_k - U Sigma_k V_k^T||^2 / 2 )
##
##  Lien avec l'OMCOA : les mises à jour de l'OMCOA (Procrustes sur X H^T,
##  Procrustes sur X_k^T U Sigma_k, Sigma_k = diag(U^T X_k V_k)) sont les MODES
##  de ces lois conditionnelles. Le Gibbs est la version stochastique de l'OMCOA :
##  il donne en plus la loi a posteriori (intervalles de crédibilité) des Sigma_k.
##
##  R de base uniquement (pas de dépendance).
## =====================================================================

## ---------------------------------------------------------------- vMF vectoriel
## Tire z ~ vMF(f) sur la sphère S^{m-1} (densité ∝ exp(f^T z)) — algorithme de Wood (1994).
rvmf_vector <- function(f) {
  m <- length(f)
  kappa <- sqrt(sum(f^2))
  if (m == 1) return(sign(f) + (f == 0))
  if (kappa < 1e-12) { z <- rnorm(m); return(z / sqrt(sum(z^2))) }
  mu <- f / kappa
  b  <- (m - 1) / (2 * kappa + sqrt(4 * kappa^2 + (m - 1)^2))   # forme stable
  x0 <- (1 - b) / (1 + b)
  cc <- kappa * x0 + (m - 1) * log(1 - x0^2)
  repeat {
    z <- rbeta(1, (m - 1) / 2, (m - 1) / 2)
    w <- (1 - (1 + b) * z) / (1 - (1 - b) * z)
    if (kappa * w + (m - 1) * log(1 - x0 * w) - cc >= log(runif(1))) break
  }
  v <- rnorm(m); v <- v - sum(v * mu) * mu; v <- v / sqrt(sum(v^2))
  w * mu + sqrt(max(0, 1 - w^2)) * v
}

## ---------------------------------------------------------------- vMF matriciel
## Un balayage de Gibbs colonne par colonne pour W ~ vMF(F) sur la variété de Stiefel
## (Hoff, 2009, §2.2) : w_j | W_{-j} = N_j z,  z ~ vMF(N_j^T f_j),
## N_j base orthonormée de l'orthogonal des autres colonnes.
rmf_matrix_sweep <- function(F, W) {
  m <- nrow(W); h <- ncol(W)
  for (j in sample(h)) {
    if (h > 1) {
      Q  <- qr.Q(qr(W[, -j, drop = FALSE]), complete = TRUE)
      Nj <- Q[, h:m, drop = FALSE]                 # orthogonal de span(W[,-j])
    } else Nj <- diag(m)
    W[, j] <- Nj %*% rvmf_vector(crossprod(Nj, F[, j]))
  }
  W
}

## ---------------------------------------------------------------- loi normale tronquée > 0
rtnorm_pos <- function(n, mean, sd) {
  p0 <- pnorm(0, mean, sd)
  u  <- runif(n, p0, 1)
  pmax(qnorm(pmin(u, 1 - 1e-16), mean, sd), 1e-12)
}

## ---------------------------------------------------------------- échantillonneur
## Xs    : liste de K blocs centrés (N x p_k)
## H     : dimension de l'espace latent
## init  : liste(U, V, S) (par ex. solution MCOA ou OMCOA) ; NULL = aléatoire
gibbs_omcoa <- function(Xs, H, n_iter = 3000, burn = 1000, thin = 1, init = NULL,
                        a0 = 1e-3, b0 = 1e-3, lambda = 1e-2, verbose = TRUE) {
  K <- length(Xs); N <- nrow(Xs[[1]]); P <- sapply(Xs, ncol)
  # base orthonormée de 1^perp : U = Q %*% Ut, avec Ut ((N-1) x H) sur Stiefel
  Q <- qr.Q(qr(scale(diag(N), scale = FALSE)))[, 1:(N - 1)]

  if (is.null(init)) {
    Ut <- qr.Q(qr(matrix(rnorm((N - 1) * H), N - 1, H)))
    V  <- lapply(P, function(p) qr.Q(qr(matrix(rnorm(p * H), p, H))))
    U  <- Q %*% Ut
    S  <- lapply(seq_len(K), function(k) pmax(abs(diag(t(U) %*% Xs[[k]] %*% V[[k]])), 1e-3))
  } else {
    U  <- init$U
    s  <- svd(crossprod(Q, U)); Ut <- s$u %*% t(s$v); U <- Q %*% Ut   # projection sur Stiefel
    V  <- lapply(init$V, function(v) { s <- svd(v); s$u %*% t(s$v) })
    S  <- lapply(init$S, function(s) pmax(abs(s), 1e-3))
  }
  rss <- function(k) sum((Xs[[k]] - U %*% diag(S[[k]], H) %*% t(V[[k]]))^2)
  tau <- sapply(seq_len(K), function(k) ((N - 1) * P[k]) / rss(k))

  n_keep <- floor((n_iter - burn) / thin)
  draws <- list(U = array(NA, c(n_keep, N, H)),
                V = lapply(P, function(p) array(NA, c(n_keep, p, H))),
                S = array(NA, c(n_keep, K, H)),
                tau = matrix(NA, n_keep, K))
  trace_S <- array(NA, c(n_iter, K, H))      # trace complète (burn-in inclus)
  i_keep <- 0

  for (it in seq_len(n_iter)) {
    ## 1) U | .  ~ vMF(sum_k tau_k X_k V_k Sigma_k), restreint à 1^perp
    FU <- Reduce(`+`, lapply(seq_len(K), function(k) tau[k] * Xs[[k]] %*% V[[k]] %*% diag(S[[k]], H)))
    Ut <- rmf_matrix_sweep(crossprod(Q, FU), Ut)
    U  <- Q %*% Ut
    for (k in seq_len(K)) {
      ## 2) V_k | . ~ vMF(tau_k X_k^T U Sigma_k)
      V[[k]] <- rmf_matrix_sweep(tau[k] * crossprod(Xs[[k]], U) %*% diag(S[[k]], H), V[[k]])
      ## 3) sigma_k | . ~ N+ (conjuguée)
      d <- diag(t(U) %*% Xs[[k]] %*% V[[k]])
      prec <- tau[k] + lambda
      S[[k]] <- rtnorm_pos(H, tau[k] * d / prec, 1 / sqrt(prec))
      ## 4) tau_k | . ~ Gamma (conjuguée)
      tau[k] <- rgamma(1, a0 + (N - 1) * P[k] / 2, b0 + rss(k) / 2)
    }
    trace_S[it, , ] <- do.call(rbind, S)
    if (it > burn && (it - burn) %% thin == 0) {
      i_keep <- i_keep + 1
      draws$U[i_keep, , ] <- U
      for (k in seq_len(K)) draws$V[[k]][i_keep, , ] <- V[[k]]
      draws$S[i_keep, , ] <- do.call(rbind, S)
      draws$tau[i_keep, ] <- tau
    }
    if (verbose && it %% 500 == 0) cat("itération", it, "/", n_iter, "\n")
  }
  draws$trace_S <- trace_S
  draws
}

## ---------------------------------------------------------------- post-traitement
## La loi a posteriori est invariante par permutation des H dimensions (et par
## changement de signe conjoint de u_h et des v_k^(h)). On aligne chaque tirage
## sur une référence U_ref (par ex. la solution OMCOA) : permutation + signe.
perms <- function(n) {
  if (n == 1) return(matrix(1))
  p <- perms(n - 1)
  do.call(rbind, lapply(seq_len(n), function(i) cbind(i, ifelse(p >= i, p + 1, p))))
}

align_draws <- function(draws, U_ref) {
  H <- ncol(U_ref); PP <- perms(H); n <- dim(draws$S)[1]; K <- dim(draws$S)[2]
  for (i in seq_len(n)) {
    C <- crossprod(U_ref, draws$U[i, , ])
    sc <- apply(PP, 1, function(p) sum(abs(C[cbind(seq_len(H), p)])))
    p  <- PP[which.max(sc), ]
    sg <- sign(C[cbind(seq_len(H), p)])
    draws$U[i, , ] <- sweep(draws$U[i, , p, drop = FALSE][1, , ], 2, sg, "*")
    for (k in seq_len(K))
      draws$V[[k]][i, , ] <- sweep(draws$V[[k]][i, , p, drop = FALSE][1, , ], 2, sg, "*")
    draws$S[i, , ] <- draws$S[i, , p]
  }
  draws
}

## Résumé a posteriori de Sigma_k : moyenne et intervalle de crédibilité
summary_sigma <- function(draws, level = 0.95) {
  a <- (1 - level) / 2
  K <- dim(draws$S)[2]; H <- dim(draws$S)[3]
  do.call(rbind, lapply(seq_len(K), function(k) do.call(rbind, lapply(seq_len(H), function(h) {
    x <- draws$S[, k, h]
    data.frame(block = k, dim = h, mean = mean(x),
               lower = unname(quantile(x, a)), upper = unname(quantile(x, 1 - a)))
  }))))
}

## =====================================================================
##  Démonstration : profils permutés (t = 1), même simulation que
##  omcoa_vs_mcoa.R ; le Gibbs est volontairement initialisé sur la MCOA.
## =====================================================================
if (sys.nframe() == 0) {
  set.seed(2026)
  N <- 100; P <- c(10, 5, 15); H <- 3; K <- 3
  sig_true <- rbind(c(.7, .5, .3), c(.5, .3, .7), c(.3, .7, .5))
  orth <- function(n, h) qr.Q(qr(matrix(rnorm(n * h), n, h)))
  U0 <- qr.Q(qr(scale(matrix(rnorm(N * H), N, H), scale = FALSE)))
  Es <- vector("list", K)
  Xs <- lapply(seq_len(K), function(k) {
    S <- U0 %*% diag(sig_true[k, ]) %*% t(orth(P[k], H))
    E <- scale(matrix(rnorm(N * P[k]), N, P[k]), scale = FALSE)
    Es[[k]] <<- E * sqrt(0.10 / 0.90 * sum(S^2) / sum(E^2))
    S + Es[[k]]
  })

  ## MCOA (initialisation)
  mcoa <- function(Xs) {
    Xt <- Xs; U <- matrix(0, N, H)
    V <- lapply(P, function(p) matrix(0, p, H)); S <- lapply(1:K, function(k) numeric(H))
    for (h in 1:H) {
      u <- svd(do.call(cbind, Xt), nu = 1, nv = 0)$u[, 1]; U[, h] <- u
      for (k in 1:K) {
        w <- crossprod(Xt[[k]], u); S[[k]][h] <- sqrt(sum(w^2)); V[[k]][, h] <- w / S[[k]][h]
        Xt[[k]] <- Xt[[k]] - Xt[[k]] %*% V[[k]][, h] %*% t(V[[k]][, h])
      }
    }
    list(U = U, V = V, S = S)
  }
  m <- mcoa(Xs)

  t0 <- Sys.time()
  dr <- gibbs_omcoa(Xs, H, n_iter = 6000, burn = 2000, init = m)
  cat("Durée :", format(Sys.time() - t0), "\n")

  dr <- align_draws(dr, U0)          # alignement sur le vrai U (ou sur l'OMCOA en pratique)
  res <- summary_sigma(dr)
  res$truth <- sig_true[cbind(res$block, res$dim)]
  res$covered <- res$lower <= res$truth & res$truth <= res$upper
  print(res, digits = 3)
  cat("Couverture des IC 95 % :", mean(res$covered), "\n")
  cat("Bruit estimé (1/tau_k) vs vrai :\n")
  print(rbind(estime = colMeans(1 / dr$tau),
              vrai = sapply(1:K, function(k) sum(Es[[k]]^2) / ((N - 1) * P[k]))))

  ## Trace des sigma_k^(h) : départ MCOA -> convergence vers la solution OMCOA
  library(ggplot2)
  tr <- dr$trace_S
  d_tr <- do.call(rbind, lapply(1:K, function(k) do.call(rbind, lapply(1:H, function(h)
    data.frame(iter = seq_len(dim(tr)[1]), value = tr[, k, h],
               block = paste("Bloc", k), comp = paste0("comp", h))))))
  g <- ggplot(d_tr, aes(iter, value, colour = comp)) +
    geom_line(linewidth = 0.4) +
    geom_vline(xintercept = 2000, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ block) +
    scale_colour_manual(values = c("#2a78d6", "#eb6834", "#1baf7a"), name = NULL) +
    labs(x = "Itération (pointillés : fin du burn-in)", y = expression(sigma[k]^{(h)}),
         title = "Gibbs initialisé sur la MCOA, profils permutés (t = 1)",
         subtitle = "Composantes dans l'ordre du tirage (avant alignement)") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  ggsave("gibbs_trace.png", g, width = 10, height = 4, dpi = 200)
}
