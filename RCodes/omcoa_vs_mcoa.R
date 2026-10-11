## =====================================================================
##  MCOA vs OMCOA : un exemple où les deux méthodes divergent (version R)
## ---------------------------------------------------------------------
##  Modèle : X_k = U Sigma_k V_k^T + E_k,  U (N x H) orthonormé commun,
##           V_k orthonormés, Sigma_k diagonales.
##  Profils permutés : Sigma_1 = diag(a,b,c), Sigma_2 = diag(b,c,a),
##           Sigma_3 = diag(c,a,b)  =>  sum_k Sigma_k^2 = (a²+b²+c²) I.
##  Les valeurs singulières de X sont égales : la MCOA (1re étape = SVD
##  de X) ne peut pas identifier U ; l'OMCOA le peut.
##  t in [0,1] fait passer d'un profil commun (t=0) aux profils permutés (t=1).
## =====================================================================

library(ggplot2)
set.seed(2026)

N <- 100; P <- c(10, 5, 15); H <- 3; K <- length(P)
a <- 0.7; b <- 0.5; c <- 0.3
BASE <- c(a, b, c)
PERM <- list(c(a, b, c), c(b, c, a), c(c, a, b))
NOISE <- 0.10        # part de variance due au bruit dans chaque bloc
R_EXT <- 200         # répétitions pour t = 0 et t = 1
R_MID <- 60          # répétitions pour les t intermédiaires
N_STARTS <- 10       # départs aléatoires de l'OMCOA (en plus de l'init MCOA)

## ---------------------------------------------------------------- outils
orth <- function(n, h, center = FALSE) {
  M <- matrix(rnorm(n * h), n, h)
  if (center) M <- scale(M, scale = FALSE)
  qr.Q(qr(M))
}

procrustes <- function(M) {          # argmax_{Q^T Q = I} trace(Q^T M)
  s <- svd(M)
  s$u %*% t(s$v)
}

crit <- function(Xs, U, S, V) {      # sum_k ||X_k - U Sigma_k V_k^T||^2
  sum(sapply(seq_len(K), function(k)
    sum((Xs[[k]] - U %*% diag(S[[k]], H) %*% t(V[[k]]))^2)))
}

simulate <- function(t, noise = NOISE) {
  U <- orth(N, H, center = TRUE)
  sig <- lapply(PERM, function(s) (1 - t) * BASE + t * s)
  Xs <- vector("list", K)
  for (k in seq_len(K)) {
    V <- orth(P[k], H)
    S <- U %*% diag(sig[[k]]) %*% t(V)
    E <- scale(matrix(rnorm(N * P[k]), N, P[k]), scale = FALSE)
    E <- E * sqrt(noise / (1 - noise) * sum(S^2) / sum(E^2))
    Xs[[k]] <- S + E
  }
  list(Xs = Xs, U = U, sig = sig)
}

## ---------------------------------------------------------------- MCOA
## Approche séquentielle : rang 1 sur X déflaté par les loadings v_k.
mcoa <- function(Xs) {
  Xt <- Xs
  U <- matrix(0, N, H)
  V <- lapply(P, function(p) matrix(0, p, H))
  S <- lapply(seq_len(K), function(k) numeric(H))
  for (h in seq_len(H)) {
    u <- svd(do.call(cbind, Xt), nu = 1, nv = 0)$u[, 1]
    U[, h] <- u
    for (k in seq_len(K)) {
      w <- crossprod(Xt[[k]], u)
      S[[k]][h] <- sqrt(sum(w^2))
      V[[k]][, h] <- w / S[[k]][h]
      v <- V[[k]][, h, drop = FALSE]
      Xt[[k]] <- Xt[[k]] - Xt[[k]] %*% v %*% t(v)     # déflation
    }
  }
  # Sigma_k ajusté (même définition que l'OMCOA) pour comparer les critères
  S_fit <- lapply(seq_len(K), function(k) diag(t(U) %*% Xs[[k]] %*% V[[k]]))
  list(U = U, S = S, V = V, S_fit = S_fit, crit = crit(Xs, U, S_fit, V))
}

## ---------------------------------------------------------------- OMCOA
## Moindres carrés alternés avec U^T U = I et V_k^T V_k = I.
omcoa_once <- function(Xs, U0 = NULL, V0 = NULL, tol = 1e-12, maxit = 5000) {
  X <- do.call(cbind, Xs)
  if (is.null(U0)) {
    U <- orth(N, H); V <- lapply(P, orth, h = H)
  } else {
    U <- procrustes(U0); V <- lapply(V0, procrustes)
  }
  S <- lapply(seq_len(K), function(k) diag(t(U) %*% Xs[[k]] %*% V[[k]]))
  old <- Inf
  for (it in seq_len(maxit)) {
    Hm <- do.call(cbind, lapply(seq_len(K), function(k) diag(S[[k]], H) %*% t(V[[k]])))
    U <- procrustes(X %*% t(Hm))                                        # U
    V <- lapply(seq_len(K), function(k)
      procrustes(t(Xs[[k]]) %*% U %*% diag(S[[k]], H)))                 # V_k : SVD de X_k^T U Sigma_k
    S <- lapply(seq_len(K), function(k) diag(t(U) %*% Xs[[k]] %*% V[[k]]))  # Sigma_k = diag(U^T X_k V_k)
    f <- crit(Xs, U, S, V)
    if (old - f < tol * max(1, f)) break
    old <- f
  }
  for (k in seq_len(K)) {            # convention : Sigma_k >= 0
    sg <- sign(S[[k]]); sg[sg == 0] <- 1
    V[[k]] <- sweep(V[[k]], 2, sg, "*"); S[[k]] <- S[[k]] * sg
  }
  list(U = U, S = S, V = V, crit = f, iter = it)
}

omcoa <- function(Xs, init = NULL, n_starts = N_STARTS) {
  best <- if (is.null(init)) omcoa_once(Xs) else omcoa_once(Xs, init$U, init$V)
  for (s in seq_len(n_starts)) {
    r <- omcoa_once(Xs)
    if (r$crit < best$crit) best <- r
  }
  best
}

## ---------------------------------------------------------------- appariement
perms <- function(n) {
  if (n == 1) return(matrix(1))
  p <- perms(n - 1)
  do.call(rbind, lapply(seq_len(n), function(i) cbind(i, ifelse(p >= i, p + 1, p))))
}
PERMS_H <- perms(H)

## Aligne les colonnes estimées sur les vraies (max de la somme des |cos|)
match_cols <- function(Utrue, Uhat) {
  Uhat <- sweep(Uhat, 2, sqrt(colSums(Uhat^2)), "/")
  C <- abs(crossprod(Utrue, Uhat))
  score <- apply(PERMS_H, 1, function(p) sum(C[cbind(seq_len(H), p)]))
  idx <- PERMS_H[which.max(score), ]
  list(idx = idx, cos = C[cbind(seq_len(H), idx)])
}

run_one <- function(t) {
  sim <- simulate(t)
  m <- mcoa(sim$Xs)
  o <- omcoa(sim$Xs, init = m)
  tot <- sum(sapply(sim$Xs, function(X) sum(X^2)))
  sig_true <- do.call(rbind, sim$sig)                     # K x H
  out <- list()
  for (meth in c("MCOA", "OMCOA")) {
    fit <- if (meth == "MCOA") m else o
    mc <- match_cols(sim$U, fit$U)
    sig_hat <- do.call(rbind, lapply(fit$S, function(s) s[mc$idx]))
    out[[meth]] <- list(cos = mc$cos, sigma = sig_hat,
                        err = mean(abs(sig_hat - sig_true)),
                        fit = 1 - fit$crit / tot)
  }
  out$sig_true <- sig_true
  out
}

res0<-list()
B<-1000
for(b in 1:B){
  res0[[b]]<-run_one(0)
}
save(res0,file="C:/Users/j.galharret/Documents/GitHub/OMCOA/Rcodes/res0.RData")

res05<-list()
for(b in 1:B){
  res05[[b]]<-run_one(0.5)
}
save(res05,file="C:/Users/j.galharret/Documents/GitHub/OMCOA/Rcodes/res05.RData")

res1<-list()
B<-1000
for(b in 1:B){
  res1[[b]]<-run_one(1)
}
save(res1,file="C:/Users/j.galharret/Documents/GitHub/OMCOA/Rcodes/res1.RData")

## ---------------------------------------------------------------- simulation 1 : t de 0 à 1

ts <- seq(0, 1, by = 0.1)
res_t <- do.call(rbind, lapply(ts, function(t) {
  R <- if (t %in% c(0, 1)) R_EXT else R_MID
  do.call(rbind, lapply(seq_len(R), function(r) {
    o <- run_one(t)
    data.frame(t = t, rep = r, method = c("MCOA", "OMCOA"),
               min_cos = c(min(o$MCOA$cos), min(o$OMCOA$cos)),
               err = c(o$MCOA$err, o$OMCOA$err),
               fit = c(o$MCOA$fit, o$OMCOA$fit))
  }))
}))

summ <- aggregate(cbind(min_cos, err, fit) ~ t + method, res_t, median)
print(reshape(summ, idvar = "t", timevar = "method", direction = "wide"), digits = 3)

## ---------------------------------------------------------------- simulation 2 : Sigma_k à t = 1
res_sig <- do.call(rbind, lapply(seq_len(R_EXT), function(r) {
  o <- run_one(1)
  do.call(rbind, lapply(c("MCOA", "OMCOA"), function(meth)
    data.frame(rep = r, method = meth,
               block = paste("Bloc", rep(seq_len(K), H)),
               dim = paste0("Dim", rep(seq_len(H), each = K)),
               estimate = as.vector(o[[meth]]$sigma),
               truth = as.vector(o$sig_true))))
}))

## ---------------------------------------------------------------- figures
col_meth <- c(MCOA = "#eb6834", OMCOA = "#2a78d6")
col_blk  <- c("Bloc 1" = "#2a78d6", "Bloc 2" = "#eb6834", "Bloc 3" = "#1baf7a")

q10 <- function(x) quantile(x, 0.1); q90 <- function(x) quantile(x, 0.9)

# Figure 1 : estimation des Sigma_k pour t = 1
d_sig <- aggregate(estimate ~ method + block + dim + truth, res_sig,
                   function(x) c(med = median(x), lo = q10(x), hi = q90(x)))
d_sig <- cbind(d_sig[, 1:4], as.data.frame(d_sig$estimate))
names(d_sig)[5:7] <- c("med", "lo", "hi")

g1 <- ggplot(d_sig, aes(dim, med, fill = block)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.8), width = 0) +
  geom_errorbar(aes(ymin = truth, ymax = truth), position = position_dodge(0.8),
                width = 0.6, linewidth = 1) +
  facet_wrap(~ method) +
  scale_fill_manual(values = col_blk, name = NULL) +
  labs(x = NULL, y = expression(hat(Sigma)[k] ~ "(médiane, 10–90 %)"),
       title = "Profils permutés (t = 1) : Σ1 = diag(.7,.5,.3), Σ2 = diag(.5,.3,.7), Σ3 = diag(.3,.7,.5)",
       subtitle = "Trait noir horizontal = vraie valeur") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave("fig_sigma_t1_R.png", g1, width = 10, height = 4.5, dpi = 200)

# Figure 2 : évolution avec t
d_t <- do.call(rbind, lapply(c("min_cos", "err"), function(v) {
  a <- aggregate(res_t[[v]], list(t = res_t$t, method = res_t$method),
                 function(x) c(med = median(x), lo = q10(x), hi = q90(x)))
  data.frame(a[, 1:2], as.data.frame(a$x),
             metric = ifelse(v == "min_cos", "Pire |cos| entre u_h vrai et estimé",
                             "Erreur moyenne sur Σ_k"))
}))
names(d_t)[3:5] <- c("med", "lo", "hi")

g2 <- ggplot(d_t, aes(t, med, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(values = col_meth, name = NULL) +
  scale_fill_manual(values = col_meth, name = NULL) +
  labs(x = "t   (0 : profil commun  →  1 : profils permutés)", y = NULL,
       title = "Quand les valeurs singulières de X se rapprochent, la MCOA décroche") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave("fig_courbe_t_R.png", g2, width = 10, height = 4.5, dpi = 200)
