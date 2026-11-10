## =============================================================================
##  Application of the sparse Bayesian OMCOA to the breast.TCGA data (mixOmics)
##  Package BayesCanoFacto
## =============================================================================
##
##  Data: breast.TCGA$data.train, N = 150 tumours, K = 3 blocks
##        miRNA (184), mRNA (200), proteins (142); PAM50 subtype (Basal, Her2,
##        LumA), NOT used in the fit, only for the interpretation.
##
##  Sections
##    0. Settings
##    1. Data and preprocessing
##    2. Implied prior of R_eff
##    3. Frequentist references: OMCOA, MCOA, Gavish-Donoho (and Hoff)
##    4. Choice of R_max: saturation check (short chains)
##    5. Main fit: C chains from gamma = "all" and gamma = "none"
##    6. MCMC diagnostics
##    7. Posterior of R_eff and of the activity patterns of the dimensions
##    8. Block weights and inclusion probabilities (aligned draws)
##    9. Interpretation: components vs subtypes, top loadings
##   10. Sensitivity to lambda, (a_theta, b_theta) and R_max
##   11. Saving
## =============================================================================

library(BayesCanoFacto)
library(ggplot2)

## -----------------------------------------------------------------------------
## 0. Settings
## -----------------------------------------------------------------------------
QUICK     <- FALSE        # TRUE: very short chains to test the script
R_MAX     <- if (QUICK) 10 else 40   # maximal rank (see section 4)
N_CHAINS  <- 4
N_ITER    <- if (QUICK) 300 else 6000
BURN      <- if (QUICK) 100 else 2000
THIN      <- if (QUICK) 1   else 10    # V is stored: thinning limits memory
N_CORES   <- max(1, min(N_CHAINS, parallel::detectCores() - 1))
RUN_SATURATION  <- TRUE   # section 4
RUN_HOFF        <- FALSE  # Hoff (2007) on the concatenated matrix: slow
RUN_SENSITIVITY <- TRUE   # section 10
OUT_DIR   <- "results_breastTCGA"
dir.create(OUT_DIR, showWarnings = FALSE)

# parallel lapply (forks on Linux/macOS, sequential on Windows)
plapply <- function(X, FUN, ...) {
  if (.Platform$OS.type == "windows" || N_CORES == 1) lapply(X, FUN, ...)
  else parallel::mclapply(X, FUN, ..., mc.cores = N_CORES, mc.set.seed = TRUE)
}
RNGkind("L'Ecuyer-CMRG")
set.seed(2026)

## -----------------------------------------------------------------------------
## 1. Data and preprocessing
## -----------------------------------------------------------------------------
if (requireNamespace("mixOmics", quietly = TRUE)) {
  data("breast.TCGA", package = "mixOmics")
} else {
  # same data, from the GitHub repository of mixOmics
  f <- file.path(tempdir(), "breast.TCGA.rda")
  download.file("https://raw.githubusercontent.com/mixOmicsTeam/mixOmics/master/data/breast.TCGA.rda",
                f, mode = "wb")
  load(f)
}
train   <- breast.TCGA$data.train
subtype <- train$subtype
blocks  <- c("mirna", "mrna", "protein")

# Columns centred only (no reduction): the model assumes an i.i.d. noise
# within each block, and the sampler is invariant to the scale of the blocks.
Xs <- lapply(train[blocks], function(x) scale(as.matrix(x), center = TRUE, scale = FALSE))
names(Xs) <- blocks
N <- nrow(Xs[[1]]); K <- length(Xs); p <- vapply(Xs, ncol, integer(1))

prepro <- data.frame(
  block = blocks, p = p,
  norm = vapply(Xs, function(x) sqrt(sum(x^2)), numeric(1)),
  sd_min = vapply(Xs, function(x) min(apply(x, 2, sd)), numeric(1)),
  sd_max = vapply(Xs, function(x) max(apply(x, 2, sd)), numeric(1)),
  gd_rank = vapply(Xs, function(x) gd_threshold(x)$rank, numeric(1)),  # per block
  row.names = NULL)
print(prepro, digits = 3)
table(subtype)

## -----------------------------------------------------------------------------
## 2. Implied prior of R_eff
## -----------------------------------------------------------------------------
# P(dimension r inactive in the K blocks) = prod_{i=0}^{K-1} (b+i)/(a+b+i),
# hence R_eff ~ Bin(R, 1 - that); prior correlation of the gamma_k^(r) of a
# dimension: 1/(a+b+1).
prior_reff <- function(R, K, a_theta = 1, b_theta = 1) {
  q0 <- prod((b_theta + 0:(K - 1)) / (a_theta + b_theta + 0:(K - 1)))
  list(p_gamma = a_theta / (a_theta + b_theta), p_active = 1 - q0,
       mean = R * (1 - q0), pmf = stats::dbinom(0:R, R, 1 - q0),
       corr = 1 / (a_theta + b_theta + 1))
}
prior_tab <- do.call(rbind, lapply(list(c(1, 1), c(1, 3), c(0.5, 0.5)), function(ab) {
  pr <- prior_reff(R_MAX, K, ab[1], ab[2])
  data.frame(a_theta = ab[1], b_theta = ab[2], p_gamma = pr$p_gamma,
             p_active = pr$p_active, E_Reff = pr$mean,
             q05 = qbinom(0.05, R_MAX, pr$p_active), q95 = qbinom(0.95, R_MAX, pr$p_active),
             corr = pr$corr)
}))
print(prior_tab, digits = 3)

## -----------------------------------------------------------------------------
## 3. Frequentist references
## -----------------------------------------------------------------------------
# OMCOA on the blocks divided by their norms (= default init of the sampler),
# weights brought back to the original scale; used as reference for alignment.
nrm   <- prepro$norm
fit_omcoa <- omcoa_als(Map(`/`, Xs, nrm), R = R_MAX)
fit_omcoa$S <- Map(`*`, fit_omcoa$S, nrm)
fit_mcoa  <- mcoa(Map(`/`, Xs, nrm), R = R_MAX)

Xcat <- do.call(cbind, Xs)
gd   <- gd_threshold(Xcat)
cat("Gavish-Donoho rank of the concatenated matrix:", gd$rank, "\n")
if (RUN_HOFF) {
  fit_hoff <- hoff_svd(Xcat, R = min(80, N - 1), n_iter = if (QUICK) 200 else 3000,
                       burn = if (QUICK) 50 else 500, verbose = TRUE)
  cat("Hoff (2007), posterior mode of the rank:", fit_hoff$mode, "\n")
}

## -----------------------------------------------------------------------------
## 4. Choice of R_max: saturation check
## -----------------------------------------------------------------------------
# R_eff counts the dimensions active in AT LEAST ONE block: a dimension
# specific to one block counts. If P(R_eff = R_max | X) is close to 1, R_max
# is too small for R_eff, but the number of dimensions active in 2 or 3
# blocks may already be stable. Short chains for increasing R_max.
patterns <- function(G) {
  # G: B x K x R -> B x (K) matrix of the number of dims active in exactly j blocks
  nb <- apply(G, c(1, 3), sum)                          # B x R
  t(apply(nb, 1, function(z) tabulate(z + 1, nbins = K + 1)))[, -1, drop = FALSE]
}
if (RUN_SATURATION) {
  R_grid <- if (QUICK) c(10, 20) else c(10, 20, 40, 60)
  sat <- plapply(R_grid, function(R) {
    d <- gibbs_sparse_omcoa(Xs, R = R, n_iter = if (QUICK) 200 else 1500,
                            burn = if (QUICK) 100 else 750, thin = 5, verbose = FALSE)
    pt <- patterns(d$G)
    data.frame(R_max = R, E_Reff = mean(reff_posterior(d)$rank),
               P_sat = mean(reff_posterior(d)$rank == R),
               n_1block = mean(pt[, 1]), n_2blocks = mean(pt[, 2]),
               n_3blocks = mean(pt[, 3]))
  })
  sat <- do.call(rbind, sat)
  print(sat, digits = 3)
}

## -----------------------------------------------------------------------------
## 5. Main fit
## -----------------------------------------------------------------------------
init_g <- rep(c("all", "none"), length.out = N_CHAINS)
chains <- plapply(seq_len(N_CHAINS), function(c)
  gibbs_sparse_omcoa(Xs, R = R_MAX, n_iter = N_ITER, burn = BURN, thin = THIN,
                     init_gamma = init_g[c], verbose = FALSE))
names(chains) <- paste0("chain", seq_len(N_CHAINS), "_", init_g)
cat("Time per chain (min):", round(sapply(chains, `[[`, "time") / 60, 1), "\n")

## -----------------------------------------------------------------------------
## 6. MCMC diagnostics
## -----------------------------------------------------------------------------
dg <- mcmc_diagnostics(chains, U_ref = fit_omcoa$U)
print(dg)
pdf(file.path(OUT_DIR, "diagnostics.pdf"), width = 8, height = 9)
plot(dg)
dev.off()
# indicators that disagree between chains
subset(dg$gamma, max_diff > 0.2)

## -----------------------------------------------------------------------------
## 7. Posterior of R_eff and of the activity patterns
## -----------------------------------------------------------------------------
# Aligned and pooled draws
bind_draws <- function(L) {
  b3 <- function(nm) {
    A <- lapply(L, `[[`, nm); d <- dim(A[[1]])
    out <- array(NA_real_, c(sum(sapply(A, function(a) dim(a)[1])), d[-1]))
    i <- 0
    for (a in A) { n <- dim(a)[1]; out[i + seq_len(n), , ] <- a; i <- i + n }
    out
  }
  out <- list(U = b3("U"), S = b3("S"), G = b3("G"),
              tau = do.call(rbind, lapply(L, `[[`, "tau")),
              theta = do.call(rbind, lapply(L, `[[`, "theta")))
  out$V <- lapply(seq_along(L[[1]]$V), function(k) {
    A <- lapply(L, function(d) d$V[[k]]); d <- dim(A[[1]])
    o <- array(NA_real_, c(sum(sapply(A, function(a) dim(a)[1])), d[-1])); i <- 0
    for (a in A) { n <- dim(a)[1]; o[i + seq_len(n), , ] <- a; i <- i + n }
    o
  })
  names(out$V) <- names(L[[1]]$V)
  out
}
al    <- lapply(chains, align_draws, U_ref = fit_omcoa$U)
draws <- bind_draws(al)
names(draws$V) <- blocks
B_tot <- dim(draws$S)[1]

# R_eff (invariant to the labels: no alignment needed)
reff <- reff_posterior(draws)
print(round(reff$post[reff$post > 0], 3))

# number of dimensions active in exactly 1, 2, 3 blocks
pt <- patterns(draws$G)
colnames(pt) <- paste0(1:K, "_block", c("", "s", "s"))
print(t(apply(pt, 2, function(z) c(mean = mean(z), quantile(z, c(0.025, 0.5, 0.975))))))

# number of dimensions for each subset of blocks (2^K - 1 subsets)
subsets <- expand.grid(rep(list(0:1), K))[-1, ]
names(subsets) <- blocks
subset_count <- t(apply(draws$G, 1, function(Gb) {
  code <- colSums(Gb * 2^(0:(K - 1)))                  # code of each dimension
  tabulate(code, nbins = 2^K - 1)
}))
pattern_tab <- data.frame(subsets,
                          mean = colMeans(subset_count),
                          lower = apply(subset_count, 2, quantile, 0.025),
                          upper = apply(subset_count, 2, quantile, 0.975))
pattern_tab$label <- apply(subsets, 1, function(z) paste(blocks[z == 1], collapse = "+"))
print(pattern_tab[order(-pattern_tab$mean), c("label", "mean", "lower", "upper")],
      row.names = FALSE, digits = 3)

## -----------------------------------------------------------------------------
## 8. Block weights and inclusion probabilities (aligned draws)
## -----------------------------------------------------------------------------
ss <- summary_sigma(draws)                                 # block, dim, estimate, CI, pip
ss$block <- factor(blocks[ss$block], levels = blocks)
ss$ratio <- ss$estimate / nrm[as.integer(ss$block)]        # sigma / ||X_k||
# order the dimensions by total relative weight
ord <- order(-tapply(ss$ratio^2, ss$dim, sum))
ss$dim_ord <- match(ss$dim, ord)

g_pip <- ggplot(ss, aes(dim_ord, block, fill = pip)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  labs(x = "dimension (ordered by relative weight)", y = NULL,
       fill = "P(active | X)") + theme_minimal()
g_sig <- ggplot(ss, aes(dim_ord, ratio^2, colour = block)) +
  geom_point(position = position_dodge(0.6)) +
  geom_errorbar(aes(ymin = (lower / nrm[as.integer(block)])^2,
                    ymax = (upper / nrm[as.integer(block)])^2),
                width = 0, position = position_dodge(0.6)) +
  labs(x = "dimension (ordered by relative weight)",
       y = expression((sigma[k]^(r) / "||" * X[k] * "||")^2)) + theme_minimal()
ggsave(file.path(OUT_DIR, "pip_heatmap.pdf"), g_pip, width = 10, height = 2.5)
ggsave(file.path(OUT_DIR, "weights.pdf"), g_sig, width = 10, height = 4)

# comparison with OMCOA and MCOA (relative weights, original order of OMCOA)
cmp <- rbind(
  data.frame(method = "sparse Gibbs", ss[, c("block", "dim")], ratio2 = ss$ratio^2),
  data.frame(method = "OMCOA", block = rep(blocks, each = R_MAX), dim = rep(1:R_MAX, K),
             ratio2 = as.vector(t(sapply(seq_len(K), function(k) (fit_omcoa$S[[k]] / nrm[k])^2)))))
ggsave(file.path(OUT_DIR, "weights_vs_omcoa.pdf"),
       ggplot(cmp, aes(dim, ratio2, colour = method)) + geom_point() +
         facet_wrap(~block, ncol = 1) + theme_minimal(), width = 10, height = 6)

## -----------------------------------------------------------------------------
## 9. Interpretation: components vs subtypes, top loadings
## -----------------------------------------------------------------------------
# Posterior of R^2(u^(r) ~ subtype) for each aligned dimension
# u is centred: R^2 = ||projection on the group means||^2 / ||u||^2
Z  <- model.matrix(~ subtype - 1); Pz <- Z %*% solve(crossprod(Z), t(Z))
r2 <- apply(draws$U, c(1, 3), function(u) sum((Pz %*% u)^2) / sum(u^2))      # B x R
r2_tab <- data.frame(dim = 1:R_MAX, dim_ord = match(1:R_MAX, ord),
                     r2 = colMeans(r2), lower = apply(r2, 2, quantile, 0.025),
                     upper = apply(r2, 2, quantile, 0.975),
                     pattern = sapply(1:R_MAX, function(r)
                       paste(blocks[pip(draws)[, r] > 0.5], collapse = "+")))
r2_tab <- r2_tab[order(-r2_tab$r2), ]
print(head(r2_tab, 10), digits = 2, row.names = FALSE)

# Components most associated with the subtypes
top2 <- r2_tab$dim[1:2]
U_mean <- apply(draws$U, c(2, 3), mean)
g_u <- ggplot(data.frame(u1 = U_mean[, top2[1]], u2 = U_mean[, top2[2]], subtype),
              aes(u1, u2, colour = subtype)) + geom_point() +
  labs(x = paste0("dimension ", top2[1], " (", r2_tab$pattern[1], ")"),
       y = paste0("dimension ", top2[2], " (", r2_tab$pattern[2], ")")) + theme_minimal()
ggsave(file.path(OUT_DIR, "components_subtype.pdf"), g_u, width = 6, height = 5)

# Top 10 variables per active block, for the 3 dimensions most related to subtype
top_load <- do.call(rbind, lapply(r2_tab$dim[1:3], function(r) do.call(rbind,
  lapply(seq_len(K), function(k) {
    if (pip(draws)[k, r] < 0.5) return(NULL)
    v <- colMeans(draws$V[[k]][, , r])
    i <- order(-abs(v))[1:10]
    data.frame(dim = r, block = blocks[k], variable = colnames(Xs[[k]])[i],
               loading = round(v[i], 3))
  }))))
print(top_load, row.names = FALSE)

## -----------------------------------------------------------------------------
## 10. Sensitivity (one chain per setting, one parameter at a time)
## -----------------------------------------------------------------------------
if (RUN_SENSITIVITY) {
  settings <- list(
    list(name = "default",        lambda = 1,   a = 1,   b = 1,   R = R_MAX),
    list(name = "lambda=0.1",     lambda = 0.1, a = 1,   b = 1,   R = R_MAX),
    list(name = "lambda=10",      lambda = 10,  a = 1,   b = 1,   R = R_MAX),
    list(name = "a,b=1,3",        lambda = 1,   a = 1,   b = 3,   R = R_MAX),
    list(name = "a,b=0.5,0.5",    lambda = 1,   a = 0.5, b = 0.5, R = R_MAX),
    list(name = "slab=noise",     lambda = 1,   a = 1,   b = 1,   R = R_MAX, slab = "noise"))
  sens <- plapply(settings, function(s) {
    d <- gibbs_sparse_omcoa(Xs, R = s$R, lambda = s$lambda, a_theta = s$a, b_theta = s$b,
                            slab = if (is.null(s$slab)) "norm" else s$slab,
                            n_iter = if (QUICK) 200 else 3000, burn = if (QUICK) 100 else 1000,
                            thin = 10, verbose = FALSE)
    pt <- patterns(d$G)
    data.frame(setting = s$name, E_Reff = mean(reff_posterior(d)$rank),
               n_1block = mean(pt[, 1]), n_2blocks = mean(pt[, 2]),
               n_3blocks = mean(pt[, 3]),
               n_active_mirna = mean(rowSums(d$G[, 1, , drop = FALSE])),
               n_active_mrna = mean(rowSums(d$G[, 2, , drop = FALSE])),
               n_active_protein = mean(rowSums(d$G[, 3, , drop = FALSE])))
  })
  sens <- do.call(rbind, sens)
  print(sens, digits = 3, row.names = FALSE)
}

## -----------------------------------------------------------------------------
## 11. Saving
## -----------------------------------------------------------------------------
save(prepro, prior_tab, gd, fit_omcoa, fit_mcoa, dg, reff, pattern_tab, ss, r2_tab,
     top_load, list = c(if (RUN_SATURATION) "sat", if (RUN_SENSITIVITY) "sens",
                        if (RUN_HOFF) "fit_hoff"),
     file = file.path(OUT_DIR, "results.RData"))
saveRDS(chains, file.path(OUT_DIR, "chains.rds"))   # large file
