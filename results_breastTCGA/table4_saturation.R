## =============================================================================
##  Table 4 (Appendix D): saturation of the effective rank
##  Sparse Gibbs sampler on data with the dimensions of breast.TCGA
##  (N = 150, p = 184, 200, 142), from gamma = 1 ("all") and gamma = 0 ("none").
##
##  Data
##    iid     : i.i.d. N(0, 1) noise, no structure            (N_REP replications)
##    hetero  : independent columns with the TCGA column sds,  (N_REP replications)
##              no common structure
##    tcga_std: TCGA, centred and reduced                       (N_SEED chains per init)
##    tcga_raw: TCGA, centred only (data of Section 6)          (N_SEED chains per init)
##
##  Each run is saved in OUT_DIR as soon as it ends: the script can be
##  interrupted and restarted, the finished runs are not recomputed.
## =============================================================================

library(BayesCanoFacto)

## ---- settings ---------------------------------------------------------------
R_MAX   <- 20
N_ITER  <- 6000
BURN    <- 2000
THIN    <- 10
N_REP   <- 5          # replications of the simulated data sets
N_SEED  <- 2          # chains per initialisation for the TCGA data
N_CORES <- max(1, parallel::detectCores() - 1)
OUT_DIR <- "results_table4"
dir.create(OUT_DIR, showWarnings = FALSE)
# Time: about 0.15 s per iteration for R_MAX = 20, i.e. 15 min per run;
# 2 x (2 N_REP + 2 N_SEED) = 28 runs with the values above.

## ---- data -------------------------------------------------------------------
if (requireNamespace("mixOmics", quietly = TRUE)) {
  data("breast.TCGA", package = "mixOmics")
} else {
  f <- file.path(tempdir(), "breast.TCGA.rda")
  download.file("https://raw.githubusercontent.com/mixOmicsTeam/mixOmics/master/data/breast.TCGA.rda",
                f, mode = "wb")
  load(f)
}
blocks <- c("mirna", "mrna", "protein")
tr  <- breast.TCGA$data.train[blocks]
N   <- nrow(tr[[1]])
sds <- lapply(tr, function(x) apply(x, 2, sd))

make_data <- function(type, rep) {
  set.seed(1000 + rep)
  Xs <- switch(type,
    iid      = lapply(sds, function(s) matrix(rnorm(N * length(s)), N)),
    hetero   = lapply(sds, function(s) sweep(matrix(rnorm(N * length(s)), N), 2, s, "*")),
    tcga_std = lapply(tr, as.matrix),
    tcga_raw = lapply(tr, as.matrix))
  lapply(Xs, function(x) scale(x, center = TRUE, scale = (type == "tcga_std")))
}

## ---- jobs -------------------------------------------------------------------
jobs <- rbind(
  expand.grid(type = c("iid", "hetero"), rep = seq_len(N_REP), init = c("all", "none"),
              stringsAsFactors = FALSE),
  expand.grid(type = c("tcga_std", "tcga_raw"), rep = seq_len(N_SEED), init = c("all", "none"),
              stringsAsFactors = FALSE))
jobs$file <- file.path(OUT_DIR, sprintf("%s_rep%d_%s.rds", jobs$type, jobs$rep, jobs$init))

run_job <- function(i) {
  jb <- jobs[i, ]
  if (file.exists(jb$file)) return(readRDS(jb$file))
  Xs <- make_data(jb$type, jb$rep)            # same data for the two initialisations
  set.seed(2000 + i)                          # different chain for each job
  d <- gibbs_sparse_omcoa(Xs, R = R_MAX, n_iter = N_ITER, burn = BURN, thin = THIN,
                          init_gamma = jb$init, verbose = FALSE)
  rk <- reff_posterior(d)$rank
  B  <- length(rk)
  out <- data.frame(
    type = jb$type, rep = jb$rep, init = jb$init,
    reff = mean(rk),
    reff_first = mean(rk[seq_len(B %/% 4)]),           # first quarter of the stored draws
    reff_last  = mean(rk[B - seq_len(B %/% 4) + 1]),   # last quarter: drift if different
    p_sat = mean(rk == R_MAX),
    t(colMeans(apply(d$G, c(1, 2), sum))),             # active dimensions per block
    time_min = d$time / 60)
  names(out)[8:10] <- paste0("n_", blocks)
  saveRDS(out, jb$file)
  out
}

res <- if (.Platform$OS.type == "windows" || N_CORES == 1) {
  lapply(seq_len(nrow(jobs)), run_job)
} else {
  parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = N_CORES, mc.preschedule = FALSE)
}
res <- do.call(rbind, res)
saveRDS(res, file.path(OUT_DIR, "all_runs.rds"))

## ---- table -----------------------------------------------------------------
vars <- c("reff", "reff_first", "reff_last", "p_sat", paste0("n_", blocks))
tab <- aggregate(res[vars], res[c("type", "init")], mean)
tab_sd <- aggregate(res[vars], res[c("type", "init")], sd)
tab$type <- factor(tab$type, c("iid", "hetero", "tcga_std", "tcga_raw"))
tab <- tab[order(tab$type, tab$init), ]
print(tab, digits = 3, row.names = FALSE)
cat("\nStandard deviations over the replications / chains:\n")
print(tab_sd[order(match(tab_sd$type, levels(tab$type)), tab_sd$init), ], digits = 2, row.names = FALSE)

## ---- LaTeX rows of Table 4 ---------------------------------------------------
lab <- c(iid = "i.i.d.\\ Gaussian noise", hetero = "heteroscedastic noise$^a$",
         tcga_std = "TCGA, standardised", tcga_raw = "TCGA, centred")
fmt <- function(type, init) {
  r <- tab[tab$type == type & tab$init == init, ]
  s <- tab_sd[tab_sd$type == type & tab_sd$init == init, ]
  c(sprintf("%.1f (%.1f)", r$reff, s$reff),
    sprintf("%.1f / %.1f / %.1f", r$n_mirna, r$n_mrna, r$n_protein))
}
for (ty in levels(tab$type)) {
  cat(lab[ty], "&", paste(c(fmt(ty, "all"), fmt(ty, "none")), collapse = " & "), "\\\\\n")
}
