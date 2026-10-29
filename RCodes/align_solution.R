## -----------------------------------------------------------------------------
## Alignment of a multiblock solution (MCOA, OMCOA, or a Gibbs draw) on a
## reference U_ref: permutation of the dimensions + sign of the axes.
##
## Rule (Section "Post-processing" of the manuscript):
##   C        = t(U_ref) %*% U                   (cosines, columns of U are unit-norm)
##   pi_hat   = argmax_pi  sum_r |C[r, pi(r)]|
##   s_r      = sign(C[r, pi_hat(r)])
##   columns of U and V_k are reordered by pi_hat and multiplied by s_r,
##   sigma_k is reordered by pi_hat (a joint change of sign of u^(r) and of all
##   v_k^(r) leaves sigma_k^(r) = u' X_k v_k unchanged).
##
## The number of estimated dimensions R_est may differ from the number of
## reference dimensions R_ref (e.g. 5 estimated dimensions vs. a true U with 3
## columns):
##   - R_est >= R_ref : each reference axis receives one estimated axis; the
##     R_est - R_ref unmatched estimated axes are put at the end, in decreasing
##     order of sum_k sigma_k^2;
##   - R_est <  R_ref : reference axes without a match get a column of NA.
## -----------------------------------------------------------------------------

## Optimal assignment maximizing sum |C[r, pi(r)]|, with nrow(C) <= ncol(C).
## Uses clue::solve_LSAP (Hungarian algorithm) when available, otherwise an
## exhaustive search (fine up to about 8 dimensions).
.best_assignment <- function(A) {
  stopifnot(nrow(A) <= ncol(A))
  if (requireNamespace("clue", quietly = TRUE))
    return(as.integer(clue::solve_LSAP(A, maximum = TRUE)))
  n <- nrow(A); m <- ncol(A)
  if (choose(m, n) * factorial(n) > 5e6)
    stop("Too many dimensions for the exhaustive search: install.packages('clue')")
  best <- NULL; best_val <- -Inf
  rec <- function(r, used, cur, val) {
    if (r > n) {
      if (val > best_val) { best_val <<- val; best <<- cur }
      return(invisible())
    }
    for (j in setdiff(seq_len(m), used))
      rec(r + 1, c(used, j), c(cur, j), val + A[r, j])
  }
  rec(1, integer(0), integer(0), 0)
  best
}

#' Align a multiblock solution on a reference.
#'
#' @param U      N x R_est matrix of common components.
#' @param V      list of K matrices p_k x R_est (block loadings), or NULL.
#' @param sigma  K x R_est matrix of block weights sigma_k^(r)
#'               (a list of K vectors or of K diagonal matrices is also accepted),
#'               or NULL.
#' @param U_ref  N x R_ref reference matrix (true U, or the OMCOA solution).
#' @param positive_sigma if TRUE, a negative sigma_k^(r) is made positive by
#'               changing the sign of v_k^(r) (useful for MCOA, whose block
#'               weights u' X_k v_k are not constrained to be positive).
#' @param normalize if TRUE, the columns of U and U_ref are normalized before
#'               computing the cosines (MCOA components may not be unit-norm).
#' @return a list with the aligned U, V, sigma, and
#'   perm     : perm[r] = estimated dimension matched with reference axis r
#'              (NA if none), followed by the unmatched estimated dimensions;
#'   signs    : signs applied, in the output order;
#'   cosines  : |cos(u_ref^(r), u_hat^(perm(r)))| for r = 1..R_ref.
align_solution <- function(U, V = NULL, sigma = NULL, U_ref,
                           positive_sigma = TRUE, normalize = TRUE) {
  U <- as.matrix(U); U_ref <- as.matrix(U_ref)
  stopifnot(nrow(U) == nrow(U_ref))
  R_est <- ncol(U); R_ref <- ncol(U_ref)

  ## block weights as a K x R_est matrix
  if (!is.null(sigma)) {
    if (is.list(sigma))
      sigma <- do.call(rbind, lapply(sigma, function(s)
        if (is.matrix(s)) diag(s) else as.numeric(s)))
    sigma <- as.matrix(sigma)
    stopifnot(ncol(sigma) == R_est)
  }

  ## optional: positive block weights (sign absorbed in v_k^(r))
  if (positive_sigma && !is.null(sigma)) {
    for (k in seq_len(nrow(sigma))) {
      neg <- which(sigma[k, ] < 0)
      if (length(neg)) {
        if (!is.null(V)) V[[k]][, neg] <- -V[[k]][, neg]
        sigma[k, neg] <- -sigma[k, neg]
      }
    }
  }

  ## cosines between reference and estimated axes
  nrm <- function(M) sweep(M, 2, pmax(sqrt(colSums(M^2)), .Machine$double.eps), "/")
  C <- if (normalize) crossprod(nrm(U_ref), nrm(U)) else crossprod(U_ref, U)  # R_ref x R_est
  A <- abs(C)

  ## optimal assignment
  if (R_ref <= R_est) {
    match_ref <- .best_assignment(A)                       # ref r -> est match_ref[r]
  } else {
    est_to_ref <- .best_assignment(t(A))                   # est j -> ref est_to_ref[j]
    match_ref <- rep(NA_integer_, R_ref)
    match_ref[est_to_ref] <- seq_len(R_est)
  }

  ## unmatched estimated dimensions, by decreasing explained variance
  extra <- setdiff(seq_len(R_est), match_ref)
  if (length(extra) && !is.null(sigma))
    extra <- extra[order(colSums(sigma[, extra, drop = FALSE]^2), decreasing = TRUE)]
  perm <- c(match_ref, extra)

  ## signs: matched axes take the sign of the cosine, extra axes are left as is
  signs <- rep(1, length(perm))
  for (r in seq_len(R_ref))
    if (!is.na(perm[r])) signs[r] <- ifelse(C[r, perm[r]] < 0, -1, 1)

  ## reorder (NA columns for unmatched reference axes)
  take <- function(M) {
    out <- matrix(NA_real_, nrow(M), length(perm))
    ok <- !is.na(perm)
    out[, ok] <- sweep(M[, perm[ok], drop = FALSE], 2, signs[ok], "*")
    out
  }
  U_al <- take(U)
  V_al <- if (is.null(V)) NULL else lapply(V, take)
  sigma_al <- NULL
  if (!is.null(sigma)) {
    sigma_al <- matrix(NA_real_, nrow(sigma), length(perm))
    ok <- !is.na(perm)
    sigma_al[, ok] <- sigma[, perm[ok], drop = FALSE]     # no sign change
    rownames(sigma_al) <- rownames(sigma)
  }

  cosines <- sapply(seq_len(R_ref), function(r)
    if (is.na(perm[r])) NA_real_ else A[r, perm[r]])

  list(U = U_al, V = V_al, sigma = sigma_al,
       perm = perm, signs = signs, cosines = cosines)
}

