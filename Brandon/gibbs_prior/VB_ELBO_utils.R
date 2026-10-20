# Local Gibbs VI helpers for VB_master.R. No installed packages are modified.

vb_support_cache <- function(F_mats, exact = FALSE) {
  stopifnot(length(F_mats) > 0L)
  template <- F_mats[[1L]]
  list(values = do.call(rbind, lapply(F_mats, as.vector)),
       dim = dim(template), exact = exact,
       n_shapes = phylodyn:::zigzag(nrow(template) + 1L))
}

vb_gibbs_state <- function(M, beta, cache) {
  stopifnot(is.finite(beta), beta > 0, identical(dim(M), cache$dim))
  delta <- sweep(cache$values, 2L, as.vector(M), "-")
  d2 <- rowSums(delta^2)
  log_mass <- -beta * d2
  offset <- max(log_mass)
  mass <- exp(log_mass - offset)
  probs <- mass / sum(mass)
  # Enumeration gives an exact sum; uniform proposals estimate that sum.
  log_Z <- offset + log(sum(mass))
  if (!cache$exact) {
    log_Z <- log_Z - log(length(mass)) + log(cache$n_shapes)
  }
  list(log_Z = log_Z, probs = probs,
       mean_F = matrix(as.numeric(crossprod(probs, cache$values)), cache$dim[1L]),
       mean_d2 = sum(probs * d2))
}

vb_mean_log_likelihood <- function(F, coal_times, sequences, R, rate) {
  stopifnot(R >= 1L, R == as.integer(R))
  n <- nrow(F) + 1L
  taxa <- names(sequences)
  stopifnot(length(taxa) == n, !anyDuplicated(taxa))
  lls <- replicate(R, {
    tree <- fmatrix::mytree_from_F(F, coal_times)
    # Explicit uniform assignments also ensure labels match the alignment.
    tree$tip.label <- sample(taxa, n, replace = FALSE)
    phangorn::pml(tree, sequences, bf = rep(0.25, 4L),
                  Q = c(1, 2, 1, 1, 2, 1), rate = rate)$logLik
  })
  mean(lls)
}

vb_gibbs_elbo_from_loglik <- function(samples, mean_loglik, M, beta, state,
                                     weights = NULL, baseline = 0) {
  S <- length(samples)
  stopifnot(S > 0L, length(mean_loglik) == S, all(is.finite(mean_loglik)))
  if (is.null(weights)) weights <- rep(1 / S, S)
  stopifnot(length(weights) == S, all(is.finite(weights)),
            all(weights >= 0), sum(weights) > 0)
  weights <- weights / sum(weights)
  values <- do.call(rbind, lapply(samples, as.vector))
  d2 <- rowSums(sweep(values, 2L, as.vector(M), "-")^2)
  cherries <- vapply(samples, function(F) {
    sum(diff(c(0, F[nrow(F), ])) == 2)
  }, numeric(1L))
  n <- nrow(M) + 1L
  log_q_shape <- -beta * d2 - state$log_Z
  log_q_labels <- cherries * log(2) - lfactorial(n)
  # Uniform prior over labeled ranked trees: 2^(n-1)/(n! (n-1)!).
  log_prior <- (n - 1L) * log(2) - lfactorial(n) - lfactorial(n - 1L)
  reward <- mean_loglik + log_prior - log_q_shape - log_q_labels
  score_M <- 2 * beta * sweep(values, 2L, as.vector(state$mean_F), "-")
  score_g <- beta * (state$mean_d2 - d2) # g = log(beta)
  # The baseline must be independent of this sampled batch. It changes only
  # gradient variance because the exact Gibbs score has expectation zero.
  advantage <- reward - baseline
  list(elbo = sum(weights * reward),
       grad_M = matrix(as.numeric(crossprod(weights * advantage, score_M)), nrow(M)),
       grad_g = sum(weights * advantage * score_g),
       rewards = reward)
}

vb_estimate_gibbs_elbo <- function(data, coal_times, samples, M, beta, state,
                                  num_tip_label_iters, rate = 1, baseline = 0) {
  mean_loglik <- vapply(samples, vb_mean_log_likelihood, numeric(1L),
                       coal_times = coal_times, sequences = data,
                       R = num_tip_label_iters, rate = rate)
  vb_gibbs_elbo_from_loglik(samples, mean_loglik, M, beta, state, baseline = baseline)
}
