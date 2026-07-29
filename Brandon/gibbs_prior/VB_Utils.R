sample_from_prior <- function(M,
                              b,
                              num_samps,
                              all_Fmats,
                              seed = NULL,
                              return_info = FALSE,
                              prob_tol = 1e-8) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.matrix(M) || nrow(M) != ncol(M)) {
    stop("M must be a square F-matrix-shaped matrix.")
  }
  if (!is.numeric(b) || length(b) != 1 || !is.finite(b) || b <= 0) {
    stop("b must be a single positive finite beta value.")
  }
  if (length(num_samps) != 1 || is.na(num_samps) || num_samps < 1 ||
      num_samps != as.integer(num_samps)) {
    stop("num_samps must be a positive integer.")
  }
  if (!is.list(all_Fmats) || length(all_Fmats) == 0) {
    stop("all_Fmats must be a non-empty list of F-matrices.")
  }

  expected_dim <- dim(M)
  valid_dims <- vapply(
    all_Fmats,
    function(Fmat) is.matrix(Fmat) && identical(dim(Fmat), expected_dim),
    logical(1)
  )
  if (!all(valid_dims)) {
    stop("Every F-matrix in all_Fmats must have the same dimensions as M.")
  }

  cache <- phylodyn:::precompute_tree_chain_distance_cache(all_Fmats)
  distances <- phylodyn:::tree_chain_squared_distances(cache, M)
  log_weights <- -b * distances
  log_Z <- phylodyn:::compute_log_Z_est(beta = b, M = M, cache = cache)
  probs <- exp(log_weights - log_Z)

  prob_sum <- sum(probs)
  if (!isTRUE(all.equal(prob_sum, 1, tolerance = prob_tol))) {
    stop("Prior probabilities do not sum to 1; sum = ", prob_sum)
  }

  sample_indices <- sample.int(
    length(all_Fmats),
    size = num_samps,
    replace = TRUE,
    prob = probs
  )
  samples <- all_Fmats[sample_indices]

  if (!return_info) {
    return(samples)
  }

  list(
    samples = samples,
    indices = sample_indices,
    probs = probs,
    log_Z = log_Z,
    prob_sum = prob_sum,
    ess = 1 / sum(probs^2),
    support_size = length(all_Fmats)
  )
}
