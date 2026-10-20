compute_prior_info <- function(M,
                               b,
                               all_Fmats,
                               cache = NULL,
                               prob_tol = 1e-8) {
  if (!is.matrix(M) || nrow(M) != ncol(M)) {
    stop("M must be a square F-matrix-shaped matrix.")
  }
  if (!is.numeric(b) || length(b) != 1 || !is.finite(b) || b <= 0) {
    stop("b must be a single positive finite beta value.")
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

  if (is.null(cache)) {
    cache <- phylodyn:::precompute_tree_chain_distance_cache(all_Fmats)
  }
  distances <- phylodyn:::tree_chain_squared_distances(cache, M)
  log_weights <- -b * distances
  log_Z <- phylodyn:::compute_log_Z_est(beta = b, M = M, cache = cache)
  probs <- exp(log_weights - log_Z)

  prob_sum <- sum(probs)
  if (!isTRUE(all.equal(prob_sum, 1, tolerance = prob_tol))) {
    stop("Prior probabilities do not sum to 1; sum = ", prob_sum)
  }

  list(
    probs = probs,
    log_Z = log_Z,
    prob_sum = prob_sum,
    ess = 1 / sum(probs^2),
    support_size = length(all_Fmats),
    cache = cache
  )
}

sample_from_prior <- function(M,
                              b,
                              num_samps,
                              all_Fmats,
                              cache = NULL,
                              seed = NULL,
                              return_info = FALSE,
                              prob_tol = 1e-8) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (length(num_samps) != 1 || is.na(num_samps) || num_samps < 1 ||
      num_samps != as.integer(num_samps)) {
    stop("num_samps must be a positive integer.")
  }

  prior_info <- compute_prior_info(
    M = M,
    b = b,
    all_Fmats = all_Fmats,
    cache = cache,
    prob_tol = prob_tol
  )

  sample_indices <- sample.int(
    length(all_Fmats),
    size = num_samps,
    replace = TRUE,
    prob = prior_info$probs
  )
  samples <- all_Fmats[sample_indices]

  if (!return_info) {
    return(samples)
  }

  list(
    samples = samples,
    indices = sample_indices,
    probs = prior_info$probs,
    log_Z = prior_info$log_Z,
    prob_sum = prior_info$prob_sum,
    ess = prior_info$ess,
    support_size = prior_info$support_size
  )
}

write_prior_trees_to_file <- function(M,
                                      b,
                                      num_trees,
                                      all_Fmats,
                                      coal_times,
                                      file,
                                      cache = NULL,
                                      chunk_size = 5000,
                                      seed = NULL,
                                      indices_file = NULL,
                                      prob_tol = 1e-8) {
  if (length(num_trees) != 1 || is.na(num_trees) || num_trees < 1 ||
      num_trees != as.integer(num_trees)) {
    stop("num_trees must be a positive integer.")
  }
  if (length(chunk_size) != 1 || is.na(chunk_size) || chunk_size < 1 ||
      chunk_size != as.integer(chunk_size)) {
    stop("chunk_size must be a positive integer.")
  }
  if (length(file) != 1 || is.na(file) || file == "") {
    stop("file must be a non-empty output path.")
  }

  prior_info <- compute_prior_info(
    M = M,
    b = b,
    all_Fmats = all_Fmats,
    cache = cache,
    prob_tol = prob_tol
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (file.exists(file)) {
    unlink(file)
  }

  sample_indices <- integer(num_trees)
  start_idx <- 1L

  while (start_idx <= num_trees) {
    end_idx <- min(start_idx + chunk_size - 1L, num_trees)
    this_n <- end_idx - start_idx + 1L
    chunk_indices <- sample.int(
      length(all_Fmats),
      size = this_n,
      replace = TRUE,
      prob = prior_info$probs
    )
    sample_indices[start_idx:end_idx] <- chunk_indices

    chunk_trees <- lapply(
      all_Fmats[chunk_indices],
      function(Fmat) mytree_from_F(Fmat, coal_times)
    )
    class(chunk_trees) <- "multiPhylo"
    ape::write.tree(chunk_trees, file = file, append = start_idx > 1L)

    start_idx <- end_idx + 1L
  }

  if (!is.null(indices_file)) {
    saveRDS(sample_indices, file = indices_file)
  }

  list(
    file = file,
    indices_file = indices_file,
    num_trees = num_trees,
    beta = b,
    log_beta = log(b),
    log_Z = prior_info$log_Z,
    prob_sum = prior_info$prob_sum,
    ess = prior_info$ess,
    support_size = prior_info$support_size,
    unique_fmats = length(unique(sample_indices))
  )
}
