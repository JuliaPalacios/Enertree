build_fmatrix_from_categoricals <- function(n, categorical_vectors, seed = NULL) {
  if (length(n) != 1 || is.na(n) || n != as.integer(n) || n < 1) {
    stop("`n` must be a positive integer giving the F-matrix dimension.")
  }

  if (!is.list(categorical_vectors)) {
    stop("`categorical_vectors` must be a list of numeric vectors.")
  }

  expected_vectors <- max(0, n - 2)
  if (length(categorical_vectors) != expected_vectors) {
    stop(
      "`categorical_vectors` must have length ", expected_vectors,
      " for an ", n, " x ", n, " F-matrix."
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  F_mat <- matrix(0, nrow = n, ncol = n)
  diag(F_mat) <- seq_len(n) + 1

  if (n >= 2) {
    F_mat[2, 1] <- 1
  }

  breakpoints <- integer(expected_vectors)

  if (n <= 2) {
    attr(F_mat, "breakpoints") <- breakpoints
    return(F_mat)
  }

  for (row_idx in 3:n) {
    prev_row <- F_mat[row_idx - 1, seq_len(row_idx - 1)]
    probs <- categorical_vectors[[row_idx - 2]]

    if (!is.numeric(probs) || anyNA(probs) || length(probs) != (row_idx - 1)) {
      stop(
        "Categorical vector for row ", row_idx,
        " must be numeric, have no missing values, and have length ",
        row_idx - 1, "."
      )
    }

    if (any(probs < 0)) {
      stop("Categorical vector for row ", row_idx, " cannot contain negatives.")
    }

    valid_spots <- c(prev_row[1], diff(prev_row)) > 0
    masked_probs <- ifelse(valid_spots, probs, 0)

    if (sum(masked_probs) <= 0) {
      stop(
        "Categorical vector for row ", row_idx,
        " assigns zero mass to all valid breakpoint locations."
      )
    }

    masked_probs <- masked_probs / sum(masked_probs)
    breakpoint <- sample.int(row_idx - 1, size = 1, prob = masked_probs)

    F_mat[row_idx, seq_len(row_idx - 1)] <-
      prev_row - as.integer(seq_len(row_idx - 1) >= breakpoint)

    breakpoints[row_idx - 2] <- breakpoint
  }

  attr(F_mat, "breakpoints") <- breakpoints
  F_mat
}


sample_true_categorical_vectors <- function(n, alpha = 1, seed = NULL) {
  if (length(n) != 1 || is.na(n) || n != as.integer(n) || n < 3) {
    stop("`n` must be an integer at least 3.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  categorical_vectors <- vector("list", n - 2)

  for (row_idx in 3:n) {
    row_len <- row_idx - 1

    if (length(alpha) == 1) {
      alpha_row <- rep(alpha, row_len)
    } else if (length(alpha) == row_len) {
      alpha_row <- alpha
    } else {
      stop(
        "`alpha` must be a scalar or have length ", row_len,
        " for row ", row_idx, "."
      )
    }

    if (any(alpha_row <= 0)) {
      stop("All Dirichlet concentration parameters must be positive.")
    }

    draws <- rgamma(row_len, shape = alpha_row, rate = 1)
    categorical_vectors[[row_idx - 2]] <- draws / sum(draws)
  }

  categorical_vectors
}


sample_fmatrix_batch <- function(n,
                                 categorical_vectors,
                                 n_samples,
                                 seed = NULL) {
  if (length(n_samples) != 1 || is.na(n_samples) ||
      n_samples != as.integer(n_samples) || n_samples < 1) {
    stop("`n_samples` must be a positive integer.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  F_mats <- vector("list", n_samples)
  breakpoints <- matrix(NA_integer_, nrow = n_samples, ncol = max(0, n - 2))

  for (sample_idx in seq_len(n_samples)) {
    F_mats[[sample_idx]] <- build_fmatrix_from_categoricals(
      n = n,
      categorical_vectors = categorical_vectors
    )
    breakpoints[sample_idx, ] <- attr(F_mats[[sample_idx]], "breakpoints")
  }

  list(
    F_mats = F_mats,
    breakpoints = breakpoints
  )
}


extract_breakpoint_data_from_fmatrix <- function(F_mat) {
  if (!is.matrix(F_mat) || nrow(F_mat) != ncol(F_mat)) {
    stop("`F_mat` must be a square matrix.")
  }

  n <- nrow(F_mat)
  if (n < 3) {
    return(list(
      breakpoints = integer(0),
      valid_masks = list()
    ))
  }

  breakpoints <- integer(n - 2)
  valid_masks <- vector("list", n - 2)

  for (row_idx in 3:n) {
    prev_row <- F_mat[row_idx - 1, seq_len(row_idx - 1)]
    curr_row <- F_mat[row_idx, seq_len(row_idx - 1)]
    row_diff <- prev_row - curr_row

    if (!all(row_diff %in% c(0, 1))) {
      stop("Observed row differences are not compatible with a suffix drop.")
    }

    breakpoint <- match(TRUE, row_diff == 1)
    if (is.na(breakpoint)) {
      stop("Could not recover a breakpoint for row ", row_idx, ".")
    }

    expected_diff <- as.integer(seq_len(row_idx - 1) >= breakpoint)
    if (!all(row_diff == expected_diff)) {
      stop("Observed row differences are not a valid suffix pattern in row ", row_idx, ".")
    }

    valid_masks[[row_idx - 2]] <- c(prev_row[1], diff(prev_row)) > 0
    breakpoints[row_idx - 2] <- breakpoint
  }

  list(
    breakpoints = breakpoints,
    valid_masks = valid_masks
  )
}


extract_breakpoint_data_from_batch <- function(F_mats) {
  if (!is.list(F_mats) || length(F_mats) < 1) {
    stop("`F_mats` must be a non-empty list of F-matrices.")
  }

  n <- nrow(F_mats[[1]])
  if (!all(vapply(F_mats, function(F_mat) is.matrix(F_mat) &&
                  nrow(F_mat) == n && ncol(F_mat) == n, logical(1)))) {
    stop("All F-matrices must be square and have the same dimension.")
  }

  row_breakpoints <- vector("list", max(0, n - 2))
  row_valid_masks <- vector("list", max(0, n - 2))

  if (n <= 2) {
    return(list(
      n = n,
      n_samples = length(F_mats),
      row_breakpoints = row_breakpoints,
      row_valid_masks = row_valid_masks
    ))
  }

  extracted <- lapply(F_mats, extract_breakpoint_data_from_fmatrix)

  for (row_idx in 3:n) {
    row_breakpoints[[row_idx - 2]] <- vapply(
      extracted,
      function(x) x$breakpoints[[row_idx - 2]],
      integer(1)
    )
    row_valid_masks[[row_idx - 2]] <- do.call(
      rbind,
      lapply(extracted, function(x) x$valid_masks[[row_idx - 2]])
    )
  }

  list(
    n = n,
    n_samples = length(F_mats),
    row_breakpoints = row_breakpoints,
    row_valid_masks = row_valid_masks
  )
}


softmax_with_baseline_zero <- function(eta) {
  theta <- c(eta, 0)
  theta <- theta - max(theta)
  weights <- exp(theta)
  weights / sum(weights)
}


fit_row_categorical_mle <- function(observed_breakpoints,
                                    valid_mask,
                                    init_eta = NULL,
                                    method = "BFGS") {
  if (!is.numeric(observed_breakpoints) || anyNA(observed_breakpoints)) {
    stop("`observed_breakpoints` must be a numeric vector without missing values.")
  }

  if (!is.matrix(valid_mask)) {
    stop("`valid_mask` must be a logical matrix.")
  }

  valid_mask <- matrix(as.logical(valid_mask), nrow = nrow(valid_mask))
  n_obs <- length(observed_breakpoints)
  K <- ncol(valid_mask)

  if (n_obs != nrow(valid_mask)) {
    stop("`observed_breakpoints` and `valid_mask` must have matching numbers of rows.")
  }

  if (K < 2) {
    stop("Each row-wise categorical must have at least two entries.")
  }

  if (any(observed_breakpoints < 1 | observed_breakpoints > K)) {
    stop("Observed breakpoints fall outside the support of the categorical.")
  }

  observed_breakpoints <- as.integer(observed_breakpoints)

  if (any(!valid_mask[cbind(seq_len(n_obs), observed_breakpoints)])) {
    stop("Some observed breakpoints are not contained in their valid sets.")
  }

  if (is.null(init_eta)) {
    init_eta <- rep(0, K - 1)
  }

  if (length(init_eta) != (K - 1)) {
    stop("`init_eta` must have length ", K - 1, ".")
  }

  neg_loglik <- function(eta) {
    theta <- c(eta, 0)
    exp_theta <- exp(theta)
    loglik <- 0

    for (obs_idx in seq_len(n_obs)) {
      valid_idx <- valid_mask[obs_idx, ]
      loglik <- loglik +
        theta[observed_breakpoints[obs_idx]] -
        log(sum(exp_theta[valid_idx]))
    }

    -loglik
  }

  neg_grad <- function(eta) {
    theta <- c(eta, 0)
    exp_theta <- exp(theta)
    grad <- numeric(K)

    for (obs_idx in seq_len(n_obs)) {
      valid_idx <- valid_mask[obs_idx, ]
      denom <- sum(exp_theta[valid_idx])
      grad[observed_breakpoints[obs_idx]] <- grad[observed_breakpoints[obs_idx]] + 1
      grad[valid_idx] <- grad[valid_idx] - exp_theta[valid_idx] / denom
    }

    -grad[seq_len(K - 1)]
  }

  opt <- optim(
    par = init_eta,
    fn = neg_loglik,
    gr = neg_grad,
    method = method
  )

  list(
    probabilities = softmax_with_baseline_zero(opt$par),
    eta = opt$par,
    logLik = -opt$value,
    convergence = opt$convergence,
    counts = tabulate(observed_breakpoints, nbins = K)
  )
}


fit_categorical_vectors_from_fmats <- function(F_mats, method = "BFGS") {
  extracted <- extract_breakpoint_data_from_batch(F_mats)
  n <- extracted$n

  fitted_vectors <- vector("list", max(0, n - 2))
  fit_details <- vector("list", max(0, n - 2))

  if (n <= 2) {
    return(list(
      fitted_vectors = fitted_vectors,
      fit_details = fit_details,
      extracted_data = extracted
    ))
  }

  for (row_idx in 3:n) {
    fit_details[[row_idx - 2]] <- fit_row_categorical_mle(
      observed_breakpoints = extracted$row_breakpoints[[row_idx - 2]],
      valid_mask = extracted$row_valid_masks[[row_idx - 2]],
      method = method
    )
    fitted_vectors[[row_idx - 2]] <- fit_details[[row_idx - 2]]$probabilities
  }

  list(
    fitted_vectors = fitted_vectors,
    fit_details = fit_details,
    extracted_data = extracted
  )
}


loglik_fmatrix_under_categoricals <- function(F_mat, categorical_vectors) {
  extracted <- extract_breakpoint_data_from_fmatrix(F_mat)
  n <- nrow(F_mat)

  if (length(categorical_vectors) != max(0, n - 2)) {
    stop("`categorical_vectors` does not match the dimension of `F_mat`.")
  }

  if (n <= 2) {
    return(0)
  }

  loglik <- 0

  for (row_idx in 3:n) {
    probs <- categorical_vectors[[row_idx - 2]]
    valid_mask <- extracted$valid_masks[[row_idx - 2]]
    breakpoint <- extracted$breakpoints[[row_idx - 2]]
    denom <- sum(probs[valid_mask])

    if (denom <= 0) {
      stop("Encountered zero probability mass on the valid set in row ", row_idx, ".")
    }

    if (probs[breakpoint] <= 0) {
      return(-Inf)
    }

    loglik <- loglik + log(probs[breakpoint]) - log(denom)
  }

  loglik
}


loglik_batch_under_categoricals <- function(F_mats, categorical_vectors) {
  if (!is.list(F_mats) || length(F_mats) < 1) {
    stop("`F_mats` must be a non-empty list of F-matrices.")
  }

  vapply(
    F_mats,
    loglik_fmatrix_under_categoricals,
    categorical_vectors = categorical_vectors,
    numeric(1)
  )
}


mean_fmatrix <- function(F_mats, project_to_fspace = TRUE) {
  if (!is.list(F_mats) || length(F_mats) < 1) {
    stop("`F_mats` must be a non-empty list of F-matrices.")
  }

  mean_raw <- Reduce(`+`, F_mats) / length(F_mats)

  if (project_to_fspace && exists("nearby_Fmat", mode = "function")) {
    mean_projected <- nearby_Fmat(mean_raw)
  } else {
    mean_projected <- mean_raw
  }

  list(
    raw = mean_raw,
    projected = mean_projected
  )
}


count_unique_fmats <- function(F_mats) {
  if (!is.list(F_mats) || length(F_mats) < 1) {
    stop("`F_mats` must be a non-empty list of F-matrices.")
  }

  fmat_key <- function(F_mat) paste(F_mat[lower.tri(F_mat, diag = TRUE)], collapse = "|")
  length(unique(vapply(F_mats, fmat_key, character(1))))
}


flatten_categorical_vectors <- function(categorical_vectors) {
  if (!is.list(categorical_vectors)) {
    stop("`categorical_vectors` must be a list.")
  }

  do.call(
    rbind,
    lapply(seq_along(categorical_vectors), function(idx) {
      probs <- categorical_vectors[[idx]]
      data.frame(
        row = idx + 2L,
        category = seq_along(probs),
        probability = probs
      )
    })
  )
}


fit_and_generate_from_fmats <- function(F_mats,
                                        generated_sample_size = length(F_mats),
                                        seed = 1,
                                        method = "BFGS") {
  if (!is.list(F_mats) || length(F_mats) < 1) {
    stop("`F_mats` must be a non-empty list of F-matrices.")
  }

  n <- nrow(F_mats[[1]])
  fit <- fit_categorical_vectors_from_fmats(F_mats = F_mats, method = method)
  generated <- sample_fmatrix_batch(
    n = n,
    categorical_vectors = fit$fitted_vectors,
    n_samples = generated_sample_size,
    seed = seed
  )

  beast_mean <- mean_fmatrix(F_mats)
  generated_mean <- mean_fmatrix(generated$F_mats)

  beast_loglik <- loglik_batch_under_categoricals(F_mats, fit$fitted_vectors)
  generated_loglik <- loglik_batch_under_categoricals(generated$F_mats, fit$fitted_vectors)

  summary <- data.frame(
    n = n,
    beast_sample_size = length(F_mats),
    generated_sample_size = generated_sample_size,
    unique_beast_fmats = count_unique_fmats(F_mats),
    unique_generated_fmats = count_unique_fmats(generated$F_mats),
    mean_loglik_beast = mean(beast_loglik),
    sd_loglik_beast = sd(beast_loglik),
    mean_loglik_generated = mean(generated_loglik),
    sd_loglik_generated = sd(generated_loglik),
    mean_fmat_l2_distance = if (exists("distance_Fmat", mode = "function")) {
      distance_Fmat(beast_mean$projected, generated_mean$projected, dist = "l2")
    } else {
      sqrt(sum((beast_mean$projected - generated_mean$projected)^2))
    }
  )

  list(
    fitted_vectors = fit$fitted_vectors,
    fit_details = fit$fit_details,
    extracted_data = fit$extracted_data,
    generated_fmats = generated$F_mats,
    generated_breakpoints = generated$breakpoints,
    beast_mean = beast_mean,
    generated_mean = generated_mean,
    beast_loglik = beast_loglik,
    generated_loglik = generated_loglik,
    summary = summary
  )
}


summarize_categorical_recovery <- function(true_vectors, estimated_vectors) {
  if (!is.list(true_vectors) || !is.list(estimated_vectors) ||
      length(true_vectors) != length(estimated_vectors)) {
    stop("`true_vectors` and `estimated_vectors` must be lists of equal length.")
  }

  summary_rows <- lapply(seq_along(true_vectors), function(idx) {
    true_vec <- true_vectors[[idx]]
    est_vec <- estimated_vectors[[idx]]

    if (length(true_vec) != length(est_vec)) {
      stop("True and estimated vectors must match in length row by row.")
    }

    data.frame(
      row = idx + 2L,
      dimension = length(true_vec),
      l1_error = sum(abs(true_vec - est_vec)),
      l2_error = sqrt(sum((true_vec - est_vec)^2)),
      max_abs_error = max(abs(true_vec - est_vec))
    )
  })

  do.call(rbind, summary_rows)
}


run_categorical_recovery_experiment <- function(n = 10,
                                                n_samples = 10000,
                                                alpha = 1,
                                                seed = 1,
                                                true_categorical_vectors = NULL,
                                                method = "BFGS") {
  if (is.null(true_categorical_vectors)) {
    true_categorical_vectors <- sample_true_categorical_vectors(
      n = n,
      alpha = alpha,
      seed = seed
    )
    simulation_seed <- seed + 1
  } else {
    simulation_seed <- seed
  }

  simulated <- sample_fmatrix_batch(
    n = n,
    categorical_vectors = true_categorical_vectors,
    n_samples = n_samples,
    seed = simulation_seed
  )

  fit <- fit_categorical_vectors_from_fmats(
    F_mats = simulated$F_mats,
    method = method
  )

  comparison <- summarize_categorical_recovery(
    true_vectors = true_categorical_vectors,
    estimated_vectors = fit$fitted_vectors
  )

  list(
    true_categorical_vectors = true_categorical_vectors,
    simulated = simulated,
    fit = fit,
    comparison = comparison
  )
}


# Example:
# result <- run_categorical_recovery_experiment(n = 10, n_samples = 10000, seed = 1)
# result$comparison
# result$true_categorical_vectors[[1]]
# result$fit$fitted_vectors[[1]]
