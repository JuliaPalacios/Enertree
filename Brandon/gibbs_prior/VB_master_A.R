# Monte Carlo VI on labeled ranked trees represented by A-matrices.
#
# q_Theta(A) proportional to exp(<Theta,A>), A a valid binary A-matrix.
# This is the real-center squared-distance family with its redundant beta
# absorbed into Theta. beta_init controls initialization only.
#
# Gradient: Cov_q(A, log p(Y|A,tau) - <Theta,A>).
# RNNI Metropolis-Hastings supplies the draws; no enumeration or log Z is
# used in the gradient. Likelihoods keep the labels encoded by A fixed.
#
# The displayed ELBO/entropy use a uniform-reference estimate of log Z.
# They are Monte Carlo diagnostics, not certified lower bounds. Inspect
# reference ESS, within-chain ESS and between-chain differences together.
#
# From the repository root:
# source("Brandon/gibbs_prior/VB_master_A.R")
# fit <- run_VB_master_A(num_tips=6, num_grad_desc_steps=20,
#                       output_dir="/private/tmp/VB_master_A_test")
# Rscript Brandon/gibbs_prior/VB_master_A.R [output_dir] runs defaults.

vb_A_source_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
if (is.null(vb_A_source_file)) {
  vb_A_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(vb_A_file_arg)) vb_A_source_file <- sub("^--file=", "", vb_A_file_arg[1L])
}
# A-matrix validation, RNNI sampling and gradients are provided by phylodyn.
# Install phylodyn with the updated R/VB_Utils.R before running this script.

vb_A_from_tree <- function(tree, leaf_order = tree$tip.label) {
  if (anyDuplicated(leaf_order) || !setequal(leaf_order, tree$tip.label)) {
    stop("leaf_order must contain each tree label exactly once.")
  }
  A <- phylodyn:::gen_Amat(tree)
  A <- A[, match(leaf_order, tree$tip.label), drop = FALSE]
  dimnames(A) <- list(paste0("rank_", seq_len(nrow(A))), leaf_order)
  if (!phylodyn:::a_matrix_valid(A)) stop("The tree did not produce a valid ranked A-matrix.")
  A
}

vb_A_to_tree <- function(A, coal_times, leaf_order = colnames(A), check = TRUE) {
  if (check && !phylodyn:::a_matrix_valid(A)) stop("A must be a valid binary A-matrix.")
  n <- ncol(A)
  if (length(leaf_order) != n || anyNA(leaf_order) || anyDuplicated(leaf_order)) {
    stop("leaf_order must contain one unique label per column of A.")
  }
  if (length(coal_times) != n - 1L || any(!is.finite(coal_times)) ||
      any(coal_times <= 0) || is.unsorted(coal_times, strictly = TRUE)) {
    stop("coal_times must contain n-1 strictly increasing positive times.")
  }
  blocks <- lapply(seq_len(n), function(j) as.numeric(seq_len(n) == j))
  nodes <- seq_len(n)
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  cursor <- 0L
  for (i in seq.int(n - 1L, 1L)) {
    selected <- which(vapply(blocks, function(B) any(B > 0 & A[i, ] > 0), logical(1L)))
    if (length(selected) != 2L ||
        !all(blocks[[selected[1L]]] + blocks[[selected[2L]]] == A[i, ])) {
      stop("An A-matrix row must merge exactly two active clades.")
    }
    parent <- n + i
    edge[cursor + 1:2, ] <- cbind(parent, nodes[selected])
    cursor <- cursor + 2L
    blocks <- c(blocks[-selected], list(as.numeric(A[i, ])))
    nodes <- c(nodes[-selected], parent)
  }
  heights <- c(rep(0, n), rev(coal_times))
  tree <- list(edge = edge, edge.length = heights[edge[, 1L]] - heights[edge[, 2L]],
               tip.label = as.character(leaf_order), Nnode = n - 1L)
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

vb_A_uniform_draw <- function(leaf_order) {
  # Uniform labeled ranked genealogy: choose each coalescing pair uniformly.
  n <- length(leaf_order)
  if (n < 2L || anyNA(leaf_order) || anyDuplicated(leaf_order)) stop("Invalid leaf labels.")
  blocks <- lapply(seq_len(n), function(j) as.integer(seq_len(n) == j))
  A <- matrix(0L, nrow = n - 1L, ncol = n,
              dimnames = list(paste0("rank_", seq_len(n - 1L)), leaf_order))
  for (i in seq.int(n - 1L, 1L)) {
    selected <- sample.int(length(blocks), 2L, replace = FALSE)
    A[i, ] <- blocks[[selected[1L]]] + blocks[[selected[2L]]]
    blocks <- c(blocks[-selected], list(A[i, ]))
  }
  A
}

run_VB_master_A <- function(num_tips = 10L, rate = 0.01, step_size = 0.1,
                            num_samps = 200L, num_unif_samples = 10000L,
                            num_grad_desc_steps = 5000L, take_every = 20L,
                            burnin = 1000L, num_chains = 2L, beta_init = 0.01,
                            init_method = c("upgma", "caterpillar"), seed = 42L,
                            sequences = NULL, coal_times = NULL, true_tree = NULL,
                            log_every = 25L, max_likelihood_cache = 100000L,
                            output_dir = NULL, make_plots = interactive()) {
  suppressPackageStartupMessages({
    library(ape); library(phylodyn); library(phyclust)
    library(phangorn); library(phylotools); library(fmatrix)
  })
  init_method <- match.arg(init_method)
  counts <- c(num_tips, num_samps, num_unif_samples, num_grad_desc_steps,
              take_every, burnin, num_chains, log_every, max_likelihood_cache)
  if (any(!is.finite(counts)) || any(counts != as.integer(counts)) ||
      num_tips < 3L || num_samps < 2L || num_unif_samples < 2L ||
      num_grad_desc_steps < 1L || take_every < 1L || burnin < 0L ||
      num_chains < 1L || log_every < 1L || max_likelihood_cache < 0L) {
    stop("Invalid integer configuration values.")
  }
  if (!is.finite(rate) || rate <= 0 || !is.finite(step_size) || step_size <= 0 ||
      !is.finite(beta_init) || beta_init < 0) stop("Invalid rate, step_size or beta_init.")
  if (!is.null(seed)) set.seed(seed)
  if (!make_plots) {
    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  start_time <- proc.time()[["elapsed"]]
  if (is.null(sequences)) {
    init_results <- phylodyn:::generate_true_M_and_data(num_tips, rate = rate, seq_len = 1000)
    sequences <- init_results$sequences
    true_tree <- init_results$M_true_tree
  }
  if (is.null(coal_times)) {
    if (is.null(true_tree)) stop("Supply fixed coal_times when providing sequences without a true_tree.")
    intervals <- ape::coalescent.intervals(true_tree)$interval.length
    if (any(!is.finite(intervals)) || any(intervals <= 0)) stop("Invalid coalescent intervals.")
    coal_times <- cumsum(intervals)
  }
  leaf_order <- names(sequences)
  if (length(leaf_order) != num_tips || anyNA(leaf_order) || anyDuplicated(leaf_order)) {
    stop("sequences must contain num_tips uniquely named taxa.")
  }
  if (length(coal_times) != num_tips - 1L || any(!is.finite(coal_times)) ||
      any(coal_times <= 0) || is.unsorted(coal_times, strictly = TRUE)) stop("Invalid fixed coal_times.")
  if (init_method == "upgma") {
    initial_tree <- phangorn::upgma(phangorn::dist.hamming(sequences))
    initial_tree <- phylodyn:::update_time(initial_tree, coal_times)
  } else {
    initial_tree <- fmatrix::mytree_from_F(phylodyn:::gen_caterpillar(num_tips), coal_times)
    initial_tree$tip.label <- sample(leaf_order, replace = FALSE)
  }
  A0 <- vb_A_from_tree(initial_tree, leaf_order)
  A_true <- if (!is.null(true_tree)) vb_A_from_tree(true_tree, leaf_order) else NULL
  theta_est <- beta_init * (2 * A0 - 1)
  # Root-row logits and a common shift of the bottom-row logits cancel.
  theta_est[1L, ] <- 0
  theta_est[nrow(theta_est), ] <- theta_est[nrow(theta_est), ] - mean(theta_est[nrow(theta_est), ])
  theta_init <- theta_est
  theta_accumulator <- matrix(step_size, nrow = nrow(A0), ncol = ncol(A0))
  chain_states <- c(list(A0), replicate(num_chains - 1L, vb_A_uniform_draw(leaf_order), simplify = FALSE))

  cat("Building", num_unif_samples, "uniform labeled-ranked A reference draws...\n")
  uniform_values <- t(vapply(seq_len(num_unif_samples), function(i) {
    as.vector(vb_A_uniform_draw(leaf_order))
  }, numeric(length(A0))))
  log_n_trees <- lfactorial(num_tips) + lfactorial(num_tips - 1L) - (num_tips - 1L) * log(2)
  log_prior <- -log_n_trees
  likelihood_cache <- new.env(hash = TRUE, parent = emptyenv())
  cache_size <- 0L
  likelihood <- function(A) {
    key <- phylodyn:::a_matrix_key(A)
    if (exists(key, likelihood_cache, inherits = FALSE)) return(get(key, likelihood_cache, inherits = FALSE))
    tree <- vb_A_to_tree(A, coal_times, leaf_order, check = FALSE)
    ll <- phangorn::pml(tree, sequences, bf = rep(0.25, 4L),
                        Q = c(1, 2, 1, 1, 2, 1), rate = rate)$logLik
    if (!is.finite(ll)) stop("Non-finite sequence likelihood.")
    if (cache_size < max_likelihood_cache) {
      assign(key, ll, likelihood_cache)
      cache_size <<- cache_size + 1L
    }
    ll
  }
  draw_batch <- function(theta, starts) {
    chains <- lapply(starts, function(A) phylodyn:::a_rnni_sample(
      theta, A, n_samples = num_samps, burnin = burnin, thin = take_every
    ))
    log_likelihoods <- lapply(chains, function(chain) vapply(chain$samples, likelihood, numeric(1L)))
    gradients <- lapply(seq_along(chains), function(k) {
      phylodyn:::a_vi_covariance_gradient(chains[[k]]$samples, log_likelihoods[[k]], theta)
    })
    all_samples <- unlist(lapply(chains, `[[`, "samples"), recursive = FALSE)
    values <- do.call(rbind, lapply(all_samples, as.vector))
    ll <- unlist(log_likelihoods, use.names = FALSE)
    mean_A <- matrix(colMeans(values), nrow(theta), ncol(theta), dimnames = dimnames(theta))
    linear <- as.numeric(values %*% as.vector(theta))
    reference_logits <- as.numeric(uniform_values %*% as.vector(theta))
    offset <- max(reference_logits)
    reference_mass <- exp(reference_logits - offset)
    reference_weights <- reference_mass / sum(reference_mass)
    log_Z_mc <- log_n_trees + offset + log(mean(reference_mass))
    entropy_mc <- log_Z_mc - mean(linear)
    ess <- vapply(log_likelihoods, function(v) {
      if (length(unique(v)) < 2L) return(0)
      if (!requireNamespace("coda", quietly = TRUE)) return(NA_real_)
      tryCatch(as.numeric(coda::effectiveSize(v)), error = function(e) NA_real_)
    }, numeric(1L))
    chain_means <- vapply(chains, function(chain) {
      colMeans(do.call(rbind, lapply(chain$samples, as.vector)))
    }, numeric(length(theta)))
    disagreement <- if (num_chains > 1L) {
      max(apply(chain_means, 1L, function(v) max(v) - min(v)))
    } else NA_real_
    finite_ess <- ess[is.finite(ess)]
    list(samples = all_samples, log_likelihood = ll,
         gradient = Reduce("+", gradients) / length(gradients),
         mean_A = mean_A, entropy_mc = entropy_mc,
         elbo_mc = mean(ll) + log_prior + entropy_mc,
         log_Z_mc = log_Z_mc, reference_ess = 1 / sum(reference_weights^2),
         A_variance = sum(mean_A * (1 - mean_A)),
         min_likelihood_ess = if (length(finite_ess)) min(finite_ess) else NA_real_,
         chain_mean_difference = disagreement,
         acceptance = mean(vapply(chains, `[[`, numeric(1L), "acceptance_rate")),
         chains = chains)
  }
  history <- data.frame(iteration = seq_len(num_grad_desc_steps),
                        elbo_mc = NA_real_, mean_log_likelihood = NA_real_,
                        entropy_mc = NA_real_, reference_ess = NA_real_,
                        A_variance = NA_real_, expected_pairwise_A_distance = NA_real_,
                        acceptance = NA_real_, min_likelihood_ess = NA_real_,
                        chain_mean_difference = NA_real_, gradient_norm = NA_real_,
                        step_norm = NA_real_)
  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(num_grad_desc_steps)) {
    batch <- draw_batch(theta_est, chain_states)
    chain_states <- lapply(batch$chains, `[[`, "final")
    gradient <- batch$gradient
    gradient[1L, ] <- 0
    if (any(!is.finite(gradient))) stop("Non-finite Monte Carlo gradient.")
    theta_accumulator <- theta_accumulator + gradient^2
    step <- step_size / sqrt(theta_accumulator) * gradient
    theta_est <- theta_est + step
    theta_est[1L, ] <- 0
    theta_est[nrow(theta_est), ] <- theta_est[nrow(theta_est), ] - mean(theta_est[nrow(theta_est), ])
    history[i, -1L] <- c(batch$elbo_mc, mean(batch$log_likelihood), batch$entropy_mc,
                         batch$reference_ess, batch$A_variance, 2 * batch$A_variance,
                         batch$acceptance, batch$min_likelihood_ess,
                         batch$chain_mean_difference, sqrt(sum(gradient^2)), sqrt(sum(step^2)))
    if (i == 1L || i %% log_every == 0L || i == num_grad_desc_steps) {
      cat(sprintf(paste0("iter %4d | estimated ELBO %.3f | mean log likelihood %.3f | ",
                         "accept %.3f | min ESS %.1f | chain difference %.3f | reference ESS %.1f\n"),
                  i, batch$elbo_mc, mean(batch$log_likelihood), batch$acceptance,
                  batch$min_likelihood_ess, batch$chain_mean_difference, batch$reference_ess))
      if (!is.null(output_dir)) utils::write.csv(history[seq_len(i), ], file.path(output_dir, "trace.csv"), row.names = FALSE)
    }
  }
  # Re-equilibrate and summarize the final parameter, after the last update.
  final <- draw_batch(theta_est, chain_states)
  energy <- vapply(final$samples, function(A) sum(theta_est * A), numeric(1L))
  estimated_A <- final$samples[[which.max(energy)]] # Mode among sampled states.
  estimated_tree <- vb_A_to_tree(estimated_A, coal_times, leaf_order)
  center <- (theta_est + 1) / 2
  center[1L, ] <- 1 # Root-center choice is density-invariant.
  result <- list(theta = theta_est, center_at_beta_1 = center, theta_initial = theta_init,
                 A0 = A0, A_true = A_true, mean_A = final$mean_A,
                 estimated_A = estimated_A, estimated_tree = estimated_tree,
                 samples = final$samples, log_likelihood = final$log_likelihood,
                 history = history, initial_tree = initial_tree, true_tree = true_tree,
                 sequences = sequences, coal_times = coal_times,
                 diagnostics = final[setdiff(names(final), c("samples", "chains", "log_likelihood", "gradient", "mean_A"))],
                 chain_states = lapply(final$chains, `[[`, "final"),
                 config = list(num_tips = num_tips, rate = rate, step_size = step_size,
                               num_samps_per_chain = num_samps, num_chains = num_chains,
                               num_unif_samples = num_unif_samples, take_every = take_every,
                               burnin = burnin, num_grad_desc_steps = num_grad_desc_steps,
                               beta_init = beta_init, init_method = init_method, seed = seed,
                               log_prior = log_prior, likelihood_cache_size = cache_size,
                               elapsed_seconds = proc.time()[["elapsed"]] - start_time,
                               sampler = "rank/NNI Metropolis-Hastings; persistent chains",
                               elbo_method = "uniform-reference log-normalizer estimate; not a certified bound"))
  if (!is.null(output_dir)) {
    saveRDS(result, file.path(output_dir, "result.rds"))
    utils::write.csv(theta_est, file.path(output_dir, "theta.csv"))
    utils::write.csv(center, file.path(output_dir, "center_at_beta_1.csv"))
    ape::write.tree(estimated_tree, file = file.path(output_dir, "sampled_mode.newick"))
  }
  if (make_plots) {
    old_par <- par(mfrow = c(2, 2)); on.exit(par(old_par), add = TRUE)
    plot(history$iteration, history$elbo_mc, type = "l", xlab = "Iteration",
         ylab = "Estimated ELBO", main = "Monte Carlo ELBO diagnostic")
    plot(history$iteration, history$expected_pairwise_A_distance, type = "l",
         xlab = "Iteration", ylab = "Expected differing A entries", main = "Structural spread")
    plot(history$iteration, history$reference_ess, type = "l", xlab = "Iteration",
         ylab = "Reference importance ESS", main = "Normalizer diagnostic")
    plot(estimated_tree, main = "Highest-density sampled genealogy")
  }
  if (final$reference_ess < 10 ||
      (is.finite(final$min_likelihood_ess) && final$min_likelihood_ess < 20) ||
      (!is.na(final$chain_mean_difference) && final$chain_mean_difference > 0.25)) {
    message("Final diagnostics indicate limited precision or chain disagreement; inspect result$diagnostics.")
  }
  cat("Final estimates use", length(final$samples), "MCMC draws; elapsed",
      round(result$config$elapsed_seconds, 2), "seconds.\n")
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  script_dir <- if (!is.null(vb_A_source_file)) dirname(normalizePath(vb_A_source_file)) else getwd()
  output_dir <- if (length(args)) args[1L] else file.path(script_dir, "logs", "VB_master_A", format(Sys.time(), "%Y%m%d_%H%M%S"))
  fit <- run_VB_master_A(output_dir = output_dir)
}
