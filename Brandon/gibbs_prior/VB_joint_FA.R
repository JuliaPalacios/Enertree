# Joint F-geometry / labeled A variational prototype. Source without running.
# Rscript Brandon/gibbs_prior/VB_joint_FA.R runs a 10-tip example in memory.

fa_source_dir <- local({
  file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(file)) {
    arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(arg)) file <- sub("^--file=", "", arg[1])
  }
  if (is.null(file)) file <- "Brandon/gibbs_prior/VB_joint_FA.R"
  dirname(normalizePath(file, mustWork = TRUE))
})
source(file.path(fa_source_dir, "VB_master_A.R"), local = TRUE)
suppressPackageStartupMessages(library(ape)) # phylodyn::update_time uses write.tree()

fa_log_count <- function(n) lfactorial(n) + lfactorial(n - 1L) - (n - 1L) * log(2)
fa_logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
fa_key <- function(A) phylodyn:::a_matrix_key(A)

fa_F <- function(A, coal_times = seq_len(ncol(A) - 1L)) {
  tree <- vb_A_to_tree(A, coal_times, check = FALSE)
  round(phylodyn:::gen_Fmat(tree, tol = 8), 0)
}
fa_cherries <- function(A) sum(rowSums(A) == 2L)
fa_log_multiplicity <- function(A) lfactorial(ncol(A)) - fa_cherries(A) * log(2)
fa_features <- function(A, M, coal_times) {
  F <- fa_F(A, coal_times)
  list(F = F, d2 = sum((F - M)^2))
}

# The offset divides out n!/2^cherries, preserving the F-only Gibbs marginal
# at Theta=0. The state itself remains a single labeled ranked history.
fa_score <- function(A, M, g, theta, coal_times, features = NULL) {
  if (is.null(features)) features <- fa_features(A, M, coal_times)
  -exp(g) * features$d2 + sum(theta * A) - fa_log_multiplicity(A)
}

fa_propose <- function(A, p_label = 0.5) {
  if (runif(1) < p_label) {
    pair <- sample.int(ncol(A), 2L)
    B <- A
    B[, pair] <- A[, rev(pair), drop = FALSE]
    list(A = B, move = "label")
  } else {
    proposal <- phylodyn:::a_rnni_propose(A, check = FALSE)
    list(A = proposal$A, move = proposal$move)
  }
}

fa_sample <- function(initial, M, g, theta, coal_times, n_samples,
                      burnin = 200L, thin = 5L, p_label = 0.5) {
  stopifnot(phylodyn:::a_matrix_valid(initial), identical(dim(theta), dim(initial)),
            n_samples >= 1L, burnin >= 0L, thin >= 1L,
            p_label > 0, p_label < 1)
  current <- initial
  current_score <- fa_score(current, M, g, theta, coal_times)
  draws <- vector("list", n_samples)
  proposals <- accepted <- c(label = 0L, rank = 0L, NNI = 0L)
  cursor <- 0L
  total <- burnin + n_samples * thin
  for (iteration in seq_len(total)) {
    candidate <- fa_propose(current, p_label)
    proposals[candidate$move] <- proposals[candidate$move] + 1L
    candidate_score <- fa_score(candidate$A, M, g, theta, coal_times)
    if (log(runif(1)) < min(0, candidate_score - current_score)) {
      current <- candidate$A
      current_score <- candidate_score
      accepted[candidate$move] <- accepted[candidate$move] + 1L
    }
    if (iteration > burnin && (iteration - burnin) %% thin == 0L) {
      cursor <- cursor + 1L
      draws[[cursor]] <- current
    }
  }
  list(samples = draws, final = current, proposals = proposals, accepted = accepted,
       acceptance = sum(accepted) / total)
}

# Cov_q(h, T) is the ELBO gradient; log Z is absent from h because covariance
# with a constant is zero. M is a continuous lower-triangular F center.
fa_gradient <- function(samples, loglik, M, g, theta, coal_times,
                        log_prior = -fa_log_count(ncol(theta))) {
  stopifnot(length(samples) > 1L, length(loglik) == length(samples),
            all(is.finite(loglik)))
  n <- length(samples)
  F_values <- matrix(0, n, length(M))
  A_values <- matrix(0, n, length(theta))
  d2 <- scores <- numeric(n)
  for (i in seq_len(n)) {
    features <- fa_features(samples[[i]], M, coal_times)
    F_values[i, ] <- as.vector(features$F)
    A_values[i, ] <- as.vector(samples[[i]])
    d2[i] <- features$d2
    scores[i] <- fa_score(samples[[i]], M, g, theta, coal_times,
                          features = features)
  }
  h <- loglik + log_prior - scores
  grad_M <- matrix(2 * exp(g) * as.numeric(cov(F_values, h)), nrow(M))
  grad_g <- -exp(g) * cov(d2, h)
  grad_theta <- matrix(as.numeric(cov(A_values, h)), nrow(theta))
  grad_M[upper.tri(grad_M)] <- 0
  grad_theta[1L, ] <- 0
  list(M = grad_M, g = grad_g, theta = grad_theta,
       mean_loglik = mean(loglik), mean_score = mean(scores))
}

fa_uniform_draw <- function(labels) vb_A_uniform_draw(labels)

fa_logZ_reference <- function(reference, M, g, theta, coal_times) {
  scores <- vapply(reference, fa_score, numeric(1), M = M, g = g,
                   theta = theta, coal_times = coal_times)
  w <- exp(scores - max(scores))
  list(logZ = fa_log_count(ncol(theta)) + fa_logsumexp(scores) - log(length(scores)),
       reference_ess = sum(w)^2 / sum(w^2))
}

# AIS is for an absolute ELBO only. Each intermediate MH kernel preserves
# exp(t * score), starting from an exactly sampled uniform labeled history.
fa_logZ_ais <- function(M, g, theta, coal_times, particles = 128L,
                        temperatures = 100L, moves_per_temperature = 2L) {
  stopifnot(particles >= 2L, temperatures >= 2L,
            moves_per_temperature >= 1L)
  logw <- numeric(particles)
  for (j in seq_len(particles)) {
    A <- fa_uniform_draw(colnames(theta))
    current <- fa_score(A, M, g, theta, coal_times)
    for (k in seq_len(temperatures)) {
      t <- k / temperatures
      logw[j] <- logw[j] + current / temperatures
      for (move in seq_len(moves_per_temperature)) {
        B <- fa_propose(A)$A
        proposed <- fa_score(B, M, g, theta, coal_times)
        if (log(runif(1)) < min(0, t * (proposed - current))) {
          A <- B
          current <- proposed
        }
      }
    }
  }
  w <- exp(logw - max(logw))
  list(logZ = fa_log_count(ncol(theta)) + fa_logsumexp(logw) - log(particles),
       importance_ess = sum(w)^2 / sum(w^2),
       particles = particles, temperatures = temperatures)
}

fa_enumerate <- function(labels) {
  n <- length(labels)
  if (n > 6L) stop("Exact labeled enumeration is limited to n <= 6.")
  blocks <- lapply(seq_len(n), function(j) as.integer(seq_len(n) == j))
  out <- list()
  visit <- function(active, A, row) {
    if (row == 0L) { out[[length(out) + 1L]] <<- A; return() }
    for (pair in combn(seq_along(active), 2L, simplify = FALSE)) {
      B <- A
      B[row, ] <- active[[pair[1]]] + active[[pair[2]]]
      visit(c(active[-pair], list(B[row, ])), B, row - 1L)
    }
  }
  visit(blocks, matrix(0L, n - 1L, n, dimnames = list(NULL, labels)), n - 1L)
  out
}

fa_likelihood <- function(A, sequences, coal_times, rate) {
  tree <- vb_A_to_tree(A, coal_times, names(sequences), check = FALSE)
  as.numeric(phangorn::pml(tree, sequences, bf = rep(0.25, 4L),
                           Q = c(1, 2, 1, 1, 2, 1), rate = rate)$logLik)
}

fa_fit <- function(sequences, coal_times, rate = 0.01, steps = 10L,
                   samples_per_chain = 40L, chains = 2L, burnin = 300L,
                   thin = 10L, step_size = 0.1, seed = 42L,
                   reference_n = 2000L, init_beta = 0.25,
                   init_theta_scale = 0.3, ais_particles = 128L) {
  set.seed(seed)
  n <- length(sequences)
  stopifnot(n >= 3L, length(coal_times) == n - 1L,
            !is.unsorted(coal_times, strictly = TRUE),
            length(names(sequences)) == n, !anyDuplicated(names(sequences)))
  upgma <- phangorn::upgma(phangorn::dist.hamming(sequences))
  initial <- vb_A_from_tree(phylodyn:::update_time(upgma, coal_times), names(sequences))
  M <- fa_F(initial, coal_times)
  g <- log(init_beta)
  theta <- init_theta_scale * (2 * initial - 1)
  theta[1L, ] <- 0
  M_acc <- matrix(1, nrow(M), ncol(M))
  g_acc <- 1
  theta_acc <- matrix(1, nrow(theta), ncol(theta))
  states <- c(list(initial), replicate(chains - 1L,
                                      fa_uniform_draw(names(sequences)), simplify = FALSE))
  likelihood_cache <- new.env(hash = TRUE, parent = emptyenv())
  likelihood <- function(A) {
    key <- fa_key(A)
    if (!exists(key, likelihood_cache, inherits = FALSE))
      assign(key, fa_likelihood(A, sequences, coal_times, rate), likelihood_cache)
    get(key, likelihood_cache, inherits = FALSE)
  }
  trace <- vector("list", steps)
  for (iteration in seq_len(steps)) {
    batches <- lapply(states, fa_sample, M = M, g = g, theta = theta,
                      coal_times = coal_times, n_samples = samples_per_chain,
                      burnin = burnin, thin = thin)
    states <- lapply(batches, `[[`, "final")
    draws <- unlist(lapply(batches, `[[`, "samples"), recursive = FALSE)
    ll <- vapply(draws, likelihood, numeric(1))
    grad <- fa_gradient(draws, ll, M, g, theta, coal_times)
    M_acc <- M_acc + grad$M^2
    g_acc <- g_acc + grad$g^2
    theta_acc <- theta_acc + grad$theta^2
    M <- M + step_size * grad$M / sqrt(M_acc)
    M <- matrix(pmax(0, pmin(n, M)), n - 1L, n - 1L)
    M[upper.tri(M)] <- 0
    g <- max(log(0.005), min(log(3), g + step_size * grad$g / sqrt(g_acc)))
    theta <- theta + step_size * grad$theta / sqrt(theta_acc)
    theta[1L, ] <- 0
    theta[nrow(theta), ] <- theta[nrow(theta), ] - mean(theta[nrow(theta), ])
    trace[[iteration]] <- data.frame(iteration = iteration,
                                    mean_loglik = mean(ll),
                                    acceptance = mean(vapply(batches, `[[`, numeric(1), "acceptance")),
                                    beta = exp(g))
  }
  final <- lapply(states, fa_sample, M = M, g = g, theta = theta,
                  coal_times = coal_times, n_samples = samples_per_chain,
                  burnin = burnin, thin = thin)
  reference <- replicate(reference_n, fa_uniform_draw(names(sequences)), simplify = FALSE)
  normalization <- fa_logZ_reference(reference, M, g, theta, coal_times)
  ais <- fa_logZ_ais(M, g, theta, coal_times, particles = ais_particles)
  list(M = M, g = g, theta = theta, chains = final, trace = do.call(rbind, trace),
       normalization = normalization, ais = ais, coal_times = coal_times, rate = rate,
       labels = names(sequences), seed = seed, likelihood_cache_size = length(ls(likelihood_cache)))
}

if (sys.nframe() == 0L) {
  data_dir <- file.path(fa_source_dir, "A matrix analysis", "outputs", "20261024_184601")
  previous <- readRDS(file.path(fa_source_dir, "logs", "VB_master_A", "20261024_184601", "result.rds"))
  sequences <- phangorn::as.phyDat(ape::read.FASTA(file.path(data_dir, "sequences.fasta")))
  stopifnot(identical(sequences, previous$sequences))
  fit <- fa_fit(sequences, previous$coal_times, rate = previous$config$rate,
                steps = 50L, samples_per_chain = 100L, chains = 2L,
                burnin = 400L, thin = 10L, reference_n = 3000L)
  print(tail(fit$trace, 1L), row.names = FALSE)
  cat("Normalizer diagnostics (not an ELBO):\n")
  print(list(reference = fit$normalization, ais = fit$ais))
}
