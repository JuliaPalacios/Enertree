suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

# Estimate the rooted RF Boltzmann temperature beta when the center tree M is known.
#
# Variational family:
#   q_beta(T) is proportional to exp(-beta * rooted_RF(T, M)).
#
# We use unsquared RF distance here. Squaring is natural for Euclidean
# F-matrix distances, but RF is already a discrete tree-edit distance.
#
# The sampler below is Metropolis-Hastings with random NNI proposals.
# phangorn already provides the NNI proposal through rNNI(), so we wrap it
# instead of writing our own tree surgery.

set.seed(42)

true_tree_file <- "rf_true_tree_10tip.newick"
sequence_file <- "rf_sequences_10tip.fasta"

mutation_rate <- 0.01

initial_beta <- 1
n_iterations <- 1000
n_samples <- 10000
burnin <- 1000
thin <- 1

learning_rate <- 0.001
max_gamma_step <- 0.25
min_acceptance_rate_for_positive_step <- 0.02
min_rf_sd_for_positive_step <- 0.10

nni_move_options <- c(1, 3, 4)
nni_move_probabilities <- c(0.5, 0.3, 0.2)

output_history_file <- "rooted_rf_beta_history.csv"
output_fit_file <- "rooted_rf_beta_fit.rds"
output_sample_tree_file <- "rooted_rf_boltzmann_estimated_beta_100000_thinned_trees.newick"

final_sample_size <- 100000
final_sample_burnin <- 1000
final_sample_thin <- 10

read_inputs <- function(tree_file, fasta_file) {
  center_tree <- read.tree(tree_file)
  sequences <- read.phyDat(fasta_file, format = "fasta", type = "DNA")

  if (!setequal(center_tree$tip.label, names(sequences))) {
    stop("The tree tips and sequence names do not match.")
  }

  if (is.null(center_tree$edge.length)) {
    stop("The true tree does not have branch lengths.")
  }

  if (!is.rooted(center_tree)) {
    stop("The true tree must be rooted for rooted RF beta estimation.")
  }

  list(center_tree = center_tree, sequences = sequences)
}

rf_energy <- function(tree, center_tree) {
  as.numeric(RF.dist(
    tree,
    center_tree,
    normalize = FALSE,
    check.labels = TRUE,
    rooted = TRUE
  ))
}

choose_nni_move_count <- function(move_options, move_probabilities) {
  sample(move_options, size = 1, prob = move_probabilities)
}

propose_nni <- function(tree, n_moves) {
  proposal <- rNNI(tree, moves = n_moves)

  if (!is.rooted(proposal)) {
    stop("The NNI proposal returned an unrooted tree.")
  }

  proposal
}

accept_proposal <- function(current_energy, proposed_energy, beta) {
  log_acceptance_ratio <- -beta * (proposed_energy - current_energy)
  log(runif(1)) < min(0, log_acceptance_ratio)
}

sample_rf_boltzmann <- function(beta,
                                center_tree,
                                n_samples,
                                burnin,
                                thin,
                                nni_move_options,
                                nni_move_probabilities) {
  current_tree <- center_tree
  current_energy <- rf_energy(current_tree, center_tree)

  total_steps <- burnin + n_samples * thin
  samples <- vector("list", n_samples)
  energies <- numeric(n_samples)

  accepted <- 0
  saved <- 0

  for (step in seq_len(total_steps)) {
    proposed_move_count <- choose_nni_move_count(
      nni_move_options,
      nni_move_probabilities
    )
    proposed_tree <- propose_nni(current_tree, proposed_move_count)
    proposed_energy <- rf_energy(proposed_tree, center_tree)

    if (accept_proposal(current_energy, proposed_energy, beta)) {
      current_tree <- proposed_tree
      current_energy <- proposed_energy
      accepted <- accepted + 1
    }

    if (step > burnin && (step - burnin) %% thin == 0) {
      saved <- saved + 1
      samples[[saved]] <- current_tree
      energies[saved] <- current_energy
    }
  }

  list(
    trees = samples,
    energies = energies,
    acceptance_rate = accepted / total_steps
  )
}

log_likelihood_tree <- function(tree, sequences, mutation_rate) {
  if (is.null(tree$edge.length) || length(tree$edge.length) != nrow(tree$edge)) {
    stop("Every tree scored by the likelihood must have one length per edge.")
  }

  # Match the data-generating model: HKY with equal base frequencies and
  # transition/transversion ratio 2. The mutation rate is passed through rate.
  fit <- pml(
    tree,
    sequences,
    bf = rep(0.25, 4),
    Q = c(1, 2, 1, 1, 2, 1),
    rate = mutation_rate,
    model = "USER"
  )

  as.numeric(logLik(fit))
}

estimate_gradient <- function(energies, log_likelihoods, beta) {
  mean_energy <- mean(energies)
  score <- beta * (mean_energy - energies)

  # The unknown log normalizing constant drops out because mean(score) = 0.
  elbo_part_without_constant <- log_likelihoods + beta * energies

  mean(score * elbo_part_without_constant)
}

estimate_beta <- function(center_tree,
                          sequences,
                          initial_beta,
                          n_iterations,
                          n_samples,
                          burnin,
                          thin,
                          mutation_rate,
                          learning_rate,
                          max_gamma_step,
                          min_acceptance_rate_for_positive_step,
                          min_rf_sd_for_positive_step,
                          nni_move_options,
                          nni_move_probabilities) {
  gamma <- log(initial_beta)
  history <- data.frame()

  for (iteration in seq_len(n_iterations)) {
    beta <- exp(gamma)

    cat(sprintf(
      "\nIteration %d/%d | beta = %.6f\n",
      iteration, n_iterations, beta
    ))

    sampled <- sample_rf_boltzmann(
      beta = beta,
      center_tree = center_tree,
      n_samples = n_samples,
      burnin = burnin,
      thin = thin,
      nni_move_options = nni_move_options,
      nni_move_probabilities = nni_move_probabilities
    )

    log_likelihoods <- vapply(
      sampled$trees,
      log_likelihood_tree,
      numeric(1),
      sequences = sequences,
      mutation_rate = mutation_rate
    )

    gradient_gamma <- estimate_gradient(
      energies = sampled$energies,
      log_likelihoods = log_likelihoods,
      beta = beta
    )

    raw_gamma_step <- learning_rate * gradient_gamma
    gamma_step <- max(-max_gamma_step, min(max_gamma_step, raw_gamma_step))

    sampler_is_too_sticky <- sampled$acceptance_rate < min_acceptance_rate_for_positive_step ||
      sd(sampled$energies) < min_rf_sd_for_positive_step

    if (sampler_is_too_sticky && gamma_step > 0) {
      gamma_step <- 0
    }

    gamma <- gamma + gamma_step

    mean_energy <- mean(sampled$energies)
    var_energy <- mean((sampled$energies - mean_energy)^2)
    cov_energy_loglik <- mean(
      (sampled$energies - mean_energy) *
        (log_likelihoods - mean(log_likelihoods))
    )
    beta_covariance_estimate <- if (var_energy > 0) {
      -cov_energy_loglik / var_energy
    } else {
      NA_real_
    }

    row <- data.frame(
      iteration = iteration,
      beta = beta,
      gamma = log(beta),
      gradient_gamma = gradient_gamma,
      raw_gamma_step = raw_gamma_step,
      gamma_step = gamma_step,
      sampler_is_too_sticky = sampler_is_too_sticky,
      acceptance_rate = sampled$acceptance_rate,
      mean_rf = mean_energy,
      sd_rf = sd(sampled$energies),
      mean_log_likelihood = mean(log_likelihoods),
      cov_rf_loglik = cov_energy_loglik,
      var_rf = var_energy,
      beta_covariance_estimate = beta_covariance_estimate
    )

    history <- rbind(history, row)
    write.csv(history, output_history_file, row.names = FALSE)

    cat(sprintf(
      "  accept = %.3f | mean RF = %.3f | mean logLik = %.2f | grad = %.3f | step = %.3f | sticky = %s\n",
      row$acceptance_rate,
      row$mean_rf,
      row$mean_log_likelihood,
      row$gradient_gamma,
      row$gamma_step,
      row$sampler_is_too_sticky
    ))

    if (!is.na(row$beta_covariance_estimate)) {
      cat(sprintf(
        "  covariance beta estimate from this batch = %.6f\n",
        row$beta_covariance_estimate
      ))
    }
  }

  list(
    beta = exp(gamma),
    gamma = gamma,
    history = history,
    settings = list(
      initial_beta = initial_beta,
      n_iterations = n_iterations,
      n_samples = n_samples,
      burnin = burnin,
      thin = thin,
      mutation_rate = mutation_rate,
      branch_length_mode = "true tree lengths inherited under NNI proposals; mutation rate passed through pml(rate)",
      learning_rate = learning_rate,
      max_gamma_step = max_gamma_step,
      min_acceptance_rate_for_positive_step = min_acceptance_rate_for_positive_step,
      min_rf_sd_for_positive_step = min_rf_sd_for_positive_step,
      energy = "unsquared rooted RF distance",
      proposal = "phangorn::rNNI with mixture of NNI move counts",
      nni_move_options = nni_move_options,
      nni_move_probabilities = nni_move_probabilities
    )
  )
}

write_tree_sample <- function(trees, file) {
  class(trees) <- "multiPhylo"
  write.tree(trees, file = file)
}

inputs <- read_inputs(true_tree_file, sequence_file)

fit <- estimate_beta(
  center_tree = inputs$center_tree,
  sequences = inputs$sequences,
  initial_beta = initial_beta,
  n_iterations = n_iterations,
  n_samples = n_samples,
  burnin = burnin,
  thin = thin,
  mutation_rate = mutation_rate,
  learning_rate = learning_rate,
  max_gamma_step = max_gamma_step,
  min_acceptance_rate_for_positive_step = min_acceptance_rate_for_positive_step,
  min_rf_sd_for_positive_step = min_rf_sd_for_positive_step,
  nni_move_options = nni_move_options,
  nni_move_probabilities = nni_move_probabilities
)

saveRDS(fit, output_fit_file)

final_sample <- sample_rf_boltzmann(
  beta = fit$beta,
  center_tree = inputs$center_tree,
  n_samples = final_sample_size,
  burnin = final_sample_burnin,
  thin = final_sample_thin,
  nni_move_options = nni_move_options,
  nni_move_probabilities = nni_move_probabilities
)

write_tree_sample(final_sample$trees, output_sample_tree_file)

cat(sprintf(
  "\nFinal beta estimate: %.6f\nHistory written to %s\nFit written to %s\nSampled trees written to %s\n",
  fit$beta,
  output_history_file,
  output_fit_file,
  output_sample_tree_file
))
