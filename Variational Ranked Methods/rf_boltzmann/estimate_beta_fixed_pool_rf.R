suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

# Fixed-pool beta estimation for the rooted RF Boltzmann family.
#
# Variational family on a fixed empirical support:
#   q_beta(T_i) is proportional to exp(-beta * rooted_RF(T_i, M)).
#
# The pool is generated in several complementary ways:
#   1. M itself
#   2. direct one-NNI proposals from M
#   3. random rooted trees from ape::rtree()
#   4. low-beta rooted RF NNI chains around M
#
# After the pool is built, beta estimation uses deterministic reweighting of
# the same scored trees. This avoids the failure mode where a beta-dependent
# NNI sampler freezes at M.

set.seed(42)

true_tree_file <- "rf_true_tree_10tip.newick"
sequence_file <- "rf_sequences_10tip.fasta"

mutation_rate <- 0.01

n_random_pool <- 50000
n_direct_nni_pool <- 5000
n_local_pool_per_beta <- 25000
local_pool_betas <- c(0.1, 0.25, 0.5, 1)
local_pool_burnin <- 1000
local_pool_thin <- 1

nni_move_options <- c(1, 3, 4)
nni_move_probabilities <- c(0.5, 0.3, 0.2)

beta_search_min <- 1e-4
beta_search_max <- 50
beta_grid_size <- 300

final_sample_size <- 100000

output_pool_tree_file <- "fixed_pool_rooted_trees.newick"
output_pool_score_file <- "fixed_pool_rooted_tree_scores.csv"
output_topology_count_file <- "fixed_pool_rooted_topology_counts.csv"
output_grid_file <- "fixed_pool_beta_grid.csv"
output_fit_file <- "fixed_pool_beta_fit.rds"
output_sample_tree_file <- "fixed_pool_estimated_beta_100000_trees.newick"

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
    stop("The true tree must be rooted.")
  }

  list(center_tree = center_tree, sequences = sequences)
}

set_fixed_branch_lengths <- function(tree, branch_length) {
  tree$edge.length <- rep(branch_length, nrow(tree$edge))
  tree
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

propose_nni <- function(tree, n_moves, branch_length) {
  proposal <- rNNI(tree, moves = n_moves)

  if (!is.rooted(proposal)) {
    stop("The NNI proposal returned an unrooted tree.")
  }

  set_fixed_branch_lengths(proposal, branch_length)
}

accept_proposal <- function(current_energy, proposed_energy, beta) {
  log_acceptance_ratio <- -beta * (proposed_energy - current_energy)
  log(runif(1)) < min(0, log_acceptance_ratio)
}

sample_low_beta_pool <- function(center_tree,
                                 beta,
                                 n_samples,
                                 burnin,
                                 thin,
                                 branch_length,
                                 nni_move_options,
                                 nni_move_probabilities) {
  current_tree <- set_fixed_branch_lengths(center_tree, branch_length)
  current_energy <- rf_energy(current_tree, center_tree)

  total_steps <- burnin + n_samples * thin
  samples <- vector("list", n_samples)

  accepted <- 0
  saved <- 0

  for (step in seq_len(total_steps)) {
    proposed_move_count <- choose_nni_move_count(
      nni_move_options,
      nni_move_probabilities
    )
    proposed_tree <- propose_nni(
      current_tree,
      n_moves = proposed_move_count,
      branch_length = branch_length
    )
    proposed_energy <- rf_energy(proposed_tree, center_tree)

    if (accept_proposal(current_energy, proposed_energy, beta)) {
      current_tree <- proposed_tree
      current_energy <- proposed_energy
      accepted <- accepted + 1
    }

    if (step > burnin && (step - burnin) %% thin == 0) {
      saved <- saved + 1
      samples[[saved]] <- current_tree
    }
  }

  list(
    trees = samples,
    acceptance_rate = accepted / total_steps
  )
}

generate_random_rooted_pool <- function(n_trees, center_tree, branch_length) {
  n_tips <- length(center_tree$tip.label)
  trees <- vector("list", n_trees)

  for (i in seq_len(n_trees)) {
    tree <- rtree(
      n = n_tips,
      rooted = TRUE,
      tip.label = center_tree$tip.label,
      br = NULL,
      equiprob = TRUE
    )
    trees[[i]] <- set_fixed_branch_lengths(tree, branch_length)

    if (i %% 5000 == 0) {
      cat("Generated", i, "random rooted trees\n")
    }
  }

  trees
}

generate_direct_nni_pool <- function(n_trees, center_tree, branch_length) {
  trees <- vector("list", n_trees)
  center_tree <- set_fixed_branch_lengths(center_tree, branch_length)

  for (i in seq_len(n_trees)) {
    trees[[i]] <- propose_nni(
      center_tree,
      n_moves = 1,
      branch_length = branch_length
    )

    if (i %% 1000 == 0) {
      cat("Generated", i, "direct one-NNI trees from M\n")
    }
  }

  trees
}

rooted_topology_key <- function(tree) {
  tree <- reorder(tree, "postorder")
  tip_labels <- tree$tip.label
  n_tips <- length(tip_labels)
  all_tips <- sort(tip_labels)
  children_by_parent <- split(tree$edge[, 2], tree$edge[, 1])

  descendants <- function(node) {
    if (node <= n_tips) {
      return(tip_labels[node])
    }

    children <- children_by_parent[[as.character(node)]]
    sort(unlist(lapply(children, descendants), use.names = FALSE))
  }

  internal_nodes <- sort(unique(tree$edge[, 1]))
  clades <- character(0)

  for (node in internal_nodes) {
    clade <- descendants(node)
    if (length(clade) > 1 && length(clade) < n_tips) {
      clades <- c(clades, paste(clade, collapse = ","))
    }
  }

  paste(sort(clades), collapse = "|")
}

summarize_topology_counts <- function(trees, sources) {
  keys <- vapply(trees, rooted_topology_key, character(1))
  count_table <- sort(table(keys), decreasing = TRUE)
  source_table <- split(sources, keys)

  data.frame(
    topology_key = names(count_table),
    count = as.integer(count_table),
    sources = vapply(
      names(count_table),
      function(key) {
        source_counts <- sort(table(source_table[[key]]), decreasing = TRUE)
        paste(paste(names(source_counts), source_counts, sep = ":"), collapse = ";")
      },
      character(1)
    )
  )
}

log_likelihood_tree <- function(tree, sequences, mutation_rate) {
  if (is.null(tree$edge.length) || length(tree$edge.length) != nrow(tree$edge)) {
    stop("Every tree scored by the likelihood must have one length per edge.")
  }

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

score_pool <- function(trees, sources, center_tree, sequences, mutation_rate) {
  n_trees <- length(trees)
  scores <- data.frame(
    tree_index = seq_len(n_trees),
    source = sources,
    rooted_rf = numeric(n_trees),
    log_likelihood = numeric(n_trees)
  )

  for (i in seq_len(n_trees)) {
    scores$rooted_rf[i] <- rf_energy(trees[[i]], center_tree)
    scores$log_likelihood[i] <- log_likelihood_tree(
      trees[[i]],
      sequences = sequences,
      mutation_rate = mutation_rate
    )

    if (i %% 1000 == 0) {
      cat("Scored", i, "of", n_trees, "pool trees\n")
    }
  }

  scores
}

logsumexp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

softmax <- function(log_weights) {
  exp(log_weights - logsumexp(log_weights))
}

pool_summary_at_beta <- function(beta, rooted_rf, log_likelihood) {
  log_weights <- -beta * rooted_rf
  weights <- softmax(log_weights)

  mean_rf <- sum(weights * rooted_rf)
  mean_log_likelihood <- sum(weights * log_likelihood)
  entropy <- -sum(weights * log(weights))
  elbo <- mean_log_likelihood + entropy

  score <- beta * (mean_rf - rooted_rf)
  elbo_part_without_constant <- log_likelihood + beta * rooted_rf
  gradient_gamma <- sum(weights * score * elbo_part_without_constant)

  ess <- 1 / sum(weights^2)

  data.frame(
    beta = beta,
    gamma = log(beta),
    elbo = elbo,
    gradient_gamma = gradient_gamma,
    mean_rf = mean_rf,
    sd_rf = sqrt(sum(weights * (rooted_rf - mean_rf)^2)),
    mean_log_likelihood = mean_log_likelihood,
    ess = ess,
    max_weight = max(weights)
  )
}

scan_beta_grid <- function(rooted_rf, log_likelihood, beta_min, beta_max, grid_size) {
  beta_grid <- exp(seq(log(beta_min), log(beta_max), length.out = grid_size))
  rows <- lapply(
    beta_grid,
    pool_summary_at_beta,
    rooted_rf = rooted_rf,
    log_likelihood = log_likelihood
  )

  do.call(rbind, rows)
}

fit_beta_by_pool_elbo <- function(rooted_rf, log_likelihood, beta_min, beta_max) {
  objective <- function(gamma) {
    beta <- exp(gamma)
    -pool_summary_at_beta(beta, rooted_rf, log_likelihood)$elbo
  }

  opt <- optimize(
    f = objective,
    interval = c(log(beta_min), log(beta_max))
  )

  beta_hat <- exp(opt$minimum)
  summary <- pool_summary_at_beta(beta_hat, rooted_rf, log_likelihood)

  list(
    beta = beta_hat,
    gamma = opt$minimum,
    summary = summary,
    objective_value = -opt$objective,
    hit_lower_boundary = beta_hat <= beta_min * 1.001,
    hit_upper_boundary = beta_hat >= beta_max / 1.001
  )
}

write_tree_sample <- function(trees, file) {
  class(trees) <- "multiPhylo"
  write.tree(trees, file = file)
}

inputs <- read_inputs(true_tree_file, sequence_file)
center_tree <- inputs$center_tree
sequences <- inputs$sequences

fixed_branch_length <- mean(center_tree$edge.length)
center_tree_fixed <- set_fixed_branch_lengths(center_tree, fixed_branch_length)

cat("Using fixed branch length", fixed_branch_length, "for every edge in every pool tree\n")

random_trees <- generate_random_rooted_pool(
  n_trees = n_random_pool,
  center_tree = center_tree_fixed,
  branch_length = fixed_branch_length
)

direct_nni_trees <- generate_direct_nni_pool(
  n_trees = n_direct_nni_pool,
  center_tree = center_tree_fixed,
  branch_length = fixed_branch_length
)

local_trees <- list()
local_sources <- character(0)
local_acceptance_rates <- numeric(length(local_pool_betas))

for (i in seq_along(local_pool_betas)) {
  beta <- local_pool_betas[i]
  source_name <- paste0("local_beta_", beta)

  local_pool <- sample_low_beta_pool(
    center_tree = center_tree_fixed,
    beta = beta,
    n_samples = n_local_pool_per_beta,
    burnin = local_pool_burnin,
    thin = local_pool_thin,
    branch_length = fixed_branch_length,
    nni_move_options = nni_move_options,
    nni_move_probabilities = nni_move_probabilities
  )

  local_trees <- c(local_trees, local_pool$trees)
  local_sources <- c(local_sources, rep(source_name, length(local_pool$trees)))
  local_acceptance_rates[i] <- local_pool$acceptance_rate

  cat("Local pool beta", beta, "acceptance rate:", local_pool$acceptance_rate, "\n")
}

all_trees <- c(
  list(center_tree_fixed),
  direct_nni_trees,
  random_trees,
  local_trees
)
all_sources <- c(
  "center_tree",
  rep("one_nni_from_center", length(direct_nni_trees)),
  rep("random_rooted", length(random_trees)),
  local_sources
)

topology_counts <- summarize_topology_counts(all_trees, all_sources)

cat("Pool size including duplicate proposal draws:", length(all_trees), "\n")
cat("Unique rooted topologies in pool:", nrow(topology_counts), "\n")

pool_scores <- score_pool(
  trees = all_trees,
  sources = all_sources,
  center_tree = center_tree_fixed,
  sequences = sequences,
  mutation_rate = mutation_rate
)

write_tree_sample(all_trees, output_pool_tree_file)
write.csv(pool_scores, output_pool_score_file, row.names = FALSE)
write.csv(topology_counts, output_topology_count_file, row.names = FALSE)

beta_grid <- scan_beta_grid(
  rooted_rf = pool_scores$rooted_rf,
  log_likelihood = pool_scores$log_likelihood,
  beta_min = beta_search_min,
  beta_max = beta_search_max,
  grid_size = beta_grid_size
)

write.csv(beta_grid, output_grid_file, row.names = FALSE)

fit <- fit_beta_by_pool_elbo(
  rooted_rf = pool_scores$rooted_rf,
  log_likelihood = pool_scores$log_likelihood,
  beta_min = beta_search_min,
  beta_max = beta_search_max
)

weights_at_fit <- softmax(-fit$beta * pool_scores$rooted_rf)
sample_indices <- sample.int(
  n = length(all_trees),
  size = final_sample_size,
  replace = TRUE,
  prob = weights_at_fit
)
sampled_trees <- all_trees[sample_indices]
write_tree_sample(sampled_trees, output_sample_tree_file)

fit$settings <- list(
  true_tree_file = true_tree_file,
  sequence_file = sequence_file,
  mutation_rate = mutation_rate,
  n_random_pool = n_random_pool,
  n_direct_nni_pool = n_direct_nni_pool,
  n_local_pool_per_beta = n_local_pool_per_beta,
  local_pool_betas = local_pool_betas,
  local_pool_acceptance_rates = local_acceptance_rates,
  fixed_branch_length = fixed_branch_length,
  branch_length_mode = "all pool trees use the mean true-tree edge length",
  beta_search_min = beta_search_min,
  beta_search_max = beta_search_max,
  beta_grid_size = beta_grid_size,
  nni_move_options = nni_move_options,
  nni_move_probabilities = nni_move_probabilities,
  pool_size_including_duplicate_proposal_draws = length(all_trees),
  unique_rooted_topologies = nrow(topology_counts),
  fitting_pool_uses_duplicate_proposal_draws = TRUE
)

saveRDS(fit, output_fit_file)

cat("\nFixed-pool beta estimate:", fit$beta, "\n")
cat("ELBO at beta:", fit$summary$elbo, "\n")
cat("Gradient at beta:", fit$summary$gradient_gamma, "\n")
cat("ESS at beta:", fit$summary$ess, "of", nrow(pool_scores), "\n")

if (fit$hit_lower_boundary || fit$hit_upper_boundary) {
  warning("The fitted beta is on the search boundary; widen beta_search_min/beta_search_max.")
}

cat("Wrote pool trees to", output_pool_tree_file, "\n")
cat("Wrote pool scores to", output_pool_score_file, "\n")
cat("Wrote topology counts to", output_topology_count_file, "\n")
cat("Wrote beta grid to", output_grid_file, "\n")
cat("Wrote fit to", output_fit_file, "\n")
cat("Wrote sampled fitted-pool trees to", output_sample_tree_file, "\n")
