setwd("/Users/brandonmarks/Desktop/Research/Palacios Lab/Untitled/Enertree/Brandon/gibbs_prior")

library(ape)
library(phylodyn)
library(phangorn)
library(phyclust)
library(phylotools)
library(fmatrix)
library(parallel)
source("utils.R")

log_posterior_joint_timed <- function(
  beta,
  M,
  coal_times,
  prior_sd,
  lambda = 1,
  start_tree = NULL,
  mode = "average",
  R = 10
) {
  timings <- list()

  t0 <- proc.time()[["elapsed"]]
  tree_fmat <- sampleF(M, beta, iter = 2, diam = 1, startF = start_tree, take_every = 1)$chainF[[2]]
  timings$sampleF <- proc.time()[["elapsed"]] - t0

  t1 <- proc.time()[["elapsed"]]
  log_botlzman_dist <- Logdistance_prob(tree_fmat, M, beta = beta, diam = 1)
  timings$logdistance <- proc.time()[["elapsed"]] - t1

  t2 <- proc.time()[["elapsed"]]
  if (mode == "average") {
    log_likelihood <- 0
  } else if (mode == "max") {
    log_likelihood <- -Inf
  } else {
    stop("Invalid mode. Choose 'average' or 'max'.")
  }

  pml_times <- numeric(R)
  for (i in 1:R) {
    ti <- proc.time()[["elapsed"]]
    rooted_tree <- mytree_from_F(tree_fmat, coal_times)
    ll <- pml(rooted_tree, sequences)$log
    pml_times[i] <- proc.time()[["elapsed"]] - ti
    if (mode == "average") {
      log_likelihood <- log_likelihood + ll
    } else {
      log_likelihood <- max(log_likelihood, ll)
    }
  }
  if (mode == "average") {
    log_likelihood <- log_likelihood / R
  }
  timings$tree_and_likelihood <- proc.time()[["elapsed"]] - t2
  timings$mytree_from_F_avg <- mean(pml_times)
  timings$pml_avg <- mean(pml_times)
  timings$pml_total <- sum(pml_times)

  t3 <- proc.time()[["elapsed"]]
  n_tips <- ncol(M_true_fmat) + 1
  log_rooted_given_unrooted <- log(2^(num_cherries(tree_fmat)) / factorial(n_tips))
  log_prior <- dnorm(log(beta), mean = 0, sd = prior_sd, log = TRUE)
  log_posterior <- log_likelihood / lambda + log_rooted_given_unrooted + log_botlzman_dist + log_prior
  timings$bookkeeping <- proc.time()[["elapsed"]] - t3

  list(log_posterior = log_posterior, tree_fmat = tree_fmat, timings = timings)
}

log_posterior_joint_timed_parallel <- function(
  beta,
  M,
  coal_times,
  prior_sd,
  lambda = 1,
  start_tree = NULL,
  mode = "average",
  R = 10,
  cores = 2
) {
  timings <- list()

  t0 <- proc.time()[["elapsed"]]
  tree_fmat <- sampleF(M, beta, iter = 2, diam = 1, startF = start_tree, take_every = 1)$chainF[[2]]
  timings$sampleF <- proc.time()[["elapsed"]] - t0

  t1 <- proc.time()[["elapsed"]]
  log_botlzman_dist <- Logdistance_prob(tree_fmat, M, beta = beta, diam = 1)
  timings$logdistance <- proc.time()[["elapsed"]] - t1

  t2 <- proc.time()[["elapsed"]]
  rooted_tree <- mytree_from_F(tree_fmat, coal_times)
  ll_vec <- unlist(parallel::mclapply(
    seq_len(R),
    function(i) pml(rooted_tree, sequences)$log,
    mc.cores = cores,
    mc.preschedule = FALSE
  ))
  log_likelihood <- if (mode == "average") mean(ll_vec) else max(ll_vec)
  timings$tree_and_likelihood <- proc.time()[["elapsed"]] - t2
  timings$pml_total <- NA_real_
  timings$pml_avg <- NA_real_
  timings$mytree_from_F_avg <- NA_real_

  t3 <- proc.time()[["elapsed"]]
  n_tips <- ncol(M_true_fmat) + 1
  log_rooted_given_unrooted <- log(2^(num_cherries(tree_fmat)) / factorial(n_tips))
  log_prior <- dnorm(log(beta), mean = 0, sd = prior_sd, log = TRUE)
  log_posterior <- log_likelihood / lambda + log_rooted_given_unrooted + log_botlzman_dist + log_prior
  timings$bookkeeping <- proc.time()[["elapsed"]] - t3

  list(log_posterior = log_posterior, tree_fmat = tree_fmat, timings = timings)
}

run_timed_joint_chain <- function(
  M_init,
  beta_init,
  prior_sd,
  n_iter,
  coal_times,
  proposal_sd,
  lambda = 1,
  mode = "average",
  R = 10
) {
  chain_timings <- vector("list", n_iter)
  current_beta <- beta_init
  current_g <- log(current_beta)
  current_M <- M_init

  t_iter <- proc.time()[["elapsed"]]
  current_res <- log_posterior_joint_timed(
    current_beta, current_M, coal_times, prior_sd,
    lambda = lambda, start_tree = current_M, mode = mode, R = R
  )
  current_log_post <- current_res$log_posterior
  current_tree <- current_res$tree_fmat
  chain_timings[[1]] <- current_res$timings
  chain_timings[[1]]$iteration_total <- proc.time()[["elapsed"]] - t_iter
  chain_timings[[1]]$iteration_overhead <- chain_timings[[1]]$iteration_total - chain_timings[[1]]$tree_and_likelihood - chain_timings[[1]]$sampleF - chain_timings[[1]]$logdistance - chain_timings[[1]]$bookkeeping

  for (iter in 2:n_iter) {
    t_iter <- proc.time()[["elapsed"]]
    proposed_g <- rnorm(1, mean = current_g, sd = proposal_sd)
    proposed_beta <- exp(proposed_g)

    beta_res <- log_posterior_joint_timed(
      proposed_beta, current_M, coal_times, prior_sd,
      lambda = lambda, start_tree = current_tree, mode = mode, R = R
    )
    proposed_log_post <- beta_res$log_posterior

    if (log(runif(1)) < (proposed_log_post - current_log_post)) {
      current_beta <- proposed_beta
      current_g <- proposed_g
      current_log_post <- proposed_log_post
      current_tree <- beta_res$tree_fmat
    }

    current_encod <- my_encod(current_M)
    proposed_M <- Fmat_from_myencod(proposal_myencod(current_encod))
    M_res <- log_posterior_joint_timed(
      current_beta, proposed_M, coal_times, prior_sd,
      lambda = lambda, start_tree = current_tree, mode = mode, R = R
    )
    proposed_log_post_M <- M_res$log_posterior

    if (log(runif(1)) < (proposed_log_post_M - current_log_post)) {
      current_M <- proposed_M
      current_log_post <- proposed_log_post_M
      current_tree <- M_res$tree_fmat
    }

    chain_timings[[iter]] <- M_res$timings
    chain_timings[[iter]]$iteration_total <- proc.time()[["elapsed"]] - t_iter
    chain_timings[[iter]]$iteration_overhead <- chain_timings[[iter]]$iteration_total - chain_timings[[iter]]$tree_and_likelihood - chain_timings[[iter]]$sampleF - chain_timings[[iter]]$logdistance - chain_timings[[iter]]$bookkeeping
    cat("Iter", iter, "log_post", current_log_post, "\n")
  }

  chain_timings
}

run_timed_joint_chain_parallel <- function(
  M_init,
  beta_init,
  prior_sd,
  n_iter,
  coal_times,
  proposal_sd,
  lambda = 1,
  mode = "average",
  R = 10,
  cores = 2
) {
  chain_timings <- vector("list", n_iter)
  current_beta <- beta_init
  current_g <- log(current_beta)
  current_M <- M_init

  t_iter <- proc.time()[["elapsed"]]
  current_res <- log_posterior_joint_timed_parallel(
    current_beta, current_M, coal_times, prior_sd,
    lambda = lambda, start_tree = current_M, mode = mode, R = R, cores = cores
  )
  current_log_post <- current_res$log_posterior
  current_tree <- current_res$tree_fmat
  chain_timings[[1]] <- current_res$timings
  chain_timings[[1]]$iteration_total <- proc.time()[["elapsed"]] - t_iter
  chain_timings[[1]]$iteration_overhead <- chain_timings[[1]]$iteration_total - chain_timings[[1]]$sampleF - chain_timings[[1]]$logdistance - chain_timings[[1]]$tree_and_likelihood - chain_timings[[1]]$bookkeeping

  for (iter in 2:n_iter) {
    t_iter <- proc.time()[["elapsed"]]
    proposed_g <- rnorm(1, mean = current_g, sd = proposal_sd)
    proposed_beta <- exp(proposed_g)

    beta_res <- log_posterior_joint_timed_parallel(
      proposed_beta, current_M, coal_times, prior_sd,
      lambda = lambda, start_tree = current_tree, mode = mode, R = R, cores = cores
    )
    proposed_log_post <- beta_res$log_posterior

    if (log(runif(1)) < (proposed_log_post - current_log_post)) {
      current_beta <- proposed_beta
      current_g <- proposed_g
      current_log_post <- proposed_log_post
      current_tree <- beta_res$tree_fmat
    }

    current_encod <- my_encod(current_M)
    proposed_M <- Fmat_from_myencod(proposal_myencod(current_encod))
    M_res <- log_posterior_joint_timed_parallel(
      current_beta, proposed_M, coal_times, prior_sd,
      lambda = lambda, start_tree = current_tree, mode = mode, R = R, cores = cores
    )
    proposed_log_post_M <- M_res$log_posterior

    if (log(runif(1)) < (proposed_log_post_M - current_log_post)) {
      current_M <- proposed_M
      current_log_post <- proposed_log_post_M
      current_tree <- M_res$tree_fmat
    }

    chain_timings[[iter]] <- M_res$timings
    chain_timings[[iter]]$iteration_total <- proc.time()[["elapsed"]] - t_iter
    chain_timings[[iter]]$iteration_overhead <- chain_timings[[iter]]$iteration_total - chain_timings[[iter]]$sampleF - chain_timings[[iter]]$logdistance - chain_timings[[iter]]$tree_and_likelihood - chain_timings[[iter]]$bookkeeping
    cat("Iter", iter, "log_post", current_log_post, "\n")
  }

  chain_timings
}

n_tips <- 25
set.seed(123)
sequences <- as.phyDat(read.FASTA("sequences.fasta"))
M_true_tree <- read.tree("true_tree.newick")
inter_coal_times <- coalescent.intervals(M_true_tree)$interval.length
inter_coal_times[which(inter_coal_times <= 0.000001)] <- 0.0000001
true_coalescent_times <- cumsum(inter_coal_times)
M_true_fmat <- gen_Fmat(M_true_tree, tol = 8)

dm <- dist.hamming(sequences)
upgma_tree <- upgma(dm)
upgma_tree <- update_time(upgma_tree, true_coalescent_times)
M_upgma_fmat <- gen_Fmat(upgma_tree, tol = 8)

timings <- run_timed_joint_chain(
  M_init = M_upgma_fmat,
  beta_init = 1,
  prior_sd = 1,
  n_iter = 20,
  coal_times = true_coalescent_times,
  proposal_sd = 0.5,
  lambda = 1,
  mode = "max",
  R = 10
)

all_timings <- do.call(rbind, lapply(timings, function(x) unlist(x)))
print(colMeans(all_timings, na.rm = TRUE))

posterior_cols <- c("iteration_total", "sampleF", "tree_and_likelihood", "bookkeeping")
print(colMeans(all_timings[, posterior_cols, drop = FALSE], na.rm = TRUE))

timings_parallel <- run_timed_joint_chain_parallel(
  M_init = M_upgma_fmat,
  beta_init = 1,
  prior_sd = 1,
  n_iter = 20,
  coal_times = true_coalescent_times,
  proposal_sd = 0.5,
  lambda = 1,
  mode = "max",
  R = 10,
  cores = 2
)

all_timings_parallel <- do.call(rbind, lapply(timings_parallel, function(x) unlist(x)))
print(colMeans(all_timings_parallel[, posterior_cols, drop = FALSE], na.rm = TRUE))
