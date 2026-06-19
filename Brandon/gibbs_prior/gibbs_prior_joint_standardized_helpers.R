source("utils.R")
make_experiment_dirs <- function(experiment_id) {
  base <- file.path(output_root, experiment_id)
  dir.create(file.path(base, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base, "summaries"), recursive = TRUE, showWarnings = FALSE)
  base
}

save_path_for <- function(experiment_id, filename) {
  root <- get0(
    "output_root",
    ifnotfound = file.path("plots", "joint_standardized"),
    inherits = TRUE
  )
  root_relative <- sub("^plots/", "", root)
  file.path(root_relative, experiment_id, "figures", filename)
}

result_value_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_real_)
  }
  x
}

generate_proposal_tree_chain <- function(
    start_tree,
    n_trees,
    burn_in = 0,
    take_every = 1,
    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.matrix(start_tree)) {
    stop("start_tree must be an F-matrix.")
  }
  if (length(dim(start_tree)) != 2 || nrow(start_tree) != ncol(start_tree)) {
    stop("start_tree must be a square matrix.")
  }
  if (length(n_trees) != 1 || is.na(n_trees) || n_trees < 1 || n_trees != as.integer(n_trees)) {
    stop("n_trees must be a positive integer.")
  }
  if (length(burn_in) != 1 || is.na(burn_in) || burn_in < 0 || burn_in != as.integer(burn_in)) {
    stop("burn_in must be a non-negative integer.")
  }
  if (length(take_every) != 1 || is.na(take_every) || take_every < 1 || take_every != as.integer(take_every)) {
    stop("take_every must be a positive integer.")
  }

  current_tree <- start_tree
  tree_chain <- vector("list", n_trees)

  if (burn_in > 0) {
    for (step in seq_len(burn_in)) {
      current_tree <- Fmat_from_myencod(proposal_myencod(my_encod(current_tree)))
    }
  }

  tree_chain[[1]] <- current_tree

  if (n_trees > 1) {
    for (idx in 2:n_trees) {
      for (step in seq_len(take_every)) {
        current_tree <- Fmat_from_myencod(proposal_myencod(my_encod(current_tree)))
      }
      tree_chain[[idx]] <- current_tree
    }
  }

  tree_chain
}

generate_proposal_tree_chain_cache <- function(
    n_trees, n_tips = 10) {
  
  encoding_chain <- rEncod(m = n_trees, n = n_tips, distr="uniform")
  tree_chain <- vector("list", n_trees)
  for (i in 1:n_trees) {
    tree_chain[[i]] <- Fmat_from_myencod(encoding_chain[[i]])
  }


  precompute_tree_chain_distance_cache(tree_chain)
}

precompute_tree_chain_distance_cache <- function(tree_chain) {
  valid_tree_chain <- Filter(Negate(is.null), tree_chain)
  if (length(valid_tree_chain) == 0) {
    stop("tree_chain must contain at least one matrix.")
  }

  template <- valid_tree_chain[[1]]
  template_dim <- dim(template)
  if (length(template_dim) != 2 || template_dim[1] != template_dim[2]) {
    stop("tree_chain entries must be square matrices.")
  }

  lower_idx <- lower.tri(template, diag = TRUE)
  feature_count <- sum(lower_idx)

  tree_matrix <- t(vapply(
    valid_tree_chain,
    function(tree) {
      if (!identical(dim(tree), template_dim)) {
        stop("All tree matrices must have the same dimensions.")
      }
      tree[lower_idx]
    },
    numeric(feature_count)
  ))

  list(
    tree_matrix = tree_matrix,
    tree_ss = rowSums(tree_matrix * tree_matrix),
    lower_idx = lower_idx,
    template_dim = template_dim,
    n_trees = nrow(tree_matrix)
  )
}

tree_chain_squared_distances <- function(cache, M) {
  if (is.null(cache$tree_matrix) || is.null(cache$tree_ss) || is.null(cache$lower_idx)) {
    stop("cache must come from precompute_tree_chain_distance_cache().")
  }

  if (!identical(dim(M), cache$template_dim)) {
    stop("M must have the same dimensions as the cached tree matrices.")
  }

  M_vec <- M[cache$lower_idx]
  M_ss <- sum(M_vec * M_vec)

  as.numeric(cache$tree_ss - 2 * as.vector(cache$tree_matrix %*% M_vec) + M_ss)
}

compute_log_Z_est <- function(beta, M, cache, diam = 1) {
  distances <- tree_chain_squared_distances(cache, M)
  

  log_terms <- -(beta / diam) * distances
  if (any(!is.finite(log_terms))) {
    stop(
      "compute_log_Z_est() produced non-finite values: ",
      "beta=", beta,
      ", diam=", diam
  
    )
  }

  m <- max(log_terms)

  log_Z_est <- m + log(mean(exp(log_terms - m)))
  if (!is.finite(log_Z_est)) {
    stop(
      "compute_log_Z_est() returned a non-finite log_Z estimate: ",
      "beta=", beta,
      ", diam=", diam
    )
  }

  log_Z_est
}

plot_beta_histogram <- function(beta_chain, title, save_path) {
  df <- data.frame(beta = beta_chain)
  p <- ggplot(df, aes(x = beta)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    labs(title = title, x = "Beta", y = "Count") +
    theme_minimal(base_size = 14)
  print(p)
  ggsave(filename = paste0("plots/", save_path), plot = p, width = 8, height = 6, dpi = 150)
}

clade_posterior_comparison_plot <- function(
    chain_trees,
    beast_trees,
    true_tree,
    threshold = 0.01,
    title = "Clade posterior comparison",
    save_path = NULL) {
  chain_trees <- Filter(Negate(is.null), chain_trees)
  if (inherits(beast_trees, "multiPhylo")) {
    beast_trees <- as.list(beast_trees)
  }
  beast_trees <- Filter(Negate(is.null), beast_trees)

  class(chain_trees) <- "multiPhylo"
  class(beast_trees) <- "multiPhylo"

  n_chain <- length(chain_trees)
  n_beast <- length(beast_trees)

  chain_pp <- prop.part(chain_trees)
  beast_pp <- prop.part(beast_trees)
  true_pp  <- prop.part(true_tree)

  chain_labels <- attr(chain_pp, "labels")
  beast_labels <- attr(beast_pp, "labels")
  true_labels  <- attr(true_pp, "labels")

  n_tips <- length(chain_labels)

  clade_key <- function(idx, labels) paste(sort(labels[idx]), collapse = "|")
  is_internal <- function(idx, n) length(idx) >= 2 && length(idx) < n

  chain_keys <- vapply(chain_pp, function(c) if (is_internal(c, n_tips)) clade_key(c, chain_labels) else NA_character_, character(1))
  beast_keys <- vapply(beast_pp, function(c) if (is_internal(c, n_tips)) clade_key(c, beast_labels) else NA_character_, character(1))
  true_keys  <- vapply(true_pp,  function(c) if (is_internal(c, n_tips)) clade_key(c, true_labels)  else NA_character_, character(1))
  true_keys  <- true_keys[!is.na(true_keys)]

  chain_counts <- attr(chain_pp, "number")
  beast_counts <- attr(beast_pp, "number")

  chain_keep <- !is.na(chain_keys)
  beast_keep <- !is.na(beast_keys)

  chain_dict <- setNames(chain_counts[chain_keep] / n_chain, chain_keys[chain_keep])
  beast_dict <- setNames(beast_counts[beast_keep] / n_beast, beast_keys[beast_keep])

  all_keys <- union(names(chain_dict), names(beast_dict))
  p_chain <- unname(chain_dict[all_keys]); p_chain[is.na(p_chain)] <- 0
  p_beast <- unname(beast_dict[all_keys]); p_beast[is.na(p_beast)] <- 0
  in_true <- all_keys %in% true_keys

  l1_distance <- sum(abs(p_chain - p_beast))

  df <- data.frame(
    clade = all_keys,
    p_chain = p_chain,
    p_beast = p_beast,
    in_true = in_true,
    stringsAsFactors = FALSE
  )

  df_plot <- df[pmax(df$p_chain, df$p_beast) >= threshold, , drop = FALSE]

  p <- ggplot(df_plot, aes(x = p_beast, y = p_chain, color = in_true)) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    geom_abline(slope = 1, intercept = 0.2, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 1, intercept = -0.2, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.75, size = 2.2) +
    scale_color_manual(
      values = c("FALSE" = "tomato", "TRUE" = "forestgreen"),
      labels = c("FALSE" = "Not in true tree", "TRUE" = "In true tree")
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = title,
      subtitle = sprintf(
        "L1 distance over %d clades = %.3f | shown: max(p) >= %g (%d/%d clades)",
        nrow(df), l1_distance, threshold, nrow(df_plot), nrow(df)
      ),
      x = "BEAST posterior probability",
      y = "Chain posterior probability",
      color = ""
    ) +
    theme_minimal(base_size = 14)

  print(p)
  if (!is.null(save_path)) {
    ggsave(filename = paste0("plots/", save_path), plot = p, width = 8, height = 8, dpi = 150)
  }

  list(
    plot = p,
    l1_distance = l1_distance,
    threshold = threshold,
    df = df,
    df_plot = df_plot
  )
}

rf_trace_plot <- function(rf_chain,
                          title = "Robinson-Foulds distance to true tree (trace)",
                          save_path = NULL) {
  df <- data.frame(iter = seq_along(rf_chain), rf = rf_chain)
  p <- ggplot(df, aes(x = iter, y = rf)) +
    geom_line(color = "steelblue", alpha = 0.7) +
    labs(title = title, x = "Chain index", y = "Robinson-Foulds distance") +
    theme_minimal(base_size = 14)
  print(p)
  if (!is.null(save_path)) {
    ggsave(filename = paste0("plots/", save_path), plot = p, width = 10, height = 5, dpi = 150)
  }
  invisible(p)
}

rf_histogram_plot <- function(chain_trees, beast_trees, true_tree,
                              title = "Robinson-Foulds distance to true tree",
                              save_path = NULL) {
  chain_trees <- Filter(Negate(is.null), chain_trees)
  if (inherits(beast_trees, "multiPhylo")) {
    beast_trees <- as.list(beast_trees)
  }
  beast_trees <- Filter(Negate(is.null), beast_trees)

  rf_chain <- vapply(chain_trees, function(t) RF.dist(t, true_tree), numeric(1))
  rf_beast <- vapply(beast_trees, function(t) RF.dist(t, true_tree), numeric(1))

  df <- rbind(
    data.frame(rf = rf_chain, group = "Chain"),
    data.frame(rf = rf_beast, group = "BEAST")
  )

  p <- ggplot(df, aes(x = rf, fill = group)) +
    geom_histogram(position = "identity", alpha = 0.5, binwidth = 1) +
    scale_fill_manual(values = c("Chain" = "steelblue", "BEAST" = "tomato")) +
    labs(title = title, x = "Robinson-Foulds distance", y = "Count", fill = "Source") +
    theme_minimal(base_size = 14)

  print(p)
  if (!is.null(save_path)) {
    ggsave(filename = paste0("plots/", save_path), plot = p, width = 8, height = 6, dpi = 150)
  }

  list(rf_chain = rf_chain, rf_beast = rf_beast, plot = p)
}

joint_log_posterior <- function(beta, M, coal_times, prior_sd, log_Z_M_B = NULL, lambda = 1, start_tree = NULL,
                                mode = "max", sequences, n_tips) {
  tree_fmat <- sampleF(
    M,
    beta,
    iter = 2,
    diam = 1,
    startF = start_tree,
    take_every = 1
  )$chainF[[2]]

  log_boltzmann_dist <- Logdistance_prob(tree_fmat, M, beta = beta, diam = 1)
  ll_res <- log_likelihood_given_tree(
    tree_fmat = tree_fmat,
    coal_times = coal_times,
    sequences = sequences,
    mode = mode,
    R = 10
  )
  log_likelihood <- ll_res$log_likelihood
  log_rooted_given_unrooted <- log(2^(num_cherries(tree_fmat)) / factorial(n_tips))
  log_prior <- dnorm(log(beta), mean = 0, sd = prior_sd, log = TRUE)
  log_posterior <- log_likelihood / lambda + log_rooted_given_unrooted + log_boltzmann_dist + log_prior
  if (!is.null(log_Z_M_B)) {
    log_posterior <- log_posterior - log_Z_M_B
  }

  list(
    log_posterior = log_posterior,
    tree_fmat = tree_fmat,
    log_likelihood = log_likelihood,
    rooted_tree = ll_res$rooted_tree
  )
}


tree_log_posterior <- function(tree_fmat, beta, M, coal_times, lambda = 1,
                                mode = "max", sequences, n_tips, prior_sd = NULL,
                                log_Z_M_B = NULL) {

  log_boltzmann_dist <- Logdistance_prob(tree_fmat, M, beta = beta, diam = 1)
  ll_res <- log_likelihood_given_tree(
    tree_fmat = tree_fmat,
    coal_times = coal_times,
    sequences = sequences,
    mode = mode,
    R = 10
  )
  log_likelihood <- ll_res$log_likelihood
  log_rooted_given_unrooted <- log(2^(num_cherries(tree_fmat)) / factorial(n_tips))
  log_posterior <- log_likelihood / lambda + log_rooted_given_unrooted + log_boltzmann_dist
  if (!is.null(prior_sd)) {
    log_posterior <- log_posterior + dnorm(log(beta), mean = 0, sd = prior_sd, log = TRUE)
  }
  if (!is.null(log_Z_M_B)) {
    log_posterior <- log_posterior - log_Z_M_B
  }

  list(
    log_posterior = log_posterior,
    log_likelihood = log_likelihood,
    rooted_tree = ll_res$rooted_tree
  )
}

run_joint_mh_experiment <- function(
    experiment_id,
    init_M,
    coal_times,
    sequences,
    n_iter = 10000,
    beta_init = 1,
    sample_beta = TRUE,
    prior_sd = 1,
    proposal_sd = 0.5,
    lambda = 1,
    mode = "max",
    take_every = 1,
    include_log_Z = FALSE,
    include_log_Z_in_beta = FALSE,
    uniform_tree_cache = NULL) {
  n_tips <- ncol(init_M) + 1

  use_logZ_for_beta <- include_log_Z && include_log_Z_in_beta

  beta_chain <- numeric(n_iter)
  log_post_chain <- numeric(n_iter)
  M_chain <- vector("list", ceiling(n_iter / take_every))
  tree_chain <- vector("list", ceiling(n_iter / take_every))
  labeled_tree_chain <- vector("list", ceiling(n_iter / take_every))

  accepted_beta <- 0
  accepted_M <- 0

  current_beta <- beta_init
  current_g <- if (sample_beta) log(current_beta) else NULL
  current_M <- init_M

  if (is.null(uniform_tree_cache) && include_log_Z) {
    uniform_tree_cache <- generate_proposal_tree_chain_cache(
      n_trees = 100000,
      n_tips = n_tips
    )
  }

  current_log_Z_M_B <- if (include_log_Z) compute_log_Z_est(current_beta, current_M, uniform_tree_cache) else NULL

  current_res <- joint_log_posterior(
    current_beta,
    current_M,
    coal_times,
    prior_sd,
    log_Z_M_B = if (use_logZ_for_beta) current_log_Z_M_B else NULL,
    lambda = lambda,
    start_tree = current_M,
    mode = mode,
    sequences = sequences,
    n_tips = n_tips
  )
  current_tree <- current_res$tree_fmat
  current_labeled_tree <- current_res$rooted_tree
  current_log_post_beta <- current_res$log_posterior
  current_log_post_M <- if (include_log_Z) {
    tree_log_posterior(
      current_tree,
      current_beta,
      current_M,
      coal_times,
      lambda = lambda,
      mode = mode,
      sequences = sequences,
      n_tips = n_tips,
      prior_sd = prior_sd,
      log_Z_M_B = current_log_Z_M_B
    )$log_posterior
  } else {
    current_log_post_beta
  }
  current_log_post <- if (include_log_Z) current_log_post_M else current_log_post_beta

  beta_chain[1] <- current_beta
  log_post_chain[1] <- current_log_post
  M_chain[[1]] <- current_M
  tree_chain[[1]] <- current_tree
  labeled_tree_chain[[1]] <- current_labeled_tree

  for (iter in 2:n_iter) {
    if (sample_beta) {
      proposed_g <- rnorm(1, mean = current_g, sd = proposal_sd)
      proposed_beta <- exp(proposed_g)

      proposed_log_Z_for_beta <- if (include_log_Z) compute_log_Z_est(proposed_beta, current_M, uniform_tree_cache) else NULL

      beta_res <- joint_log_posterior(
        proposed_beta,
        current_M,
        coal_times,
        prior_sd,
        log_Z_M_B = if (use_logZ_for_beta) proposed_log_Z_for_beta else NULL,
        lambda = lambda,
        start_tree = current_tree,
        mode = mode,
        sequences = sequences,
        n_tips = n_tips
      )
      proposed_log_post_beta <- beta_res$log_posterior

      if (log(runif(1)) < (proposed_log_post_beta - current_log_post_beta)) {
        current_beta <- proposed_beta
        current_g <- proposed_g
        current_tree <- beta_res$tree_fmat
        current_labeled_tree <- beta_res$rooted_tree
        current_log_post_beta <- proposed_log_post_beta
        if (include_log_Z) {
          current_log_Z_M_B <- proposed_log_Z_for_beta
          current_log_post_M <- tree_log_posterior(
            current_tree,
            current_beta,
            current_M,
            coal_times,
            lambda = lambda,
            mode = mode,
            sequences = sequences,
            n_tips = n_tips,
            prior_sd = prior_sd,
            log_Z_M_B = current_log_Z_M_B
          )$log_posterior
        } else {
          current_log_post_M <- current_log_post_beta
        }
        accepted_beta <- accepted_beta + 1
      }
    }

    current_encod <- my_encod(current_M)
    proposed_M <- Fmat_from_myencod(proposal_myencod(current_encod))

    proposed_log_Z_M_B <- if (include_log_Z) compute_log_Z_est(current_beta, proposed_M, uniform_tree_cache) else NULL
    M_res <- joint_log_posterior(
      current_beta,
      proposed_M,
      coal_times,
      prior_sd,
      log_Z_M_B = proposed_log_Z_M_B,
      lambda = lambda,
      start_tree = current_tree,
      mode = mode,
      sequences = sequences,
      n_tips = n_tips
    )
    proposed_log_post_M <- M_res$log_posterior

    if (log(runif(1)) < (proposed_log_post_M - current_log_post_M)) {
      current_M <- proposed_M
      current_tree <- M_res$tree_fmat
      current_labeled_tree <- M_res$rooted_tree
      current_log_post_M <- proposed_log_post_M
      if (include_log_Z) {
        current_log_Z_M_B <- proposed_log_Z_M_B
      }
      current_log_post_beta <- tree_log_posterior(
        current_tree,
        current_beta,
        current_M,
        coal_times,
        lambda = lambda,
        mode = mode,
        sequences = sequences,
        n_tips = n_tips,
        prior_sd = prior_sd,
        log_Z_M_B = if (use_logZ_for_beta) current_log_Z_M_B else NULL
      )$log_posterior
      if (!include_log_Z) {
        current_log_post_M <- current_log_post_beta
      }
      accepted_M <- accepted_M + 1
    }

    current_log_post <- if (include_log_Z) current_log_post_M else current_log_post_beta
    beta_chain[iter] <- current_beta
    log_post_chain[iter] <- current_log_post

    if (iter %% take_every == 0) {
      slot <- iter / take_every
      M_chain[[slot]] <- current_M
      tree_chain[[slot]] <- current_tree
      labeled_tree_chain[[slot]] <- current_labeled_tree
    }

    per10 <- max(1, floor(n_iter / 10))
    if (iter %% per10 == 0) {
      cat(
        "Iteration:", iter,
        "Current Beta:", current_beta,
        "Log Posterior:", current_log_post,
        if (include_log_Z) paste("Current Log Z:", current_log_Z_M_B) else "",
        "\n"
      )
    }
  }

  valid_M_chain <- Filter(Negate(is.null), M_chain)
  average_M_raw <- Reduce(`+`, valid_M_chain) / length(valid_M_chain)
  average_M <- nearby_Fmat(average_M_raw)

  list(
    experiment_id = experiment_id,
    beta_chain = beta_chain,
    log_posterior = log_post_chain,
    M_chain = M_chain,
    tree_chain = tree_chain,
    labeled_tree_chain = labeled_tree_chain,
    average_M_raw = average_M_raw,
    average_M = average_M,
    acceptance_rate_beta = if (sample_beta) accepted_beta / (n_iter - 1) else 1,
    acceptance_rate_M = accepted_M / (n_iter - 1),
    lambda = lambda
  )
}




run_mean_mh_experiment <- function(
    experiment_id,
    init_M,
    coal_times,
    sequences,
    update_time = 50,
    recent_chain_time = 1000,
    n_iter = 100000,
    beta = 1,
    sample_beta = FALSE,
    prior_sd = 1,
    proposal_sd = 0.5,
    lambda = 1,
    mode = "max",
    take_every = 1,
    include_log_Z = FALSE,
    uniform_tree_cache = NULL) {
  n_tips <- ncol(init_M) + 1

  beta_chain <- numeric(n_iter)
  log_post_chain <- numeric(n_iter)
  M_chain <- vector("list", ceiling(n_iter / take_every))
  tree_chain <- vector("list", ceiling(n_iter / take_every))
  labeled_tree_chain <- vector("list", ceiling(n_iter / take_every))

  accepted_beta <- 0
  accepted_tree <- 0

  current_beta <- beta
  current_g <- if (sample_beta) log(current_beta) else NULL
  current_M <- init_M
  current_tree <- gen_Fmat(rcoal(n_tips), tol = 8)

  if (is.null(uniform_tree_cache) && include_log_Z) {
    uniform_tree_cache <- generate_proposal_tree_chain_cache(
      n_trees = 100000,
      n_tips = n_tips
    )
  }

  current_log_Z_M_B <- if (include_log_Z) compute_log_Z_est(current_beta, current_M, uniform_tree_cache) else NULL

  current_res <- tree_log_posterior(
      current_tree,
      current_beta,
      current_M,
      coal_times,
      lambda = lambda,
      mode = mode,
      sequences = sequences,
      n_tips = n_tips,
      prior_sd = if (sample_beta) prior_sd else NULL,
      log_Z_M_B = current_log_Z_M_B
    )
  current_log_post <- current_res$log_posterior
  current_labeled_tree <- current_res$rooted_tree

  beta_chain[1] <- current_beta
  log_post_chain[1] <- current_log_post
  M_chain[[1]] <- current_M
  tree_chain[[1]] <- current_tree
  labeled_tree_chain[[1]] <- current_labeled_tree
 
  for (iter in 2:n_iter) {
    if (sample_beta) {
      proposed_g <- rnorm(1, mean = current_g, sd = proposal_sd)
      proposed_beta <- exp(proposed_g)
      proposed_log_Z_M_B <- if (include_log_Z) compute_log_Z_est(proposed_beta, current_M, uniform_tree_cache) else NULL

      beta_res <- tree_log_posterior(
        current_tree,
        proposed_beta,
        current_M,
        coal_times,
        lambda = lambda,
        mode = mode,
        sequences = sequences,
        n_tips = n_tips,
        prior_sd = prior_sd,
        log_Z_M_B = proposed_log_Z_M_B
      )
      proposed_log_post_beta <- beta_res$log_posterior

      if (log(runif(1)) < (proposed_log_post_beta - current_log_post)) {
        current_beta <- proposed_beta
        current_g <- proposed_g
        current_log_post <- proposed_log_post_beta
        current_log_Z_M_B <- proposed_log_Z_M_B
        current_labeled_tree <- beta_res$rooted_tree
        accepted_beta <- accepted_beta + 1
      }
    }

    current_encod <- my_encod(current_tree)

    proposed_tree <- Fmat_from_myencod(proposal_myencod(current_encod))


    tree_res <- tree_log_posterior(
      proposed_tree,
      current_beta,
      current_M,
      coal_times,
      lambda = lambda,
      mode = mode,
      sequences = sequences,
      n_tips = n_tips,
      prior_sd = if (sample_beta) prior_sd else NULL,
      log_Z_M_B = current_log_Z_M_B
    )
    proposed_log_post <- tree_res$log_posterior

    if (log(runif(1)) < (proposed_log_post - current_log_post)) {
      current_log_post <- proposed_log_post
      current_tree <- proposed_tree
      current_labeled_tree <- tree_res$rooted_tree
      accepted_tree <- accepted_tree + 1
    }
    if (iter %% update_time == 0) {
      valid_tree_chain <- c(Filter(Negate(is.null), tree_chain), list(current_tree))
      recent_tree_chain <- tail(valid_tree_chain, min(recent_chain_time, length(valid_tree_chain)))
      current_M <- Reduce(`+`, valid_tree_chain) / length(valid_tree_chain)
      current_log_Z_M_B <- if (include_log_Z) compute_log_Z_est(current_beta, current_M, uniform_tree_cache) else NULL
      current_res <- tree_log_posterior(
        current_tree,
        current_beta,
        current_M,
        coal_times,
        lambda = lambda,
        mode = mode,
        sequences = sequences,
        n_tips = n_tips,
        prior_sd = if (sample_beta) prior_sd else NULL,
        log_Z_M_B = current_log_Z_M_B
      )
      current_log_post <- current_res$log_posterior
      current_labeled_tree <- current_res$rooted_tree
    }

    beta_chain[iter] <- current_beta
    log_post_chain[iter] <- current_log_post

    if (iter %% take_every == 0) {
      slot <- iter / take_every
      M_chain[[slot]] <- current_M
      tree_chain[[slot]] <- current_tree
      labeled_tree_chain[[slot]] <- current_labeled_tree
    }

    per10 <- max(1, floor(n_iter / 10))
    if (iter %% per10 == 0) {
      cat(
        "Iteration:", iter,
        "Current Beta:", current_beta,
        "Log Posterior:", current_log_post,
        if (include_log_Z) paste("Current Log Z:", current_log_Z_M_B) else "",
        "\n"
      )
    }
  }

  valid_M_chain <- Filter(Negate(is.null), M_chain)
  average_M_raw <- Reduce(`+`, valid_M_chain) / length(valid_M_chain)
  average_M <- nearby_Fmat(average_M_raw)

  list(
    experiment_id = experiment_id,
    beta_chain = beta_chain,
    log_posterior = log_post_chain,
    M_chain = M_chain,
    tree_chain = tree_chain,
    labeled_tree_chain = labeled_tree_chain,
    average_M_raw = average_M_raw,
    average_M = average_M,
    acceptance_rate_beta = if (sample_beta) accepted_beta / (n_iter - 1) else 1,
    acceptance_rate_M = 1,
    acceptance_rate_tree = accepted_tree / (n_iter - 1),
    lambda = lambda
  )
}
build_prefix_summaries <- function(
    M_chain,
    F_chain,
    beta_chain,
    n_splits = 5,
    burn_in = 0,
    subsample_size = 500,
    seed = 42) {
  valid_M_chain <- Filter(Negate(is.null), M_chain)
  valid_F_chain <- Filter(Negate(is.null), F_chain)
  n <- min(length(valid_M_chain), length(valid_F_chain), length(beta_chain))

  if (n <= burn_in) {
    stop("burn_in leaves no samples in the chains.")
  }

  valid_M_chain <- valid_M_chain[(burn_in + 1):n]
  valid_F_chain <- valid_F_chain[(burn_in + 1):n]
  valid_beta_chain <- beta_chain[(burn_in + 1):n]
  n <- length(valid_beta_chain)

  split_ends <- unique(pmax(1, floor(seq_len(n_splits) * n / n_splits)))

  cat("Building prefix summaries for", length(split_ends), "splits.\n")
  pb <- progress_bar$new(
    format = "  prefix summaries [:bar] :percent eta: :eta",
    total = length(split_ends),
    clear = FALSE,
    width = 60,
    show_after = 0
  )

  summaries <- vector("list", length(split_ends))
  for (i in seq_along(split_ends)) {
    prefix_end <- split_ends[i]
    prefix_M <- valid_M_chain[seq_len(prefix_end)]
    prefix_F <- valid_F_chain[seq_len(prefix_end)]
    prefix_beta <- valid_beta_chain[seq_len(prefix_end)]

    mean_M_raw <- Reduce(`+`, prefix_M) / length(prefix_M)
    mean_M <- nearby_Fmat(mean_M_raw)
    mean_beta <- mean(prefix_beta)

    sample_n <- min(subsample_size, length(prefix_F))
    set.seed(seed + i)
    idx <- sample(length(prefix_F), sample_n)

    summaries[[i]] <- list(
      split_index = i,
      split_fraction = i / length(split_ends),
      prefix_end = prefix_end,
      mean_M = mean_M,
      mean_M_raw = mean_M_raw,
      mean_beta = mean_beta,
      empirical_F = prefix_F[idx],
      split_label = paste0(i, "/", length(split_ends))
    )

    pb$tick()
  }
  pb$terminate()

  summaries
}

generate_sampleF_reference <- function(
    mean_M,
    mean_beta,
    n_samples = 10000,
    take_every = 10,
    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  iter <- n_samples * take_every
  cat(
    "Generating SampleF reference with",
    n_samples,
    "samples (take_every =",
    take_every,
    ", beta =",
    round(mean_beta, 4),
    ").\n"
  )
  ref_chain <- sampleF(
    M = mean_M,
    beta = mean_beta,
    iter = iter,
    diam = 1,
    startF = mean_M,
    take_every = take_every
  )$chainF

  Filter(Negate(is.null), ref_chain)
}

compute_tree_mds_panel <- function(
    trees_empirical,
    trees_reference,
    true_tree,
    mean_M,
    split_label,
    reference_label) {
  trees_all <- c(trees_empirical, trees_reference, list(true_tree), list(mean_M))
  n <- length(trees_all)
  D <- matrix(0, n, n)

  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      d <- distance_Fmat(trees_all[[i]], trees_all[[j]], dist = "l2")^2
      D[i, j] <- d
      D[j, i] <- d
    }
  }

  mds <- cmdscale(as.dist(D), k = 2, eig = TRUE)
  coords <- as.data.frame(mds$points)
  colnames(coords) <- c("Dim1", "Dim2")

  labels <- c(
    rep("Empirical Chain", length(trees_empirical)),
    rep(reference_label, length(trees_reference)),
    "True M",
    "Average M"
  )
  coords$group <- factor(
    labels,
    levels = c("Empirical Chain", reference_label, "True M", "Average M")
  )
  coords$panel <- split_label
  coords$label <- c(
    rep("", length(trees_empirical) + length(trees_reference)),
    "True M",
    "Average M"
  )

  coords
}

prefix_tree_mds_comparison_plot <- function(
    prefix_summaries,
    true_tree,
    reference_mode = c("sampleF", "beast"),
    beast_fmats = NULL,
    sampleF_n_samples = 10000,
    sampleF_take_every = 10,
    plot_subsample_size = 500,
    seed = 42,
    title = "Prefix Tree MDS Comparison",
    save_path = NULL) {
  reference_mode <- match.arg(reference_mode)
  cat("Creating prefix MDS comparison against", reference_mode, "reference.\n")
  panel_coords <- vector("list", length(prefix_summaries))

  for (i in seq_along(prefix_summaries)) {
    prefix_summary <- prefix_summaries[[i]]
    cat(
      "  split",
      prefix_summary$split_label,
      "- prefix end =",
      prefix_summary$prefix_end,
      "\n"
    )

    if (reference_mode == "sampleF") {
      reference_fmats <- generate_sampleF_reference(
        mean_M = prefix_summary$mean_M,
        mean_beta = prefix_summary$mean_beta,
        n_samples = sampleF_n_samples,
        take_every = sampleF_take_every,
        seed = seed + i
      )
      reference_label <- "SampleF Reference"
    } else {
      if (is.null(beast_fmats)) {
        stop("beast_fmats must be provided when reference_mode = 'beast'.")
      }
      reference_fmats <- beast_fmats
      reference_label <- "BEAST Reference"
    }

    compare_n <- min(
      length(prefix_summary$empirical_F),
      length(reference_fmats),
      plot_subsample_size
    )

    set.seed(seed + 100 * i)
    empirical_idx <- sample(length(prefix_summary$empirical_F), compare_n)
    reference_idx <- sample(length(reference_fmats), compare_n)

    panel_coords[[i]] <- compute_tree_mds_panel(
      trees_empirical = prefix_summary$empirical_F[empirical_idx],
      trees_reference = reference_fmats[reference_idx],
      true_tree = true_tree,
      mean_M = prefix_summary$mean_M,
      split_label = prefix_summary$split_label,
      reference_label = reference_label
    )
  }

  coords <- do.call(rbind, panel_coords)

  p <- ggplot(coords, aes(x = Dim1, y = Dim2, color = group)) +
    geom_point(size = 1.8, alpha = 0.7) +
    geom_text(
      data = coords[coords$label != "", , drop = FALSE],
      aes(label = label),
      color = "black",
      vjust = -0.8,
      show.legend = FALSE
    ) +
    facet_wrap(~panel, scales = "free") +
    scale_color_manual(values = c("steelblue", "tomato", "green", "black")) +
    labs(title = title, x = "MDS 1", y = "MDS 2", color = "Group") +
    theme_minimal(base_size = 14)

  print(p)

  if (!is.null(save_path)) {
    ggsave(filename = paste0("plots/", save_path), plot = p, width = 12, height = 8, dpi = 150)
  }
}

render_mh_experiment <- function(
    result,
    experiment_id,
    comparison_matrices,
    comparison_labels,
    beast_compare_fmats,
    coal_times,
    sequences,
    sub_samp = 1000,
    seed = 42,
    prefix_splits = 5,
    prefix_burn_in = 0,
    sampleF_n_samples = 10000,
    sampleF_take_every = 10,
    true_tree = NULL,
    beast_compare_trees = NULL) {
  make_experiment_dirs(experiment_id)
  cat("Rendering experiment outputs for", experiment_id, "\n")

  plot_beta_histogram(
    result$beta_chain,
    title = paste0(experiment_id, ": posterior of beta"),
    save_path = save_path_for(experiment_id, "beta_histogram.png")
  )

  matrix_MDS_comparison_plot(
    matrices = comparison_matrices,
    labels = comparison_labels,
    chain_matrices = result$M_chain,
    sub_samp = sub_samp,
    seed = seed,
    title = paste0(experiment_id, ": matrix MDS"),
    save_path = save_path_for(experiment_id, "mds_regular.png")
  )

  matrix_MDS_comparison_plot(
    matrices = comparison_matrices,
    labels = comparison_labels,
    chain_matrices = result$M_chain,
    sub_samp = sub_samp,
    seed = seed,
    title = paste0(experiment_id, ": matrix MDS by log likelihood"),
    save_path = save_path_for(experiment_id, "mds_log_likelihood.png"),
    color_by_log_likelihood = TRUE,
    coal_times = coal_times,
    sequences = sequences,
    log_likelihood_mode = "max"
  )

  tree_chain <- Filter(Negate(is.null), result$tree_chain)
  tree_sample_n <- min(sub_samp, length(tree_chain), length(beast_compare_fmats))
  set.seed(seed)
  gen_idx <- sample(length(tree_chain), tree_sample_n)
  beast_idx <- sample(length(beast_compare_fmats), tree_sample_n)
  trees_gen <- tree_chain[gen_idx]
  trees_data <- beast_compare_fmats[beast_idx]

  tree_MDS_comparison_plot(
    trees_gen,
    trees_data,
    M_true_fmat,
    M_tree = result$average_M,
    title = paste0(experiment_id, ": tree MDS vs BEAST"),
    save_path = save_path_for(experiment_id, "tree_mds_vs_beast.png")
  )

  tree_histogram_comparison_plot(
    trees_gen,
    trees_data,
    M_true_fmat,
    title = paste0(experiment_id, ": tree histogram vs BEAST"),
    save_path = save_path_for(experiment_id, "tree_histogram_vs_beast.png")
  )

  clade_l1_distance <- NA_real_
  if (!is.null(result$labeled_tree_chain) && !is.null(true_tree) && !is.null(beast_compare_trees)) {
    rf_res <- rf_histogram_plot(
      chain_trees = result$labeled_tree_chain,
      beast_trees = beast_compare_trees,
      true_tree = true_tree,
      title = paste0(experiment_id, ": RF distance to true tree"),
      save_path = save_path_for(experiment_id, "rf_histogram_vs_beast.png")
    )
    rf_trace_plot(
      rf_chain = rf_res$rf_chain,
      title = paste0(experiment_id, ": RF distance to true tree (trace)"),
      save_path = save_path_for(experiment_id, "rf_trace.png")
    )
    clade_res <- clade_posterior_comparison_plot(
      chain_trees = result$labeled_tree_chain,
      beast_trees = beast_compare_trees,
      true_tree = true_tree,
      threshold = 0.01,
      title = paste0(experiment_id, ": clade posterior vs BEAST"),
      save_path = save_path_for(experiment_id, "clade_posterior_vs_beast.png")
    )
    clade_l1_distance <- clade_res$l1_distance
  } else {
    cat("Skipping RF / clade diagnostics for", experiment_id,
        "- need labeled_tree_chain, true_tree, and beast_compare_trees.\n")
  }

  cat("Running prefix-over-time diagnostics for", experiment_id, "\n")
  prefix_summaries <- build_prefix_summaries(
    M_chain = result$M_chain,
    F_chain = result$tree_chain,
    beta_chain = result$beta_chain,
    n_splits = prefix_splits,
    burn_in = prefix_burn_in,
    subsample_size = sub_samp,
    seed = seed
  )

  prefix_tree_mds_comparison_plot(
    prefix_summaries = prefix_summaries,
    true_tree = M_true_fmat,
    reference_mode = "sampleF",
    sampleF_n_samples = sampleF_n_samples,
    sampleF_take_every = sampleF_take_every,
    plot_subsample_size = sub_samp,
    seed = seed,
    title = paste0(experiment_id, ": prefix tree MDS vs SampleF"),
    save_path = save_path_for(experiment_id, "prefix_tree_mds_vs_samplef.png")
  )

  prefix_tree_mds_comparison_plot(
    prefix_summaries = prefix_summaries,
    true_tree = M_true_fmat,
    reference_mode = "beast",
    beast_fmats = beast_compare_fmats,
    plot_subsample_size = sub_samp,
    seed = seed,
    title = paste0(experiment_id, ": prefix tree MDS vs BEAST"),
    save_path = save_path_for(experiment_id, "prefix_tree_mds_vs_beast.png")
  )

  matrix_chain_ess_heatmap(
    result$M_chain,
    title = paste0(experiment_id, ": M-chain ESS"),
    save_path = save_path_for(experiment_id, "m_chain_ess.png")
  )

  matrix_chain_ess_heatmap(
    result$M_chain,
    normalize = TRUE,
    per_n_samples = 1000,
    title = paste0(experiment_id, ": M-chain ESS per 1000 samples"),
    save_path = save_path_for(experiment_id, "m_chain_ess_per_1000.png")
  )

  matrix_chain_ess_heatmap(
    result$tree_chain,
    title = paste0(experiment_id, ": tree-chain ESS"),
    save_path = save_path_for(experiment_id, "tree_chain_ess.png")
  )

  matrix_chain_ess_heatmap(
    result$tree_chain,
    normalize = TRUE,
    per_n_samples = 1000,
    title = paste0(experiment_id, ": tree-chain ESS per 1000 samples"),
    save_path = save_path_for(experiment_id, "tree_chain_ess_per_1000.png")
  )

  data.frame(
    experiment = experiment_id,
    lambda = result$lambda,
    beta_mean = mean(result$beta_chain),
    beta_sd = sd(result$beta_chain),
    acceptance_rate_beta = result_value_or_na(result$acceptance_rate_beta),
    acceptance_rate_M = result_value_or_na(result$acceptance_rate_M),
    average_M_distance_to_true = distance_Fmat(result$average_M, M_true_fmat),
    clade_l1_distance_vs_beast = clade_l1_distance
  )
}

render_beast_experiment <- function(
    beast_fmats,
    experiment_id,
    coal_times,
    sequences,
    sub_samp = 1000,
    seed = 42) {
  make_experiment_dirs(experiment_id)

  beast_average_M_raw <- Reduce(`+`, beast_fmats) / length(beast_fmats)
  beast_average_M <- nearby_Fmat(beast_average_M_raw)
  comparison_matrices <- list(M_true_fmat, beast_average_M)
  comparison_labels <- c("True M", "BEAST Average M")

  matrix_MDS_comparison_plot(
    matrices = comparison_matrices,
    labels = comparison_labels,
    chain_matrices = beast_fmats,
    sub_samp = sub_samp,
    seed = seed,
    title = paste0(experiment_id, ": matrix MDS"),
    save_path = save_path_for(experiment_id, "mds_regular.png")
  )

  matrix_MDS_comparison_plot(
    matrices = comparison_matrices,
    labels = comparison_labels,
    chain_matrices = beast_fmats,
    sub_samp = sub_samp,
    seed = seed,
    title = paste0(experiment_id, ": matrix MDS by log likelihood"),
    save_path = save_path_for(experiment_id, "mds_log_likelihood.png"),
    color_by_log_likelihood = TRUE,
    coal_times = coal_times,
    sequences = sequences,
    log_likelihood_mode = "max"
  )

  matrix_chain_ess_heatmap(
    beast_fmats,
    title = paste0(experiment_id, ": BEAST ESS"),
    save_path = save_path_for(experiment_id, "ess.png")
  )

  matrix_chain_ess_heatmap(
    beast_fmats,
    normalize = TRUE,
    per_n_samples = 1000,
    title = paste0(experiment_id, ": BEAST ESS per 1000 samples"),
    save_path = save_path_for(experiment_id, "ess_per_1000.png")
  )

  data.frame(
    experiment = experiment_id,
    lambda = NA_real_,
    beta_mean = NA_real_,
    beta_sd = NA_real_,
    acceptance_rate_beta = NA_real_,
    acceptance_rate_M = NA_real_,
    average_M_distance_to_true = distance_Fmat(beast_average_M, M_true_fmat),
    clade_l1_distance_vs_beast = NA_real_
  )
}
