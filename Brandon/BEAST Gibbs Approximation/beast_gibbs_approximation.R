library(ape)
library(fmatrix)
library(phylodyn)

launch_dir <- normalizePath(getwd(), mustWork = FALSE)

find_project_root <- function() {
  candidates <- unique(normalizePath(
    c(
      getwd(),
      file.path(getwd(), ".."),
      file.path(getwd(), "..", "..")
    ),
    mustWork = FALSE
  ))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "Brandon", "Fmatrix Bernoulli", "utils.R")) &&
        file.exists(file.path(candidate, "Brandon", "gibbs_prior", "VB_Utils.R"))) {
      return(candidate)
    }
  }

  stop("Could not find the Enertree project root from ", getwd())
}

project_root <- find_project_root()
if (!identical(normalizePath(getwd(), mustWork = FALSE), project_root)) {
  setwd(project_root)
}

project_path <- function(...) {
  file.path(project_root, ...)
}

is_absolute_path <- function(path) {
  grepl("^/", path) || grepl("^[A-Za-z]:[\\\\/]", path)
}

resolve_existing_path <- function(path) {
  candidates <- path
  if (!is_absolute_path(path)) {
    candidates <- c(
      candidates,
      file.path(launch_dir, path),
      project_path("Brandon", "BEAST Gibbs Approximation", path),
      project_path("Brandon", "BEAST Gibbs Approximation", "inputs", path),
      project_path("Brandon", "Fmatrix Bernoulli", path)
    )
  }

  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }
  path
}

source(project_path("Brandon", "Fmatrix Bernoulli", "utils.R"))
source(project_path("Brandon", "gibbs_prior", "utils.R"))
source(project_path("Brandon", "gibbs_prior", "gibbs_prior_joint_standardized_helpers.R"), chdir = TRUE)
source(project_path("Brandon", "gibbs_prior", "VB_Utils.R"))

gen_Fmat <- phylodyn:::gen_Fmat

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

get_exact_fmat_support <- function(num_tips) {
  flist_name <- paste0("F.list", num_tips)
  if (!exists(flist_name, mode = "list", inherits = TRUE)) {
    stop(
      "Exact F-matrix list not available: ", flist_name, ". ",
      "For this experiment, generate/use BEAST data with a tip count that has ",
      "a matching fmatrix::F.list<num_tips> object."
    )
  }
  get(flist_name, inherits = TRUE)
}

fit_beta_given_M <- function(beast_fmats,
                             M,
                             all_Fmats,
                             beta_bounds = c(1e-8, 1e4),
                             root_tol = 1e-8) {
  if (!is.list(beast_fmats) || length(beast_fmats) == 0) {
    stop("beast_fmats must be a non-empty list of F-matrices.")
  }
  if (!is.matrix(M)) {
    stop("M must be a matrix.")
  }
  if (!is.list(all_Fmats) || length(all_Fmats) == 0) {
    stop("all_Fmats must be a non-empty list of F-matrices.")
  }
  if (length(beta_bounds) != 2 || any(!is.finite(beta_bounds)) ||
      any(beta_bounds <= 0) || beta_bounds[1] >= beta_bounds[2]) {
    stop("beta_bounds must be two increasing positive finite values.")
  }

  beast_cache <- phylodyn:::precompute_tree_chain_distance_cache(beast_fmats)
  support_cache <- phylodyn:::precompute_tree_chain_distance_cache(all_Fmats)

  beast_d2 <- phylodyn:::tree_chain_squared_distances(beast_cache, M)
  support_d2 <- phylodyn:::tree_chain_squared_distances(support_cache, M)

  target_d2 <- mean(beast_d2)

  expected_d2 <- function(beta) {
    log_w <- -beta * support_d2
    m <- max(log_w)
    w <- exp(log_w - m)
    sum(w * support_d2) / sum(w)
  }

  root_fn <- function(log_beta) {
    expected_d2(exp(log_beta)) - target_d2
  }

  lower_log_beta <- log(beta_bounds[1])
  upper_log_beta <- log(beta_bounds[2])
  lower_value <- root_fn(lower_log_beta)
  upper_value <- root_fn(upper_log_beta)

  boundary <- NA_character_
  if (abs(lower_value) <= root_tol) {
    log_beta_hat <- lower_log_beta
    boundary <- "lower"
  } else if (abs(upper_value) <= root_tol) {
    log_beta_hat <- upper_log_beta
    boundary <- "upper"
  } else if (lower_value < 0) {
    log_beta_hat <- lower_log_beta
    boundary <- "lower_target_more_dispersed_than_uniform"
  } else if (upper_value > 0) {
    log_beta_hat <- upper_log_beta
    boundary <- "upper_target_tighter_than_beta_bound"
  } else {
    log_beta_hat <- uniroot(
      f = root_fn,
      interval = c(lower_log_beta, upper_log_beta),
      tol = root_tol
    )$root
  }

  beta_hat <- exp(log_beta_hat)

  list(
    beta = beta_hat,
    log_beta = log_beta_hat,
    target_d2 = target_d2,
    expected_d2 = expected_d2(beta_hat),
    uniform_expected_d2 = expected_d2(0),
    min_support_d2 = min(support_d2),
    max_support_d2 = max(support_d2),
    boundary = boundary,
    beast_d2 = beast_d2,
    support_d2 = support_d2,
    support_cache = support_cache
  )
}

log_prob_gibbs_batch <- function(F_mats, M, beta, log_Z = NULL) {
  cache <- phylodyn:::precompute_tree_chain_distance_cache(F_mats)
  d2 <- phylodyn:::tree_chain_squared_distances(cache, M)

  if (is.null(log_Z)) {
    stop("log_Z must be supplied so log probabilities use the same normalizer.")
  }

  -beta * d2 - log_Z
}

coal_times_from_tree <- function(tree, min_interval = 1e-6) {
  inter_coal_times <- coalescent.intervals(tree)$interval.length
  inter_coal_times[inter_coal_times <= min_interval] <- min_interval
  cumsum(inter_coal_times)
}

tree_from_fmat_with_true_labels <- function(Fmat, coal_times, tip_labels) {
  num_tips <- ncol(Fmat) + 1

  if (length(tip_labels) != num_tips) {
    stop("tip_labels length does not match F-matrix tip count.")
  }

  full_Fmat <- rbind(rep(0, num_tips), cbind(rep(0, num_tips - 1), Fmat))
  tree <- phylodyn:::tree_from_F(full_Fmat, coal_times)
  tree$tip.label <- tip_labels
  reorder(tree, "postorder")
}

write_fmats_as_true_labeled_trees <- function(F_mats,
                                              true_tree,
                                              output_dir,
                                              prefix = "generated_gibbs") {
  coal_times <- coal_times_from_tree(true_tree)
  tip_labels <- true_tree$tip.label

  trees <- lapply(
    F_mats,
    tree_from_fmat_with_true_labels,
    coal_times = coal_times,
    tip_labels = tip_labels
  )
  class(trees) <- "multiPhylo"

  newick_path <- file.path(output_dir, paste0(prefix, "_true_labeled.newick"))
  nexus_path <- file.path(output_dir, paste0(prefix, "_true_labeled.trees"))
  label_order_path <- file.path(output_dir, paste0(prefix, "_true_label_order.csv"))
  coal_times_path <- file.path(output_dir, paste0(prefix, "_true_coal_times.csv"))

  write.tree(trees, file = newick_path)
  write.nexus(trees, file = nexus_path)
  write.csv(
    data.frame(tip_index = seq_along(tip_labels), tip_label = tip_labels),
    file = label_order_path,
    row.names = FALSE
  )
  write.csv(
    data.frame(coalescent_index = seq_along(coal_times), coal_time = coal_times),
    file = coal_times_path,
    row.names = FALSE
  )

  list(
    trees = trees,
    coal_times = coal_times,
    tip_labels = tip_labels,
    newick_path = newick_path,
    nexus_path = nexus_path,
    label_order_path = label_order_path,
    coal_times_path = coal_times_path
  )
}

fit_and_generate_gibbs_from_fmats <- function(F_mats,
                                              generated_sample_size = length(F_mats),
                                              seed = 1,
                                              beta_bounds = c(1e-8, 1e4)) {
  if (!is.list(F_mats) || length(F_mats) < 1) {
    stop("F_mats must be a non-empty list of F-matrices.")
  }

  fmat_dim <- nrow(F_mats[[1]])
  num_tips <- fmat_dim + 1
  all_Fmats <- get_exact_fmat_support(num_tips)

  if (!all(vapply(F_mats, function(F_mat) identical(dim(F_mat), dim(all_Fmats[[1]])), logical(1)))) {
    stop("BEAST F-matrix dimensions do not match exact support for ", num_tips, " tips.")
  }

  beast_mean <- mean_fmatrix(F_mats, project_to_fspace = FALSE)
  M_hat <- beast_mean$raw

  beta_fit <- fit_beta_given_M(
    beast_fmats = F_mats,
    M = M_hat,
    all_Fmats = all_Fmats,
    beta_bounds = beta_bounds
  )

  log_Z <- phylodyn:::compute_log_Z_est(
    beta = beta_fit$beta,
    M = M_hat,
    cache = beta_fit$support_cache
  )

  generated_info <- sample_from_prior(
    M = M_hat,
    b = beta_fit$beta,
    num_samps = generated_sample_size,
    all_Fmats = all_Fmats,
    seed = seed,
    return_info = TRUE
  )

  generated_fmats <- generated_info$samples
  generated_mean <- mean_fmatrix(generated_fmats, project_to_fspace = FALSE)

  beast_logq <- log_prob_gibbs_batch(F_mats, M = M_hat, beta = beta_fit$beta, log_Z = log_Z)
  generated_logq <- log_prob_gibbs_batch(generated_fmats, M = M_hat, beta = beta_fit$beta, log_Z = log_Z)

  summary <- data.frame(
    fmat_dim = fmat_dim,
    num_tips = num_tips,
    beast_sample_size = length(F_mats),
    generated_sample_size = generated_sample_size,
    exact_support_size = length(all_Fmats),
    unique_beast_fmats = count_unique_fmats(F_mats),
    unique_generated_fmats = count_unique_fmats(generated_fmats),
    beta_hat = beta_fit$beta,
    log_beta_hat = beta_fit$log_beta,
    beta_fit_boundary = ifelse(is.na(beta_fit$boundary), "", beta_fit$boundary),
    target_mean_d2 = beta_fit$target_d2,
    fitted_expected_d2 = beta_fit$expected_d2,
    uniform_expected_d2 = beta_fit$uniform_expected_d2,
    min_support_d2 = beta_fit$min_support_d2,
    max_support_d2 = beta_fit$max_support_d2,
    exact_prior_prob_sum = generated_info$prob_sum,
    exact_prior_ess = generated_info$ess,
    mean_logq_beast = mean(beast_logq),
    sd_logq_beast = sd(beast_logq),
    mean_logq_generated = mean(generated_logq),
    sd_logq_generated = sd(generated_logq),
    mean_fmat_l2_distance = distance_Fmat(
      beast_mean$raw,
      generated_mean$raw,
      dist = "l2"
    )
  )

  list(
    M_hat = M_hat,
    beta_fit = beta_fit,
    log_Z = log_Z,
    all_Fmats = all_Fmats,
    generated_fmats = generated_fmats,
    generated_indices = generated_info$indices,
    generated_probs = generated_info$probs,
    beast_mean = beast_mean,
    generated_mean = generated_mean,
    beast_logq = beast_logq,
    generated_logq = generated_logq,
    summary = summary
  )
}

run_beast_gibbs_approximation <- function(
    trees_path = Sys.getenv(
      "GIBBS_BEAST_TREES",
      file.path(
        "Brandon", "BEAST Gibbs Approximation", "inputs",
        "beast_trees_summarized-sequences.trees"
      )
    ),
    true_tree_path = Sys.getenv(
      "GIBBS_TRUE_TREE",
      file.path("Brandon", "BEAST Gibbs Approximation", "inputs", "true_tree.newick")
    ),
    experiment_id = "beast_vs_gibbs",
    generated_sample_size = NULL,
    plot_subsample_cap = 500,
    seed = 1) {

  trees_path <- resolve_existing_path(trees_path)
  true_tree_path <- resolve_existing_path(true_tree_path)

  if (!file.exists(trees_path)) {
    stop(
      "Could not find BEAST trees file: ", trees_path, "\n",
      "Either place your 10-tip .trees file there or run with, for example:\n",
      "GIBBS_BEAST_TREES=/path/to/your/file.trees Rscript 'Brandon/BEAST Gibbs Approximation/beast_gibbs_approximation.R'"
    )
  }
  if (!file.exists(true_tree_path)) {
    stop("Could not find true tree file: ", true_tree_path)
  }

  output_root <- file.path("Brandon", "BEAST Gibbs Approximation", "results")
  figures_dir <- file.path(output_root, "figures")
  summaries_dir <- file.path(output_root, "summaries")

  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(summaries_dir, recursive = TRUE, showWarnings = FALSE)

  save_path_for_plot_function <- function(filename) {
    file.path(figures_dir, filename)
  }

  move_helper_plot_outputs <- function(filenames) {
    for (filename in filenames) {
      source_path <- file.path("plots", figures_dir, filename)
      target_path <- file.path(figures_dir, filename)

      if (file.exists(source_path)) {
        file.copy(source_path, target_path, overwrite = TRUE)
        unlink(source_path)
      }
    }
    invisible(NULL)
  }

  trees_beast <- read.nexus(trees_path)
  beast_fmats <- lapply(trees_beast, gen_Fmat, tol = 8)
  true_tree <- read.tree(true_tree_path)
  true_fmat <- gen_Fmat(true_tree, tol = 8)

  cat("Read", length(beast_fmats), "BEAST trees from", trees_path, "\n")
  cat("Read true tree from", true_tree_path, "\n")

  if (is.null(generated_sample_size)) {
    generated_sample_size <- length(beast_fmats)
  }

  fit_result <- fit_and_generate_gibbs_from_fmats(
    F_mats = beast_fmats,
    generated_sample_size = generated_sample_size,
    seed = seed
  )

  generated_fmats <- fit_result$generated_fmats
  beast_average_fmat <- fit_result$beast_mean$raw
  generated_average_fmat <- fit_result$generated_mean$raw

  generated_tree_export <- write_fmats_as_true_labeled_trees(
    F_mats = generated_fmats,
    true_tree = true_tree,
    output_dir = summaries_dir,
    prefix = "generated_gibbs"
  )

  cat("Wrote generated Gibbs trees with true labels to:\n")
  cat("  Newick:", generated_tree_export$newick_path, "\n")
  cat("  NEXUS:", generated_tree_export$nexus_path, "\n")

  fmat_res_true <- fmat_posterior_comparison_plot(
    chain_fmats = generated_fmats,
    beast_fmats = beast_fmats,
    true_fmat = true_fmat,
    threshold = 0.01,
    title = "Gibbs Approximation: F-matrix posterior vs BEAST (true tree reference)",
    save_path = save_path_for_plot_function("fmat_posterior_vs_beast_true_reference.png")
  )

  fmat_res_beast_average <- fmat_posterior_comparison_plot(
    chain_fmats = generated_fmats,
    beast_fmats = beast_fmats,
    true_fmat = beast_average_fmat,
    threshold = 0.01,
    title = "Gibbs Approximation: F-matrix posterior vs BEAST (BEAST average reference)",
    save_path = save_path_for_plot_function("fmat_posterior_vs_beast_average_reference.png")
  )

  plot_subsample_size <- min(
    plot_subsample_cap,
    length(beast_fmats),
    length(generated_fmats)
  )
  set.seed(seed + 10)
  beast_idx <- sample(length(beast_fmats), plot_subsample_size)
  set.seed(seed + 11)
  generated_idx <- sample(length(generated_fmats), plot_subsample_size)

  mds_res <- tree_MDS_comparison_plot(
    trees_gen = generated_fmats[generated_idx],
    trees_data = beast_fmats[beast_idx],
    true_tree = true_fmat,
    M_tree = beast_average_fmat,
    reference_trees = list(generated_average_fmat),
    reference_labels = c("Generated Avg"),
    true_tree_label = "True Tree",
    M_tree_label = "BEAST Avg",
    title = "Gibbs Approximation: tree MDS vs BEAST",
    save_path = save_path_for_plot_function("tree_mds_vs_beast.png")
  )

  hist_res_true <- tree_histogram_comparison_plot(
    chain_gen = generated_fmats,
    chain_data = beast_fmats,
    M_true_fmat = true_fmat,
    title = "Gibbs Approximation: distance to true tree",
    save_path = save_path_for_plot_function("tree_histogram_vs_true.png")
  )

  hist_res_beast_average <- tree_histogram_comparison_plot(
    chain_gen = generated_fmats,
    chain_data = beast_fmats,
    M_true_fmat = beast_average_fmat,
    title = "Gibbs Approximation: distance to BEAST average",
    save_path = save_path_for_plot_function("tree_histogram_vs_beast_average.png")
  )

  move_helper_plot_outputs(c(
    "fmat_posterior_vs_beast_true_reference.png",
    "fmat_posterior_vs_beast_average_reference.png",
    "tree_mds_vs_beast.png",
    "tree_histogram_vs_true.png",
    "tree_histogram_vs_beast_average.png"
  ))

  summary_df <- cbind(
    data.frame(
      experiment = experiment_id,
      beast_tree_file = basename(trees_path),
      true_tree_file = basename(true_tree_path),
      plot_subsample_size = plot_subsample_size
    ),
    fit_result$summary,
    data.frame(
      fmat_l1_distance_vs_beast = fmat_res_true$l1_distance,
      beast_average_l2_to_true = distance_Fmat(beast_average_fmat, true_fmat, dist = "l2"),
      generated_average_l2_to_true = distance_Fmat(generated_average_fmat, true_fmat, dist = "l2"),
      generated_average_l2_to_beast_average = distance_Fmat(
        generated_average_fmat,
        beast_average_fmat,
        dist = "l2"
      ),
      generated_true_labeled_newick = generated_tree_export$newick_path,
      generated_true_labeled_nexus = generated_tree_export$nexus_path
    )
  )

  write.csv(
    summary_df,
    file = file.path(summaries_dir, "experiment_summary.csv"),
    row.names = FALSE
  )
  saveRDS(
    summary_df,
    file = file.path(summaries_dir, "experiment_summary.rds")
  )
  saveRDS(
    fit_result,
    file = file.path(summaries_dir, "fit_result.rds")
  )
  saveRDS(
    fmat_res_true,
    file = file.path(summaries_dir, "fmat_posterior_comparison_true_reference.rds")
  )
  saveRDS(
    fmat_res_beast_average,
    file = file.path(summaries_dir, "fmat_posterior_comparison_beast_average_reference.rds")
  )
  saveRDS(
    mds_res,
    file = file.path(summaries_dir, "tree_mds_vs_beast.rds")
  )
  saveRDS(
    hist_res_true,
    file = file.path(summaries_dir, "tree_histogram_vs_true.rds")
  )
  saveRDS(
    hist_res_beast_average,
    file = file.path(summaries_dir, "tree_histogram_vs_beast_average.rds")
  )
  saveRDS(
    list(
      true_fmat = true_fmat,
      beast_average_fmat = beast_average_fmat,
      generated_average_fmat = generated_average_fmat
    ),
    file = file.path(summaries_dir, "reference_fmats.rds")
  )
  saveRDS(
    fmat_res_true$df,
    file = file.path(summaries_dir, "fmat_posterior_comparison_true_reference_df.rds")
  )
  saveRDS(
    fmat_res_beast_average$df,
    file = file.path(summaries_dir, "fmat_posterior_comparison_beast_average_reference_df.rds")
  )
  saveRDS(
    generated_tree_export$trees,
    file = file.path(summaries_dir, "generated_gibbs_true_labeled_trees.rds")
  )
  write.csv(
    data.frame(
      generated_index = seq_along(fit_result$generated_indices),
      support_index = fit_result$generated_indices
    ),
    file = file.path(summaries_dir, "generated_support_indices.csv"),
    row.names = FALSE
  )

  print(summary_df)
  invisible(fit_result)
}

if (sys.nframe() == 0) {
  run_beast_gibbs_approximation()
}
