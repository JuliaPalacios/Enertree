## VB Brandon - exact prior sampler with fixed M.
##
## This is based on Brandon/gibbs_prior/VB_master_exact_prior.R, but M is
## selected once and held fixed. The gradient loop estimates only beta.

suppressPackageStartupMessages({
  library(ape)
  library(phylodyn)
  library(phyclust)
  library(phangorn)
  library(phylotools)
  library(fmatrix)
})

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }

  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) {
        frame$ofile
      } else {
        NA_character_
      }
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[length(frame_files)], mustWork = FALSE)))
  }

  getwd()
})

vb_utils_paths <- unique(c(
  file.path(script_dir, "VB_Utils.R"),
  "VB_Utils.R",
  file.path("Brandon", "fixed_M_VB", "VB_Utils.R")
))
vb_utils_path <- vb_utils_paths[file.exists(vb_utils_paths)][1]
if (is.na(vb_utils_path)) {
  stop("Could not find VB_Utils.R. Run from the repo root or Brandon/fixed_M_VB.")
}
source(vb_utils_path)

# ---- Config ----
seed                <- 42L
num_tips            <- 10
rate                <- 0.01
sequence_length     <- 10000
step_size           <- 0.01
num_samps           <- 2000
num_grad_desc_steps <- 100
num_tip_label_iters <- 100
initial_log_beta    <- log(0.01)
fixed_M_sources     <- c("true")  # "true", "upgma", or c("true", "upgma")
log_stats           <- TRUE
make_plots          <- TRUE
save_input_data     <- TRUE
generate_final_trees <- TRUE
final_tree_count    <- 100000L
final_tree_chunk_size <- 5000L


set.seed(seed)

fixed_M_sources <- tolower(fixed_M_sources)
valid_fixed_M_sources <- c("true", "upgma")
bad_sources <- setdiff(fixed_M_sources, valid_fixed_M_sources)
if (length(bad_sources) > 0) {
  stop("Unknown fixed_M_sources value(s): ", paste(bad_sources, collapse = ", "))
}

base_log_dir <- file.path(script_dir, "logs", "VB_fixed_M_exact_prior")
input_data_dir <- file.path(base_log_dir, fixed_M_sources[[1]])
input_data_paths <- list(
  sequences_fasta = file.path(input_data_dir, "sequences.fasta"),
  true_tree_newick = file.path(input_data_dir, "true_tree.newick")
)
if (save_input_data) {
  dir.create(input_data_dir, recursive = TRUE, showWarnings = FALSE)
}

# ---- Exact prior support ----
flist_name <- paste0("F.list", num_tips)
if (!exists(flist_name, mode = "list")) {
  stop("Exact F-matrix list not available: ", flist_name)
}
all_Fmats <- get(flist_name)
exact_prior_cache <- phylodyn:::precompute_tree_chain_distance_cache(all_Fmats)

# ---- Generate data and fixed-M candidates ----
init_results <- phylodyn:::generate_true_M_and_data(
  num_tips,
  rate = rate,
  seq_len = sequence_length,
  write_files = save_input_data,
  name_fasta = input_data_paths$sequences_fasta,
  name_tree = input_data_paths$true_tree_newick
)

inter_coal_times <- coalescent.intervals(init_results$M_true_tree)$interval.length
inter_coal_times[inter_coal_times <= 0.001] <- 0.01
coal_times <- cumsum(inter_coal_times)

M_true_tree <- init_results$M_true_tree
M_true <- round(phylodyn:::gen_Fmat(M_true_tree, tol = 8), 0)

upgma_tree <- upgma(dist.hamming(init_results$sequences))
upgma_tree <- phylodyn:::update_time(upgma_tree, coal_times)
M_upgma <- round(phylodyn:::gen_Fmat(upgma_tree, tol = 8), 0)

fixed_M_candidates <- list(
  true = list(M = M_true, tree = M_true_tree),
  upgma = list(M = M_upgma, tree = upgma_tree)
)

# ---- Per-iteration trace logging helpers ----
lower_values <- function(M) {
  M[lower.tri(M, diag = TRUE)]
}

unique_fmat_count <- function(F_mats) {
  keys <- vapply(
    F_mats,
    function(Fmat) paste(Fmat[lower.tri(Fmat, diag = TRUE)], collapse = "|"),
    character(1)
  )
  length(unique(keys))
}

write_run_plots <- function(log_dir, run_id, fixed_M_source, fixed_tree, elbos2) {
  if (!make_plots) {
    return(invisible(NULL))
  }

  if (length(elbos2) > 0) {
    png(file.path(log_dir, paste0("elbo_", run_id, ".png")))
    plot(
      elbos2,
      type = "l",
      xlab = "Gradient descent iteration",
      ylab = "ELBO",
      main = paste("Fixed-M VB ELBO:", fixed_M_source)
    )
    dev.off()
  }

  png(file.path(log_dir, paste0("fixed_tree_", run_id, ".png")))
  plot(fixed_tree, main = paste("Fixed M tree:", fixed_M_source))
  dev.off()

  png(file.path(log_dir, paste0("true_tree_", run_id, ".png")))
  plot(M_true_tree, main = "True M tree")
  dev.off()
}

format_count <- function(x) {
  format(x, scientific = FALSE, trim = TRUE)
}

write_run_inputs <- function(log_dir, run_id, fixed_M_source, M_fixed, fixed_tree) {
  if (!save_input_data) {
    return(NULL)
  }

  prefix <- paste0(fixed_M_source, "_", run_id)
  paths <- list(
    sequences_fasta = file.path(log_dir, "sequences.fasta"),
    true_tree_newick = file.path(log_dir, "true_tree.newick"),
    upgma_tree_newick = file.path(log_dir, paste0("upgma_tree_", prefix, ".newick")),
    fixed_tree_newick = file.path(log_dir, paste0("fixed_tree_", prefix, ".newick")),
    M_true_csv = file.path(log_dir, paste0("M_true_", prefix, ".csv")),
    M_upgma_csv = file.path(log_dir, paste0("M_upgma_", prefix, ".csv")),
    M_fixed_csv = file.path(log_dir, paste0("M_fixed_", prefix, ".csv")),
    data_rds = file.path(log_dir, paste0("input_data_", prefix, ".rds"))
  )

  if (!identical(normalizePath(input_data_paths$sequences_fasta, mustWork = FALSE),
                 normalizePath(paths$sequences_fasta, mustWork = FALSE))) {
    file.copy(input_data_paths$sequences_fasta, paths$sequences_fasta, overwrite = TRUE)
    file.copy(input_data_paths$true_tree_newick, paths$true_tree_newick, overwrite = TRUE)
  }
  ape::write.tree(upgma_tree, file = paths$upgma_tree_newick)
  ape::write.tree(fixed_tree, file = paths$fixed_tree_newick)
  write.csv(M_true, paths$M_true_csv, row.names = FALSE)
  write.csv(M_upgma, paths$M_upgma_csv, row.names = FALSE)
  write.csv(M_fixed, paths$M_fixed_csv, row.names = FALSE)
  saveRDS(
    list(
      seed = seed,
      num_tips = num_tips,
      rate = rate,
      sequence_length = sequence_length,
      fixed_M_source = fixed_M_source,
      sequences = init_results$sequences,
      M_true_tree = M_true_tree,
      upgma_tree = upgma_tree,
      fixed_tree = fixed_tree,
      M_true = M_true,
      M_upgma = M_upgma,
      M_fixed = M_fixed,
      coal_times = coal_times,
      input_paths = paths
    ),
    file = paths$data_rds
  )

  paths
}

write_final_tree_sample <- function(log_dir, run_id, fixed_M_source, run_index, M_fixed, final_beta) {
  if (!generate_final_trees) {
    return(NULL)
  }

  tree_count_label <- format_count(final_tree_count)
  prefix <- paste0(fixed_M_source, "_", run_id)
  trees_path <- file.path(
    log_dir,
    paste0("generated_", tree_count_label, "_trees_final_beta_fixed_M_", prefix, ".trees")
  )
  indices_path <- file.path(
    log_dir,
    paste0("generated_", tree_count_label, "_tree_indices_final_beta_fixed_M_", prefix, ".rds")
  )
  summary_path <- file.path(
    log_dir,
    paste0("generated_", tree_count_label, "_trees_final_beta_fixed_M_", prefix, "_summary.rds")
  )
  sample_seed <- seed + 10000L + run_index

  cat("Generating ", tree_count_label, " final trees at beta = ", final_beta, "\n", sep = "")
  info <- write_prior_trees_to_file(
    M = M_fixed,
    b = final_beta,
    num_trees = final_tree_count,
    all_Fmats = all_Fmats,
    coal_times = coal_times,
    file = trees_path,
    cache = exact_prior_cache,
    chunk_size = final_tree_chunk_size,
    seed = sample_seed,
    indices_file = indices_path
  )

  info$sample_seed <- sample_seed
  info$trees_path <- trees_path
  info$indices_path <- indices_path
  info$summary_path <- summary_path
  saveRDS(info, file = summary_path)
  info
}

run_fixed_M_vb <- function(fixed_M_source, run_index) {
  fixed_M_info <- fixed_M_candidates[[fixed_M_source]]
  M_fixed <- fixed_M_info$M
  fixed_tree <- fixed_M_info$tree

  g_est <- initial_log_beta
  g_grads <- list(step_size)
  elbos2 <- c()
  stats_trace <- list()
  stats_csv <- NULL
  stats_rds <- NULL
  config_rds <- NULL
  result_rds <- NULL
  log_dir <- NULL
  run_id <- NULL
  input_paths <- NULL
  final_tree_info <- NULL

  if (log_stats) {
    log_dir <- file.path(script_dir, "logs", "VB_fixed_M_exact_prior", fixed_M_source)
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
    stats_csv <- file.path(log_dir, paste0("trace_", run_id, ".csv"))
    stats_rds <- file.path(log_dir, paste0("trace_", run_id, ".rds"))
    config_rds <- file.path(log_dir, paste0("config_", run_id, ".rds"))
    result_rds <- file.path(log_dir, paste0("result_", run_id, ".rds"))

    saveRDS(
      list(
        seed = seed,
        num_tips = num_tips,
        rate = rate,
        sequence_length = sequence_length,
        step_size = step_size,
        num_samps = num_samps,
        num_grad_desc_steps = num_grad_desc_steps,
        num_tip_label_iters = num_tip_label_iters,
        initial_log_beta = initial_log_beta,
        fixed_M_source = fixed_M_source,
        save_input_data = save_input_data,
        generate_final_trees = generate_final_trees,
        final_tree_count = final_tree_count,
        final_tree_chunk_size = final_tree_chunk_size,
        input_data_paths = input_data_paths,
        M_fixed = M_fixed,
        M_true = M_true,
        M_upgma = M_upgma,
        M_true_tree = M_true_tree,
        upgma_tree = upgma_tree,
        fixed_tree = fixed_tree,
        coal_times = coal_times,
        true_upgma_l2 = distance_Fmat(M_true, M_upgma, dist = "l2")
      ),
      file = config_rds
    )

    cat("Writing fixed-M VB trace to", stats_csv, "\n")
    input_paths <- write_run_inputs(log_dir, run_id, fixed_M_source, M_fixed, fixed_tree)
  }

  append_stats_row <- function(row) {
    if (!log_stats) {
      return(invisible(NULL))
    }

    stats_trace[[length(stats_trace) + 1]] <<- row
    write.table(
      row,
      file = stats_csv,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(stats_csv),
      append = file.exists(stats_csv)
    )
    invisible(NULL)
  }

  record_iteration_stats <- function(iteration,
                                     elbo,
                                     log_Z,
                                     g_before,
                                     g_after,
                                     g_grad,
                                     g_step,
                                     samples,
                                     sample_info,
                                     cache_gibbs,
                                     elapsed_seconds) {
    fixed_lower <- lower_values(M_fixed)
    sample_d2 <- phylodyn:::tree_chain_squared_distances(cache_gibbs, M_fixed)

    row <- data.frame(
      iteration = iteration,
      fixed_M_source = fixed_M_source,
      elbo = elbo,
      elbo_ma20 = mean(tail(elbos2, 20)),
      log_Z = log_Z,
      log_beta_before = g_before,
      log_beta_after = g_after,
      beta_before = exp(g_before),
      beta_after = exp(g_after),
      g_grad = g_grad,
      g_step = g_step,
      M_is_fixed = TRUE,
      M_fixed_lower_mean = mean(fixed_lower),
      M_fixed_lower_min = min(fixed_lower),
      M_fixed_lower_max = max(fixed_lower),
      M_fixed_l2_true = distance_Fmat(M_fixed, M_true, dist = "l2"),
      M_fixed_l2_upgma = distance_Fmat(M_fixed, M_upgma, dist = "l2"),
      true_upgma_l2 = distance_Fmat(M_true, M_upgma, dist = "l2"),
      sample_unique_fmats = unique_fmat_count(samples),
      sample_mean_d2_to_M_fixed = mean(sample_d2),
      sample_sd_d2_to_M_fixed = sd(sample_d2),
      exact_prior_prob_sum = if (!is.null(sample_info$prob_sum)) sample_info$prob_sum else NA_real_,
      exact_prior_ess = if (!is.null(sample_info$ess)) sample_info$ess else NA_real_,
      exact_prior_support_size = if (!is.null(sample_info$support_size)) sample_info$support_size else NA_integer_,
      elapsed_seconds = elapsed_seconds
    )

    append_stats_row(row)
  }

  start_time <- Sys.time()
  for (i in seq_len(num_grad_desc_steps)) {
    print(paste0("Gradient Descent Iteration ", i))
    print(paste0("Fixed M source: ", fixed_M_source))
    print(M_fixed)
    print(paste0("Estimated log beta: ", g_est))

    g_before_iter <- g_est

    sample_info <- sample_from_prior(
      M = M_fixed,
      b = exp(g_est),
      num_samps = num_samps,
      all_Fmats = all_Fmats,
      cache = exact_prior_cache,
      return_info = TRUE
    )
    samples <- sample_info$samples
    log_Z <- sample_info$log_Z
    cache_gibbs <- phylodyn:::precompute_tree_chain_distance_cache(samples)

    g_up <- phylodyn:::estimate_grad_g2(
      data                = init_results$sequences,
      coal_times          = coal_times,
      samples_gibbs       = samples,
      cache_gibbs         = cache_gibbs,
      M_est               = M_fixed,
      b_est               = exp(g_est),
      log_Z               = log_Z,
      num_tip_label_iters = num_tip_label_iters,
      rate                = rate
    )
    n_grad_samples <- length(samples)
    g_up$grad <- g_up$grad / n_grad_samples
    g_up$elbo <- g_up$elbo / n_grad_samples

    g_grad <- g_up$grad
    g_grads[[length(g_grads) + 1]] <- g_grad^2
    g_step <- step_size / sqrt(Reduce("+", g_grads)) * g_grad
    g_est <- g_est + g_step
    elbos2 <- c(elbos2, g_up$elbo)

    record_iteration_stats(
      iteration = i,
      elbo = g_up$elbo,
      log_Z = log_Z,
      g_before = g_before_iter,
      g_after = g_est,
      g_grad = g_grad,
      g_step = g_step,
      samples = samples,
      sample_info = sample_info,
      cache_gibbs = cache_gibbs,
      elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )

    print(paste0("ELBO: ", g_up$elbo))
    print(paste0("L2 distance fixed M to true M: ", distance_Fmat(M_fixed, M_true, dist = "l2")))
    print(paste0("L2 distance fixed M to UPGMA M: ", distance_Fmat(M_fixed, M_upgma, dist = "l2")))
  }
  end_time <- Sys.time()
  print(end_time - start_time)

  if (log_stats) {
    final_tree_info <- write_final_tree_sample(
      log_dir = log_dir,
      run_id = run_id,
      fixed_M_source = fixed_M_source,
      run_index = run_index,
      M_fixed = M_fixed,
      final_beta = exp(g_est)
    )
  } else if (generate_final_trees) {
    warning("generate_final_trees is TRUE, but log_stats is FALSE; no output directory was created.")
  }

  if (log_stats) {
    saveRDS(stats_trace, file = stats_rds)
    saveRDS(
      list(
        seed = seed,
        fixed_M_source = fixed_M_source,
        M_fixed = M_fixed,
        M_true = M_true,
        M_upgma = M_upgma,
        M_true_tree = M_true_tree,
        upgma_tree = upgma_tree,
        fixed_tree = fixed_tree,
        sequences = init_results$sequences,
        coal_times = coal_times,
        final_log_beta = g_est,
        final_beta = exp(g_est),
        elbos = elbos2,
        stats_trace = stats_trace,
        elapsed = end_time - start_time,
        input_paths = input_paths,
        final_tree_info = final_tree_info
      ),
      file = result_rds
    )
    write_run_plots(log_dir, run_id, fixed_M_source, fixed_tree, elbos2)

    print(paste0("Saved fixed-M VB trace CSV: ", stats_csv))
    print(paste0("Saved fixed-M VB trace RDS: ", stats_rds))
    print(paste0("Saved fixed-M VB config RDS: ", config_rds))
    print(paste0("Saved fixed-M VB result RDS: ", result_rds))
    if (!is.null(final_tree_info)) {
      print(paste0("Saved final generated trees: ", final_tree_info$trees_path))
      print(paste0("Saved final generated tree indices: ", final_tree_info$indices_path))
    }
  }

  print(paste0("Final estimated log beta: ", g_est))
  print(paste0("Final estimated beta:     ", exp(g_est)))
  print(paste0("L2 distance fixed M to true M:  ", distance_Fmat(M_fixed, M_true, dist = "l2")))
  print(paste0("L2 distance fixed M to UPGMA M: ", distance_Fmat(M_fixed, M_upgma, dist = "l2")))

  invisible(list(
    fixed_M_source = fixed_M_source,
    M_fixed = M_fixed,
    final_log_beta = g_est,
    final_beta = exp(g_est),
    elbos = elbos2,
    stats_trace = stats_trace
  ))
}

results <- Map(run_fixed_M_vb, fixed_M_sources, seq_along(fixed_M_sources))
names(results) <- fixed_M_sources
