## VB Brandon - exact prior sampler version
library(ape)
library(phylodyn)
library(phyclust)
library(phangorn)
library(phylotools)
library(fmatrix)

vb_utils_paths <- c("VB_Utils.R", file.path("Brandon", "gibbs_prior", "VB_Utils.R"))
vb_utils_path <- vb_utils_paths[file.exists(vb_utils_paths)][1]
if (is.na(vb_utils_path)) {
  stop("Could not find VB_Utils.R. Run from the repo root or Brandon/gibbs_prior.")
}
source(vb_utils_path)
set.seed(42)
# ---- Config ----
num_tips            <- 10
rate                <- 0.01
step_size           <- 0.01
num_samps           <- 2000
num_grad_desc_steps <- 1000
num_tip_label_iters <- 100
joint               <- TRUE
init_method         <- "upgma"  # "caterpillar" or "upgma"
log_stats           <- TRUE

# ---- Exact prior support ----
flist_name <- paste0("F.list", num_tips)
if (!exists(flist_name, mode = "list")) {
  stop("Exact F-matrix list not available: ", flist_name)
}
all_Fmats <- get(flist_name)
exact_prior_cache <- phylodyn:::precompute_tree_chain_distance_cache(all_Fmats)

# ---- Generate data + initialize ----
init_results <- phylodyn:::generate_true_M_and_data(num_tips, rate = rate, seq_len = 10000)

inter_coal_times <- coalescent.intervals(init_results$M_true_tree)$interval.length
inter_coal_times[inter_coal_times <= 0.001] <- 0.01
coal_times <- cumsum(inter_coal_times)

if (init_method == "upgma") {
  upgma_tree <- upgma(dist.hamming(init_results$sequences))
  M_est_tree <- phylodyn:::update_time(upgma_tree, coal_times)
} else {
  init_Fmat  <- phylodyn:::gen_caterpillar(num_tips)
  M_est_tree <- mytree_from_F(init_Fmat, coal_times)
}
M_est <- round(phylodyn:::gen_Fmat(M_est_tree, tol = 8), 0)
M_init     <- M_est
M_true     <- round(phylodyn:::gen_Fmat(init_results$M_true_tree, tol = 8), 0)
g_est      <- log(0.01)

# ---- Per-iteration trace logging ----
stats_trace <- list()
stats_csv <- NULL
stats_rds <- NULL
config_rds <- NULL

if (log_stats) {
  log_dir <- file.path("Brandon", "gibbs_prior", "logs", "VB_master_exact_prior")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  stats_csv <- file.path(log_dir, paste0("trace_", run_id, ".csv"))
  stats_rds <- file.path(log_dir, paste0("trace_", run_id, ".rds"))
  config_rds <- file.path(log_dir, paste0("config_", run_id, ".rds"))

  saveRDS(
    list(
      num_tips = num_tips,
      rate = rate,
      step_size = step_size,
      num_samps = num_samps,
      num_grad_desc_steps = num_grad_desc_steps,
      num_tip_label_iters = num_tip_label_iters,
      joint = joint,
      init_method = init_method,
      M_init = M_init,
      M_true = M_true,
      coal_times = coal_times
    ),
    file = config_rds
  )

  cat("Writing VB trace to", stats_csv, "\n")
}

lower_values <- function(M) {
  M[lower.tri(M, diag = TRUE)]
}

matrix_frob_norm <- function(M) {
  sqrt(sum(M * M))
}

unique_fmat_count <- function(F_mats) {
  keys <- vapply(
    F_mats,
    function(Fmat) paste(Fmat[lower.tri(Fmat, diag = TRUE)], collapse = "|"),
    character(1)
  )
  length(unique(keys))
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
                                   phase,
                                   elbo,
                                   log_Z,
                                   M_before,
                                   M_after,
                                   g_before,
                                   g_after,
                                   M_grad,
                                   g_grad,
                                   M_step,
                                   g_step,
                                   samples,
                                   sample_info,
                                   cache_gibbs,
                                   elapsed_seconds) {
  lower_before <- lower_values(M_before)
  lower_after <- lower_values(M_after)
  sample_d2 <- phylodyn:::tree_chain_squared_distances(cache_gibbs, M_before)

  row <- data.frame(
    iteration = iteration,
    phase = phase,
    elbo = elbo,
    elbo_ma20 = mean(tail(elbos2, 20)),
    log_Z = log_Z,
    log_beta_before = g_before,
    log_beta_after = g_after,
    beta_before = exp(g_before),
    beta_after = exp(g_after),
    g_grad = g_grad,
    g_step = g_step,
    M_grad_frob = matrix_frob_norm(M_grad),
    M_grad_mean_abs = mean(abs(M_grad)),
    M_grad_max_abs = max(abs(M_grad)),
    M_step_frob = matrix_frob_norm(M_step),
    M_step_mean_abs = mean(abs(M_step)),
    M_step_max_abs = max(abs(M_step)),
    M_lower_mean_before = mean(lower_before),
    M_lower_mean_after = mean(lower_after),
    M_lower_min_after = min(lower_after),
    M_lower_max_after = max(lower_after),
    M_l2_true_before = distance_Fmat(M_before, M_true, dist = "l2"),
    M_l2_true_after = distance_Fmat(M_after, M_true, dist = "l2"),
    M_l2_init_after = distance_Fmat(M_after, M_init, dist = "l2"),
    sample_unique_fmats = unique_fmat_count(samples),
    sample_mean_d2_to_M_before = mean(sample_d2),
    sample_sd_d2_to_M_before = sd(sample_d2),
    exact_prior_prob_sum = if (!is.null(sample_info$prob_sum)) sample_info$prob_sum else NA_real_,
    exact_prior_ess = if (!is.null(sample_info$ess)) sample_info$ess else NA_real_,
    elapsed_seconds = elapsed_seconds
  )

  append_stats_row(row)
}

# ---- AdaGrad-style gradient descent ----
M_grads <- list(matrix(step_size, nrow = nrow(M_est), ncol = nrow(M_est)))
g_grads <- list(step_size)
elbos2  <- c()

start_time <- Sys.time()
for (i in 1:num_grad_desc_steps) {
  print(paste0("Gradient Descent Iteration ", i))
  print(M_est)
  print(paste0("Estimated log beta: ", g_est))

  M_before_iter <- M_est
  g_before_iter <- g_est

  sample_info <- sample_from_prior(
    M = M_est,
    b = exp(g_est),
    num_samps = num_samps,
    all_Fmats = all_Fmats,
    return_info = TRUE
  )
  samples <- sample_info$samples
  log_Z <- phylodyn:::compute_log_Z_est(beta = exp(g_est), M = M_est, cache = exact_prior_cache)

  if (!joint) {
    cache_gibbs <- phylodyn:::precompute_tree_chain_distance_cache(samples)
    g_up <- phylodyn:::estimate_grad_g2(
      data                = init_results$sequences,
      coal_times          = coal_times,
      samples_gibbs       = samples,
      cache_gibbs         = cache_gibbs,
      M_est               = M_est,
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
    g_step <- step_size / sqrt(Reduce('+', g_grads)) * g_grad
    g_est <- g_est + g_step

    sample_info <- sample_from_prior(
      M = M_est,
      b = exp(g_est),
      num_samps = num_samps,
      all_Fmats = all_Fmats,
      return_info = TRUE
    )
    samples <- sample_info$samples
    log_Z <- phylodyn:::compute_log_Z_est(beta = exp(g_est), M = M_est, cache = exact_prior_cache)

    M_up <- phylodyn:::estimate_grad_M2(
      data                = init_results$sequences,
      coal_times          = coal_times,
      samples_gibbs       = samples,
      M_est               = M_est,
      b_est               = exp(g_est),
      log_Z               = log_Z,
      num_tip_label_iters = num_tip_label_iters,
      rate                = rate
    )
    n_grad_samples <- length(samples)
    M_up$grad <- M_up$grad / n_grad_samples
    M_up$elbo <- M_up$elbo / n_grad_samples
    M_grad <- M_up$grad
    M_grads[[length(M_grads) + 1]] <- M_grad^2
    M_step <- step_size / sqrt(Reduce('+', M_grads)) * M_grad
    M_est <- M_est + M_step
    elbos2 <- c(elbos2, M_up$elbo)
    record_iteration_stats(
      iteration = i,
      phase = "separate",
      elbo = M_up$elbo,
      log_Z = log_Z,
      M_before = M_before_iter,
      M_after = M_est,
      g_before = g_before_iter,
      g_after = g_est,
      M_grad = M_grad,
      g_grad = g_grad,
      M_step = M_step,
      g_step = g_step,
      samples = samples,
      sample_info = sample_info,
      cache_gibbs = phylodyn:::precompute_tree_chain_distance_cache(samples),
      elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )
    print(paste0("ELBO: ", M_up$elbo))
    print(paste0("L2 distance to true M: ", distance_Fmat(M_est, M_true, dist = "l2")))
  } else {
    cache_gibbs <- phylodyn:::precompute_tree_chain_distance_cache(samples)
    J_up <- phylodyn:::estimate_grad_M_g(
      data                = init_results$sequences,
      coal_times          = coal_times,
      samples_gibbs       = samples,
      cache_gibbs         = cache_gibbs,
      M_est               = M_est,
      b_est               = exp(g_est),
      log_Z               = log_Z,
      num_tip_label_iters = num_tip_label_iters,
      rate                = rate
    )
    n_grad_samples <- length(samples)
    J_up$grad_M <- J_up$grad_M / n_grad_samples
    J_up$grad_g <- J_up$grad_g / n_grad_samples
    J_up$elbo <- J_up$elbo / n_grad_samples
    M_grad <- J_up$grad_M
    M_grads[[length(M_grads) + 1]] <- M_grad^2
    M_step <- step_size / sqrt(Reduce('+', M_grads)) * M_grad
    M_est <- M_est + M_step

    g_grad <- J_up$grad_g
    g_grads[[length(g_grads) + 1]] <- g_grad^2
    g_step <- step_size / sqrt(Reduce('+', g_grads)) * g_grad
    g_est <- g_est + g_step

    elbos2 <- c(elbos2, J_up$elbo)
    record_iteration_stats(
      iteration = i,
      phase = "joint",
      elbo = J_up$elbo,
      log_Z = log_Z,
      M_before = M_before_iter,
      M_after = M_est,
      g_before = g_before_iter,
      g_after = g_est,
      M_grad = M_grad,
      g_grad = g_grad,
      M_step = M_step,
      g_step = g_step,
      samples = samples,
      sample_info = sample_info,
      cache_gibbs = cache_gibbs,
      elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )
    print(paste0("ELBO: ", J_up$elbo))
    print(paste0("L2 distance to true M: ", distance_Fmat(M_est, M_true, dist = "l2")))
  }
}
end_time <- Sys.time()
print(end_time - start_time)

if (log_stats) {
  saveRDS(stats_trace, file = stats_rds)
  print(paste0("Saved VB trace CSV: ", stats_csv))
  print(paste0("Saved VB trace RDS: ", stats_rds))
  print(paste0("Saved VB config RDS: ", config_rds))
}

plot(elbos2)

M_estimated_tree <- mytree_from_F(nearby_Fmat(M_est), coal_times)
plot(M_estimated_tree)
plot(init_results$M_true_tree)

print(paste0("L2 distance to true M:    ", distance_Fmat(M_est, M_true, dist = "l2")))
print(paste0("L2 distance to initial M: ", distance_Fmat(M_est, M_init, dist = "l2")))

print(paste0("Estimated L2 distance to true M:    ", distance_Fmat(nearby_Fmat(M_est), M_true, dist = "l2")))
print(paste0("Estimated L2 distance to initial M: ", distance_Fmat(nearby_Fmat(M_est), M_init, dist = "l2")))
print(paste0("L2 distance initial M to True M: ", distance_Fmat(M_true, M_init, dist = "l2")))
