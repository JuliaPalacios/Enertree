## Oracle-style likelihood diagnostic for unlabeled F-matrix VB.
##
## The default settings mirror VB_master_exact_prior.R.

# ---- Experiment settings ----
num_tips <- 10
rate <- 0.01
seq_len <- 1000
seed <- 42
num_replicates <- 20

# Each replicate always compares the labeled true tree, labeled UPGMA tree,
# and random-label scores for the true and UPGMA F-matrices.
random_label_R <- 500

# The exact support scan is much more expensive. When run_exact_support_scan is
# TRUE, each replicate scores every exact 10-tip F-matrix with scan_R random
# labelings, then re-scores the top candidates plus true/UPGMA with refine_R.
run_exact_support_scan <- TRUE
scan_R <- 10
refine_R <- 500
refine_top <- 100

# Use "simulate" to generate a fresh dataset from seed, or "saved" to read
# true_tree.newick and sequences.fasta from disk.
data_mode <- "simulate"
tree_path <- file.path("Brandon", "gibbs_prior", "true_tree.newick")
fasta_path <- file.path("Brandon", "gibbs_prior", "sequences.fasta")

output_dir <- file.path(
  "Brandon", "gibbs_prior", "logs", "oracle_labeled_likelihood_diagnostic"
)

suppressPackageStartupMessages({
  library(ape)
  library(phylodyn)
  library(phyclust)
  library(phangorn)
  library(phylotools)
  library(fmatrix)
})

logmeanexp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

fmat_key <- function(Fmat) {
  paste(Fmat[lower.tri(Fmat, diag = TRUE)], collapse = "|")
}

pml_loglik <- function(tree, sequences, rate) {
  phangorn::pml(
    tree,
    sequences,
    bf = rep(0.25, 4),
    Q = c(1, 2, 1, 1, 2, 1),
    rate = rate
  )$log
}

score_fmat_random_labels <- function(Fmat, coal_times, sequences, rate, R) {
  lls <- numeric(R)
  best_tree <- NULL

  for (r in seq_len(R)) {
    tree <- mytree_from_F(Fmat, coal_times)
    lls[r] <- pml_loglik(tree, sequences, rate)

    if (is.null(best_tree) || lls[r] == max(lls[seq_len(r)])) {
      best_tree <- tree
    }
  }

  list(
    logmeanexp = logmeanexp(lls),
    max = max(lls),
    meanlog = mean(lls),
    sd = sd(lls),
    min = min(lls),
    best_tree = best_tree
  )
}

score_support <- function(indices, all_Fmats, F_true, F_init, coal_times,
                          sequences, rate, R, label) {
  true_key <- fmat_key(F_true)
  init_key <- fmat_key(F_init)
  out <- vector("list", length(indices))

  for (j in seq_along(indices)) {
    i <- indices[j]
    if (j == 1 || j %% 100 == 0 || j == length(indices)) {
      cat(sprintf("[%s] scoring %d / %d (support idx %d)\n",
                  label, j, length(indices), i))
    }

    Fmat <- all_Fmats[[i]]
    score <- score_fmat_random_labels(
      Fmat = Fmat,
      coal_times = coal_times,
      sequences = sequences,
      rate = rate,
      R = R
    )

    out[[j]] <- data.frame(
      idx = i,
      R = R,
      logmeanexp = score$logmeanexp,
      max = score$max,
      meanlog = score$meanlog,
      sd = score$sd,
      min = score$min,
      dist_to_true = distance_Fmat(Fmat, F_true, dist = "l2"),
      dist_to_init = distance_Fmat(Fmat, F_init, dist = "l2"),
      num_cherries = fmatrix:::num_cherries(Fmat),
      is_true_F = identical(fmat_key(Fmat), true_key),
      is_init_F = identical(fmat_key(Fmat), init_key)
    )
  }

  do.call(rbind, out)
}

rank_desc <- function(x) {
  rank(-x, ties.method = "min")
}

cor_or_na <- function(x, y, method = "pearson") {
  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(NA_real_)
  }
  cor(x, y, method = method)
}

load_replicate_data <- function(replicate_id) {
  if (identical(data_mode, "saved")) {
    if (num_replicates != 1) {
      stop("data_mode = 'saved' only supports one replicate.")
    }
    if (!file.exists(tree_path)) {
      stop("Could not find tree_path: ", tree_path)
    }
    if (!file.exists(fasta_path)) {
      stop("Could not find fasta_path: ", fasta_path)
    }

    return(list(
      true_tree = read.tree(tree_path),
      sequences = as.phyDat(read.FASTA(fasta_path))
    ))
  }

  if (!identical(data_mode, "simulate")) {
    stop("data_mode must be either 'simulate' or 'saved'.")
  }

  set.seed(seed + replicate_id - 1)
  init_results <- phylodyn:::generate_true_M_and_data(
    num_tips = num_tips,
    rate = rate,
    seq_len = seq_len
  )

  list(
    true_tree = init_results$M_true_tree,
    sequences = init_results$sequences
  )
}

run_exact_support_diagnostic <- function(replicate_dir, replicate_id, all_Fmats,
                                         F_true, F_init, coal_times, sequences,
                                         true_idx, init_idx) {
  scan_df <- score_support(
    indices = seq_along(all_Fmats),
    all_Fmats = all_Fmats,
    F_true = F_true,
    F_init = F_init,
    coal_times = coal_times,
    sequences = sequences,
    rate = rate,
    R = scan_R,
    label = paste0("replicate ", replicate_id, " scan")
  )

  scan_df$rank_logmeanexp <- rank_desc(scan_df$logmeanexp)
  scan_df$rank_max <- rank_desc(scan_df$max)
  scan_df$rank_meanlog <- rank_desc(scan_df$meanlog)

  top_logmean_idx <- scan_df$idx[order(-scan_df$logmeanexp)][seq_len(min(refine_top, nrow(scan_df)))]
  top_max_idx <- scan_df$idx[order(-scan_df$max)][seq_len(min(refine_top, nrow(scan_df)))]
  refine_indices <- sort(unique(c(top_logmean_idx, top_max_idx, true_idx, init_idx)))
  refine_indices <- refine_indices[!is.na(refine_indices)]

  refine_df <- score_support(
    indices = refine_indices,
    all_Fmats = all_Fmats,
    F_true = F_true,
    F_init = F_init,
    coal_times = coal_times,
    sequences = sequences,
    rate = rate,
    R = refine_R,
    label = paste0("replicate ", replicate_id, " refine")
  )

  refine_df$rank_logmeanexp_within_refined <- rank_desc(refine_df$logmeanexp)
  refine_df$rank_max_within_refined <- rank_desc(refine_df$max)
  refine_df$rank_meanlog_within_refined <- rank_desc(refine_df$meanlog)

  write.csv(scan_df, file.path(replicate_dir, "support_scan.csv"), row.names = FALSE)
  write.csv(refine_df, file.path(replicate_dir, "support_refined.csv"), row.names = FALSE)
  write.csv(head(scan_df[order(-scan_df$logmeanexp), ], 50),
            file.path(replicate_dir, "top_scan_logmeanexp.csv"), row.names = FALSE)
  write.csv(head(scan_df[order(-scan_df$max), ], 50),
            file.path(replicate_dir, "top_scan_max.csv"), row.names = FALSE)
  write.csv(head(refine_df[order(-refine_df$logmeanexp), ], 50),
            file.path(replicate_dir, "top_refined_logmeanexp.csv"), row.names = FALSE)
  write.csv(head(refine_df[order(-refine_df$max), ], 50),
            file.path(replicate_dir, "top_refined_max.csv"), row.names = FALSE)

  true_random_scan <- scan_df[scan_df$is_true_F, , drop = FALSE]
  init_random_scan <- scan_df[scan_df$is_init_F, , drop = FALSE]
  true_random_refine <- refine_df[refine_df$is_true_F, , drop = FALSE]
  init_random_refine <- refine_df[refine_df$is_init_F, , drop = FALSE]

  list(
    scan = scan_df,
    refined = refine_df,
    summary = data.frame(
      refined_size = nrow(refine_df),
      true_scan_rank_logmeanexp = true_random_scan$rank_logmeanexp,
      true_scan_rank_max = true_random_scan$rank_max,
      init_scan_rank_logmeanexp = init_random_scan$rank_logmeanexp,
      init_scan_rank_max = init_random_scan$rank_max,
      true_refine_rank_logmeanexp = true_random_refine$rank_logmeanexp_within_refined,
      true_refine_rank_max = true_random_refine$rank_max_within_refined,
      init_refine_rank_logmeanexp = init_random_refine$rank_logmeanexp_within_refined,
      init_refine_rank_max = init_random_refine$rank_max_within_refined,
      scan_hit_true_label_prob = NA_real_,
      refine_hit_true_label_prob = NA_real_,
      scan_cor_logmeanexp_neg_dist = cor_or_na(scan_df$logmeanexp, -scan_df$dist_to_true),
      scan_spearman_logmeanexp_neg_dist = cor_or_na(scan_df$logmeanexp, -scan_df$dist_to_true, method = "spearman"),
      scan_cor_max_neg_dist = cor_or_na(scan_df$max, -scan_df$dist_to_true),
      scan_spearman_max_neg_dist = cor_or_na(scan_df$max, -scan_df$dist_to_true, method = "spearman"),
      refine_cor_logmeanexp_neg_dist = cor_or_na(refine_df$logmeanexp, -refine_df$dist_to_true),
      refine_spearman_logmeanexp_neg_dist = cor_or_na(refine_df$logmeanexp, -refine_df$dist_to_true, method = "spearman"),
      refine_cor_max_neg_dist = cor_or_na(refine_df$max, -refine_df$dist_to_true),
      refine_spearman_max_neg_dist = cor_or_na(refine_df$max, -refine_df$dist_to_true, method = "spearman"),
      top_scan_logmeanexp_dist_mean = mean(head(scan_df[order(-scan_df$logmeanexp), "dist_to_true"], 20)),
      top_scan_max_dist_mean = mean(head(scan_df[order(-scan_df$max), "dist_to_true"], 20)),
      top_refine_logmeanexp_dist_mean = mean(head(refine_df[order(-refine_df$logmeanexp), "dist_to_true"], 20)),
      top_refine_max_dist_mean = mean(head(refine_df[order(-refine_df$max), "dist_to_true"], 20))
    )
  )
}

run_one_replicate <- function(replicate_id, run_dir, all_Fmats, support_keys) {
  replicate_seed <- seed + replicate_id - 1
  cat("\nReplicate", replicate_id, "of", num_replicates,
      "(seed", replicate_seed, ")\n")

  replicate_data <- load_replicate_data(replicate_id)
  true_tree <- replicate_data$true_tree
  sequences <- replicate_data$sequences

  inter_coal_times <- coalescent.intervals(true_tree)$interval.length
  inter_coal_times[inter_coal_times <= 0.001] <- 0.01
  coal_times <- cumsum(inter_coal_times)

  upgma_tree <- upgma(dist.hamming(sequences))
  upgma_timed_tree <- phylodyn:::update_time(upgma_tree, coal_times)

  F_true <- round(phylodyn:::gen_Fmat(true_tree, tol = 8), 0)
  F_init <- round(phylodyn:::gen_Fmat(upgma_timed_tree, tol = 8), 0)

  true_idx <- match(fmat_key(F_true), support_keys)
  init_idx <- match(fmat_key(F_init), support_keys)

  if (is.na(true_idx)) {
    stop("True F matrix was not found in F.list", num_tips,
         " for replicate ", replicate_id)
  }
  if (is.na(init_idx)) {
    warning("UPGMA init F matrix was not found in F.list", num_tips,
            " for replicate ", replicate_id)
  }

  true_labeled_ll <- pml_loglik(true_tree, sequences, rate)
  upgma_labeled_ll <- pml_loglik(upgma_timed_tree, sequences, rate)

  true_random <- score_fmat_random_labels(
    Fmat = F_true,
    coal_times = coal_times,
    sequences = sequences,
    rate = rate,
    R = random_label_R
  )
  init_random <- score_fmat_random_labels(
    Fmat = F_init,
    coal_times = coal_times,
    sequences = sequences,
    rate = rate,
    R = random_label_R
  )

  true_cherries <- fmatrix:::num_cherries(F_true)
  single_draw_true_label_prob <- 2^true_cherries / factorial(num_tips)

  anchor_df <- data.frame(
    replicate = replicate_id,
    quantity = c(
      "true_labeled_tree",
      "upgma_labeled_tree",
      "true_F_random_labels_logmeanexp",
      "true_F_random_labels_max",
      "true_F_random_labels_meanlog",
      "init_F_random_labels_logmeanexp",
      "init_F_random_labels_max",
      "init_F_random_labels_meanlog"
    ),
    log_likelihood = c(
      true_labeled_ll,
      upgma_labeled_ll,
      true_random$logmeanexp,
      true_random$max,
      true_random$meanlog,
      init_random$logmeanexp,
      init_random$max,
      init_random$meanlog
    )
  )

  support_summary <- data.frame(
    refined_size = NA_integer_,
    true_scan_rank_logmeanexp = NA_real_,
    true_scan_rank_max = NA_real_,
    init_scan_rank_logmeanexp = NA_real_,
    init_scan_rank_max = NA_real_,
    true_refine_rank_logmeanexp = NA_real_,
    true_refine_rank_max = NA_real_,
    init_refine_rank_logmeanexp = NA_real_,
    init_refine_rank_max = NA_real_,
    scan_hit_true_label_prob = NA_real_,
    refine_hit_true_label_prob = NA_real_,
    scan_cor_logmeanexp_neg_dist = NA_real_,
    scan_spearman_logmeanexp_neg_dist = NA_real_,
    scan_cor_max_neg_dist = NA_real_,
    scan_spearman_max_neg_dist = NA_real_,
    refine_cor_logmeanexp_neg_dist = NA_real_,
    refine_spearman_logmeanexp_neg_dist = NA_real_,
    refine_cor_max_neg_dist = NA_real_,
    refine_spearman_max_neg_dist = NA_real_,
    top_scan_logmeanexp_dist_mean = NA_real_,
    top_scan_max_dist_mean = NA_real_,
    top_refine_logmeanexp_dist_mean = NA_real_,
    top_refine_max_dist_mean = NA_real_
  )

  scan_result <- NULL
  replicate_dir <- file.path(run_dir, sprintf("replicate_%03d", replicate_id))

  if (run_exact_support_scan) {
    dir.create(replicate_dir, recursive = TRUE, showWarnings = FALSE)
    scan_result <- run_exact_support_diagnostic(
      replicate_dir = replicate_dir,
      replicate_id = replicate_id,
      all_Fmats = all_Fmats,
      F_true = F_true,
      F_init = F_init,
      coal_times = coal_times,
      sequences = sequences,
      true_idx = true_idx,
      init_idx = init_idx
    )
    support_summary <- scan_result$summary
    support_summary$scan_hit_true_label_prob <- 1 - (1 - single_draw_true_label_prob)^scan_R
    support_summary$refine_hit_true_label_prob <- 1 - (1 - single_draw_true_label_prob)^refine_R

    saveRDS(
      list(
        replicate = replicate_id,
        scan = scan_result$scan,
        refined = scan_result$refined,
        F_true = F_true,
        F_init = F_init,
        true_tree = true_tree,
        upgma_timed_tree = upgma_timed_tree,
        coal_times = coal_times
      ),
      file.path(replicate_dir, "diagnostic_results.rds")
    )
  }

  summary_df <- data.frame(
    replicate = replicate_id,
    replicate_seed = replicate_seed,
    data_mode = data_mode,
    num_tips = num_tips,
    rate = rate,
    seq_len = seq_len,
    random_label_R = random_label_R,
    run_exact_support_scan = run_exact_support_scan,
    support_size = length(all_Fmats),
    scan_R = if (run_exact_support_scan) scan_R else NA_integer_,
    refine_R = if (run_exact_support_scan) refine_R else NA_integer_,
    refine_top = if (run_exact_support_scan) refine_top else NA_integer_,
    true_idx = true_idx,
    init_idx = init_idx,
    init_l2_to_true = distance_Fmat(F_init, F_true, dist = "l2"),
    true_labeled_ll = true_labeled_ll,
    upgma_labeled_ll = upgma_labeled_ll,
    upgma_minus_true_labeled_ll = upgma_labeled_ll - true_labeled_ll,
    true_minus_upgma_labeled_ll = true_labeled_ll - upgma_labeled_ll,
    upgma_beats_true_labeled = upgma_labeled_ll > true_labeled_ll,
    true_random_logmeanexp = true_random$logmeanexp,
    true_random_max = true_random$max,
    true_random_meanlog = true_random$meanlog,
    true_random_sd = true_random$sd,
    init_random_logmeanexp = init_random$logmeanexp,
    init_random_max = init_random$max,
    init_random_meanlog = init_random$meanlog,
    init_random_sd = init_random$sd,
    true_random_logmeanexp_minus_init = true_random$logmeanexp - init_random$logmeanexp,
    true_random_max_minus_init = true_random$max - init_random$max,
    true_labeled_minus_true_random_logmeanexp = true_labeled_ll - true_random$logmeanexp,
    true_num_cherries = true_cherries,
    single_draw_true_label_prob = single_draw_true_label_prob,
    random_label_hit_true_label_prob = 1 - (1 - single_draw_true_label_prob)^random_label_R
  )

  summary_df <- cbind(summary_df, support_summary)

  cat("true labeled LL:", round(true_labeled_ll, 3), "\n")
  cat("UPGMA labeled LL:", round(upgma_labeled_ll, 3), "\n")
  cat("UPGMA - true:", round(upgma_labeled_ll - true_labeled_ll, 3), "\n")
  cat("true random-label logmeanexp:", round(true_random$logmeanexp, 3), "\n")
  cat("init random-label logmeanexp:", round(init_random$logmeanexp, 3), "\n")

  list(
    summary = summary_df,
    anchors = anchor_df
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(output_dir, run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

flist_name <- paste0("F.list", num_tips)
if (!exists(flist_name, mode = "list")) {
  stop("Exact F-matrix list not available: ", flist_name)
}
all_Fmats <- get(flist_name)

support_dim <- nrow(all_Fmats[[1]])
if (support_dim != num_tips - 1) {
  stop("Support dimension mismatch: expected ", num_tips - 1, ", got ", support_dim)
}

support_keys <- vapply(all_Fmats, fmat_key, character(1))

cat("Run directory:", run_dir, "\n")
cat("Data mode:", data_mode, "\n")
cat("num_tips:", num_tips, "rate:", rate, "seq_len:", seq_len,
    "seed:", seed, "num_replicates:", num_replicates, "\n")
cat("random_label_R:", random_label_R, "\n")
cat("run_exact_support_scan:", run_exact_support_scan, "\n")
if (run_exact_support_scan) {
  cat("support size:", length(all_Fmats), "scan_R:", scan_R,
      "refine_R:", refine_R, "refine_top:", refine_top, "\n")
}

replicate_results <- lapply(
  seq_len(num_replicates),
  run_one_replicate,
  run_dir = run_dir,
  all_Fmats = all_Fmats,
  support_keys = support_keys
)

summary_df <- do.call(rbind, lapply(replicate_results, `[[`, "summary"))
anchor_df <- do.call(rbind, lapply(replicate_results, `[[`, "anchors"))

overall_summary <- data.frame(
  num_replicates = num_replicates,
  data_mode = data_mode,
  num_tips = num_tips,
  rate = rate,
  seq_len = seq_len,
  random_label_R = random_label_R,
  true_beats_upgma_count = sum(!summary_df$upgma_beats_true_labeled),
  upgma_beats_true_count = sum(summary_df$upgma_beats_true_labeled),
  upgma_beats_true_fraction = mean(summary_df$upgma_beats_true_labeled),
  mean_upgma_minus_true_labeled_ll = mean(summary_df$upgma_minus_true_labeled_ll),
  median_upgma_minus_true_labeled_ll = median(summary_df$upgma_minus_true_labeled_ll),
  sd_upgma_minus_true_labeled_ll = sd(summary_df$upgma_minus_true_labeled_ll),
  mean_init_l2_to_true = mean(summary_df$init_l2_to_true),
  mean_true_random_logmeanexp_minus_init = mean(summary_df$true_random_logmeanexp_minus_init),
  median_true_random_logmeanexp_minus_init = median(summary_df$true_random_logmeanexp_minus_init),
  mean_true_labeled_minus_true_random_logmeanexp =
    mean(summary_df$true_labeled_minus_true_random_logmeanexp)
)

write.csv(summary_df, file.path(run_dir, "summary.csv"), row.names = FALSE)
write.csv(summary_df, file.path(run_dir, "summary_by_replicate.csv"), row.names = FALSE)
write.csv(anchor_df, file.path(run_dir, "anchor_likelihoods.csv"), row.names = FALSE)
write.csv(anchor_df, file.path(run_dir, "anchor_likelihoods_by_replicate.csv"), row.names = FALSE)
write.csv(overall_summary, file.path(run_dir, "overall_summary.csv"), row.names = FALSE)

saveRDS(
  list(
    settings = list(
      num_tips = num_tips,
      rate = rate,
      seq_len = seq_len,
      seed = seed,
      num_replicates = num_replicates,
      random_label_R = random_label_R,
      run_exact_support_scan = run_exact_support_scan,
      scan_R = scan_R,
      refine_R = refine_R,
      refine_top = refine_top,
      data_mode = data_mode
    ),
    overall_summary = overall_summary,
    summary = summary_df,
    anchors = anchor_df
  ),
  file.path(run_dir, "diagnostic_results.rds")
)

cat("\nOverall summary:\n")
print(overall_summary, row.names = FALSE)

cat("\nReplicate summary:\n")
print(summary_df, row.names = FALSE)
