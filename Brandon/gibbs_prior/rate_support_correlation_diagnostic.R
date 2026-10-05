## Rate-sensitivity support scan for unlabeled F-matrix likelihood signal.
##
## For each rate, this script samples candidate F-matrices from the exact
## support, scores each candidate by random labelings, and correlates the score
## with closeness to the true F-matrix. Positive correlations mean higher
## likelihood is associated with smaller F-distance to the truth.

# ---- Experiment settings ----
num_tips <- 10
base_seed <- 42
num_replicates <- 2

experiments <- list(
  list(
    name = "higher_rates_seq1000",
    seq_len = 1000,
    rate_grid = c(0.1, 0.15, 0.2, 0.3)
  ),
  list(
    name = "longer_sequences_moderate_rates",
    seq_len = 10000,
    rate_grid = c(0.01, 0.02, 0.05, 0.1)
  )
)

support_sample_size <- 500
random_label_R <- 50

output_dir <- file.path(
  "Brandon", "gibbs_prior", "logs", "rate_support_correlation_diagnostic"
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

  for (r in seq_len(R)) {
    tree <- mytree_from_F(Fmat, coal_times)
    lls[r] <- pml_loglik(tree, sequences, rate)
  }

  list(
    logmeanexp = logmeanexp(lls),
    max = max(lls),
    meanlog = mean(lls),
    sd = sd(lls)
  )
}

cor_or_na <- function(x, y, method = "pearson") {
  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(NA_real_)
  }
  cor(x, y, method = method)
}

rank_desc <- function(x) {
  rank(-x, ties.method = "min")
}

sample_support_indices <- function(all_indices, true_idx, upgma_idx,
                                   sample_size) {
  anchor_indices <- unique(na.omit(c(true_idx, upgma_idx)))
  sample_pool <- setdiff(all_indices, anchor_indices)
  sampled <- sample(sample_pool, min(sample_size, length(sample_pool)))
  sort(unique(c(sampled, anchor_indices)))
}

score_support_subset <- function(indices, all_Fmats, F_true, F_upgma,
                                 coal_times, sequences, rate, R, label) {
  true_key <- fmat_key(F_true)
  upgma_key <- fmat_key(F_upgma)
  out <- vector("list", length(indices))

  for (j in seq_along(indices)) {
    if (j == 1 || j %% 100 == 0 || j == length(indices)) {
      cat(sprintf("[%s] scoring %d / %d\n", label, j, length(indices)))
    }

    idx <- indices[j]
    Fmat <- all_Fmats[[idx]]
    score <- score_fmat_random_labels(
      Fmat = Fmat,
      coal_times = coal_times,
      sequences = sequences,
      rate = rate,
      R = R
    )

    out[[j]] <- data.frame(
      idx = idx,
      logmeanexp = score$logmeanexp,
      max = score$max,
      meanlog = score$meanlog,
      sd = score$sd,
      dist_to_true = distance_Fmat(Fmat, F_true, dist = "l2"),
      dist_to_upgma = distance_Fmat(Fmat, F_upgma, dist = "l2"),
      is_true_F = identical(fmat_key(Fmat), true_key),
      is_upgma_F = identical(fmat_key(Fmat), upgma_key)
    )
  }

  do.call(rbind, out)
}

pairwise_hamming_mean <- function(sequences) {
  hamming_mat <- as.matrix(dist.hamming(sequences))
  mean(hamming_mat[upper.tri(hamming_mat)])
}

summarize_scores <- function(scores_df) {
  scores_df$rank_logmeanexp <- rank_desc(scores_df$logmeanexp)
  scores_df$rank_max <- rank_desc(scores_df$max)

  true_row <- scores_df[scores_df$is_true_F, , drop = FALSE]
  upgma_row <- scores_df[scores_df$is_upgma_F, , drop = FALSE]

  data.frame(
    support_scored = nrow(scores_df),
    cor_logmeanexp_neg_dist =
      cor_or_na(scores_df$logmeanexp, -scores_df$dist_to_true),
    spearman_logmeanexp_neg_dist =
      cor_or_na(scores_df$logmeanexp, -scores_df$dist_to_true,
                method = "spearman"),
    cor_max_neg_dist = cor_or_na(scores_df$max, -scores_df$dist_to_true),
    spearman_max_neg_dist =
      cor_or_na(scores_df$max, -scores_df$dist_to_true,
                method = "spearman"),
    true_rank_logmeanexp = true_row$rank_logmeanexp,
    true_rank_max = true_row$rank_max,
    upgma_rank_logmeanexp = upgma_row$rank_logmeanexp,
    upgma_rank_max = upgma_row$rank_max,
    true_logmeanexp_minus_upgma =
      true_row$logmeanexp - upgma_row$logmeanexp,
    true_max_minus_upgma = true_row$max - upgma_row$max,
    top_20_logmeanexp_mean_dist =
      mean(head(scores_df[order(-scores_df$logmeanexp), "dist_to_true"], 20)),
    top_20_max_mean_dist =
      mean(head(scores_df[order(-scores_df$max), "dist_to_true"], 20))
  )
}

run_one_dataset <- function(experiment_name, sequence_length, rate, rate_index,
                            replicate, all_Fmats, support_keys, run_dir) {
  replicate_seed <- base_seed + replicate - 1
  support_seed <- base_seed + 100000 * rate_index + replicate
  label_seed <- base_seed + 200000 * rate_index + replicate

  set.seed(replicate_seed)
  init_results <- phylodyn:::generate_true_M_and_data(
    num_tips = num_tips,
    rate = rate,
    seq_len = sequence_length
  )

  true_tree <- init_results$M_true_tree
  sequences <- init_results$sequences

  inter_coal_times <- coalescent.intervals(true_tree)$interval.length
  inter_coal_times[inter_coal_times <= 0.001] <- 0.01
  coal_times <- cumsum(inter_coal_times)

  upgma_tree <- upgma(dist.hamming(sequences))
  upgma_timed_tree <- phylodyn:::update_time(upgma_tree, coal_times)

  F_true <- round(phylodyn:::gen_Fmat(true_tree, tol = 8), 0)
  F_upgma <- round(phylodyn:::gen_Fmat(upgma_timed_tree, tol = 8), 0)

  true_idx <- match(fmat_key(F_true), support_keys)
  upgma_idx <- match(fmat_key(F_upgma), support_keys)

  if (is.na(true_idx)) {
    stop("True F matrix was not found in F.list", num_tips,
         " for rate ", rate, ", replicate ", replicate)
  }

  set.seed(support_seed)
  indices <- sample_support_indices(
    all_indices = seq_along(all_Fmats),
    true_idx = true_idx,
    upgma_idx = upgma_idx,
    sample_size = support_sample_size
  )

  set.seed(label_seed)
  scores_df <- score_support_subset(
    indices = indices,
    all_Fmats = all_Fmats,
    F_true = F_true,
    F_upgma = F_upgma,
    coal_times = coal_times,
    sequences = sequences,
    rate = rate,
    R = random_label_R,
    label = paste0(experiment_name, " rate ", rate,
                   " replicate ", replicate)
  )

  scores_df$experiment <- experiment_name
  scores_df$seq_len <- sequence_length
  scores_df$rate <- rate
  scores_df$replicate <- replicate

  scores_path <- file.path(
    run_dir,
    sprintf(
      "scores_%s_rate_%s_replicate_%03d.csv",
      experiment_name,
      rate,
      replicate
    )
  )
  write.csv(scores_df, scores_path, row.names = FALSE)

  cbind(
    data.frame(
      experiment = experiment_name,
      seq_len = sequence_length,
      rate = rate,
      replicate = replicate,
      replicate_seed = replicate_seed,
      support_seed = support_seed,
      label_seed = label_seed,
      mean_pairwise_hamming = pairwise_hamming_mean(sequences),
      upgma_l2_to_true = distance_Fmat(F_upgma, F_true, dist = "l2")
    ),
    summarize_scores(scores_df)
  )
}

summarize_rate <- function(df) {
  data.frame(
    experiment = unique(df$experiment),
    rate = unique(df$rate),
    num_replicates = nrow(df),
    seq_len = unique(df$seq_len),
    support_sample_size = support_sample_size,
    random_label_R = random_label_R,
    mean_pairwise_hamming = mean(df$mean_pairwise_hamming),
    mean_cor_logmeanexp_neg_dist = mean(df$cor_logmeanexp_neg_dist),
    median_cor_logmeanexp_neg_dist = median(df$cor_logmeanexp_neg_dist),
    mean_spearman_logmeanexp_neg_dist =
      mean(df$spearman_logmeanexp_neg_dist),
    median_spearman_logmeanexp_neg_dist =
      median(df$spearman_logmeanexp_neg_dist),
    mean_cor_max_neg_dist = mean(df$cor_max_neg_dist),
    median_cor_max_neg_dist = median(df$cor_max_neg_dist),
    mean_spearman_max_neg_dist = mean(df$spearman_max_neg_dist),
    median_spearman_max_neg_dist = median(df$spearman_max_neg_dist),
    mean_true_logmeanexp_minus_upgma =
      mean(df$true_logmeanexp_minus_upgma),
    median_true_logmeanexp_minus_upgma =
      median(df$true_logmeanexp_minus_upgma),
    mean_true_max_minus_upgma = mean(df$true_max_minus_upgma),
    median_true_max_minus_upgma = median(df$true_max_minus_upgma)
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
support_keys <- vapply(all_Fmats, fmat_key, character(1))

cat("Run directory:", run_dir, "\n")
cat("num_tips:", num_tips, "num_replicates:", num_replicates, "\n")
cat("experiments:\n")
for (experiment in experiments) {
  cat("  ", experiment$name, "seq_len:", experiment$seq_len,
      "rate_grid:", paste(experiment$rate_grid, collapse = ", "), "\n")
}
cat("support_sample_size:", support_sample_size,
    "random_label_R:", random_label_R, "\n")

rows <- list()
row_id <- 1

for (experiment in experiments) {
  cat("\nExperiment:", experiment$name, "\n")

  for (rate_index in seq_along(experiment$rate_grid)) {
    rate <- experiment$rate_grid[rate_index]
    cat("\nRate", rate, "(", rate_index, "of",
        length(experiment$rate_grid), ")\n")

    for (replicate in seq_len(num_replicates)) {
      cat("  replicate", replicate, "of", num_replicates, "\n")
      rows[[row_id]] <- run_one_dataset(
        experiment_name = experiment$name,
        sequence_length = experiment$seq_len,
        rate = rate,
        rate_index = rate_index,
        replicate = replicate,
        all_Fmats = all_Fmats,
        support_keys = support_keys,
        run_dir = run_dir
      )
      row_id <- row_id + 1
    }
  }
}

replicate_summary <- do.call(rbind, rows)
rate_summary <- do.call(
  rbind,
  lapply(
    split(replicate_summary,
          list(replicate_summary$experiment, replicate_summary$rate),
          drop = TRUE),
    summarize_rate
  )
)

write.csv(
  replicate_summary,
  file.path(run_dir, "replicate_summary.csv"),
  row.names = FALSE
)
write.csv(
  rate_summary,
  file.path(run_dir, "rate_summary.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    settings = list(
      num_tips = num_tips,
      base_seed = base_seed,
      num_replicates = num_replicates,
      experiments = experiments,
      support_sample_size = support_sample_size,
      random_label_R = random_label_R
    ),
    replicate_summary = replicate_summary,
    rate_summary = rate_summary
  ),
  file.path(run_dir, "rate_support_correlation_results.rds")
)

cat("\nRate support-correlation summary:\n")
print(rate_summary, row.names = FALSE)
