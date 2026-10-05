## Focused 100-replicate check for seq_len = 10000, rate = 0.01.
##
## Reports UPGMA - true likelihood differences. Negative values mean the
## true tree/F scored better.

# ---- Experiment settings ----
num_tips <- 10
seq_len <- 10000
rate <- 0.01
base_seed <- 42
num_replicates <- 100
random_label_R <- 500

output_dir <- file.path(
  "Brandon", "gibbs_prior", "logs", "rate001_seq10000_100rep_diagnostic"
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

pairwise_distance_summary <- function(sequences, true_tree, rate) {
  hamming_mat <- as.matrix(dist.hamming(sequences))
  hamming_vals <- hamming_mat[upper.tri(hamming_mat)]

  tree_dist <- cophenetic(true_tree) * rate
  tree_dist_vals <- tree_dist[upper.tri(tree_dist)]

  data.frame(
    mean_pairwise_hamming = mean(hamming_vals),
    median_pairwise_hamming = median(hamming_vals),
    max_pairwise_hamming = max(hamming_vals),
    mean_pairwise_subs = mean(tree_dist_vals),
    median_pairwise_subs = median(tree_dist_vals),
    max_pairwise_subs = max(tree_dist_vals),
    tree_height_subs = max(node.depth.edgelength(true_tree)) * rate
  )
}

run_one_dataset <- function(replicate) {
  replicate_seed <- base_seed + replicate - 1
  random_label_seed <- base_seed + 100000 + replicate

  set.seed(replicate_seed)
  init_results <- phylodyn:::generate_true_M_and_data(
    num_tips = num_tips,
    rate = rate,
    seq_len = seq_len
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

  true_labeled_ll <- pml_loglik(true_tree, sequences, rate)
  upgma_labeled_ll <- pml_loglik(upgma_timed_tree, sequences, rate)

  set.seed(random_label_seed)
  true_random <- score_fmat_random_labels(
    Fmat = F_true,
    coal_times = coal_times,
    sequences = sequences,
    rate = rate,
    R = random_label_R
  )

  set.seed(random_label_seed + 1)
  upgma_random <- score_fmat_random_labels(
    Fmat = F_upgma,
    coal_times = coal_times,
    sequences = sequences,
    rate = rate,
    R = random_label_R
  )

  distance_summary <- pairwise_distance_summary(sequences, true_tree, rate)

  data.frame(
    replicate = replicate,
    replicate_seed = replicate_seed,
    seq_len = seq_len,
    rate = rate,
    random_label_R = random_label_R,
    init_l2_to_true = distance_Fmat(F_upgma, F_true, dist = "l2"),
    true_labeled_ll = true_labeled_ll,
    upgma_labeled_ll = upgma_labeled_ll,
    labeled_upgma_minus_true = upgma_labeled_ll - true_labeled_ll,
    labeled_true_beats_upgma = true_labeled_ll > upgma_labeled_ll,
    true_random_logmeanexp = true_random$logmeanexp,
    upgma_random_logmeanexp = upgma_random$logmeanexp,
    unlabeled_logmeanexp_upgma_minus_true =
      upgma_random$logmeanexp - true_random$logmeanexp,
    unlabeled_logmeanexp_true_beats_upgma =
      true_random$logmeanexp > upgma_random$logmeanexp,
    true_random_max = true_random$max,
    upgma_random_max = upgma_random$max,
    unlabeled_max_upgma_minus_true = upgma_random$max - true_random$max,
    unlabeled_max_true_beats_upgma = true_random$max > upgma_random$max,
    true_random_meanlog = true_random$meanlog,
    upgma_random_meanlog = upgma_random$meanlog,
    true_random_sd = true_random$sd,
    upgma_random_sd = upgma_random$sd,
    true_labeled_minus_true_random_logmeanexp =
      true_labeled_ll - true_random$logmeanexp,
    upgma_labeled_minus_upgma_random_logmeanexp =
      upgma_labeled_ll - upgma_random$logmeanexp,
    distance_summary
  )
}

summarize_vector <- function(x, prefix) {
  out <- data.frame(
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    se = sd(x) / sqrt(length(x)),
    q05 = unname(quantile(x, 0.05)),
    q25 = unname(quantile(x, 0.25)),
    q75 = unname(quantile(x, 0.75)),
    q95 = unname(quantile(x, 0.95))
  )
  names(out) <- paste(prefix, names(out), sep = "_")
  out
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(output_dir, run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

cat("Run directory:", run_dir, "\n")
cat("num_tips:", num_tips, "seq_len:", seq_len, "rate:", rate,
    "num_replicates:", num_replicates, "\n")
cat("random_label_R:", random_label_R, "\n")

rows <- vector("list", num_replicates)

for (replicate in seq_len(num_replicates)) {
  cat("replicate", replicate, "of", num_replicates, "\n")
  rows[[replicate]] <- run_one_dataset(replicate)
}

replicate_results <- do.call(rbind, rows)

summary_df <- cbind(
  data.frame(
    num_replicates = num_replicates,
    seq_len = seq_len,
    rate = rate,
    random_label_R = random_label_R,
    mean_pairwise_hamming = mean(replicate_results$mean_pairwise_hamming),
    median_pairwise_hamming = median(replicate_results$mean_pairwise_hamming),
    mean_pairwise_subs = mean(replicate_results$mean_pairwise_subs),
    mean_init_l2_to_true = mean(replicate_results$init_l2_to_true),
    labeled_true_beats_upgma_count =
      sum(replicate_results$labeled_true_beats_upgma),
    labeled_true_beats_upgma_fraction =
      mean(replicate_results$labeled_true_beats_upgma),
    unlabeled_logmeanexp_true_beats_upgma_count =
      sum(replicate_results$unlabeled_logmeanexp_true_beats_upgma),
    unlabeled_logmeanexp_true_beats_upgma_fraction =
      mean(replicate_results$unlabeled_logmeanexp_true_beats_upgma),
    unlabeled_max_true_beats_upgma_count =
      sum(replicate_results$unlabeled_max_true_beats_upgma),
    unlabeled_max_true_beats_upgma_fraction =
      mean(replicate_results$unlabeled_max_true_beats_upgma)
  ),
  summarize_vector(
    replicate_results$labeled_upgma_minus_true,
    "labeled_upgma_minus_true"
  ),
  summarize_vector(
    replicate_results$unlabeled_logmeanexp_upgma_minus_true,
    "unlabeled_logmeanexp_upgma_minus_true"
  ),
  summarize_vector(
    replicate_results$unlabeled_max_upgma_minus_true,
    "unlabeled_max_upgma_minus_true"
  ),
  summarize_vector(
    replicate_results$true_labeled_minus_true_random_logmeanexp,
    "true_labeling_gap"
  )
)

write.csv(
  replicate_results,
  file.path(run_dir, "replicate_results.csv"),
  row.names = FALSE
)
write.csv(
  summary_df,
  file.path(run_dir, "summary.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    settings = list(
      num_tips = num_tips,
      seq_len = seq_len,
      rate = rate,
      base_seed = base_seed,
      num_replicates = num_replicates,
      random_label_R = random_label_R
    ),
    replicate_results = replicate_results,
    summary = summary_df
  ),
  file.path(run_dir, "rate001_seq10000_100rep_results.rds")
)

cat("\nSummary:\n")
print(summary_df, row.names = FALSE)
