## Oracle-style likelihood diagnostic for unlabeled F-matrix VB.
##
## Default run mirrors VB_master_exact_prior.R:
##   NUM_TIPS=10 RATE=0.01 SEQ_LEN=1000 SEED=42
##
## Optional environment overrides:
##   DATA_MODE=saved TREE_PATH=... FASTA_PATH=...
##   SCAN_R=10 REFINE_R=500 REFINE_TOP=100 OUTPUT_DIR=...

suppressPackageStartupMessages({
  library(ape)
  library(phylodyn)
  library(phyclust)
  library(phangorn)
  library(phylotools)
  library(fmatrix)
})

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  as.integer(value)
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  as.numeric(value)
}

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

num_tips <- env_int("NUM_TIPS", 10)
rate <- env_num("RATE", 0.01)
seq_len <- env_int("SEQ_LEN", 1000)
seed <- env_int("SEED", 42)
scan_R <- env_int("SCAN_R", 10)
refine_R <- env_int("REFINE_R", 500)
refine_top <- env_int("REFINE_TOP", 100)
data_mode <- Sys.getenv("DATA_MODE", unset = "simulate")

default_output_dir <- file.path(
  "Brandon", "gibbs_prior", "logs", "oracle_labeled_likelihood_diagnostic"
)
output_dir <- Sys.getenv("OUTPUT_DIR", unset = default_output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(output_dir, run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)

if (identical(data_mode, "saved")) {
  tree_path <- Sys.getenv("TREE_PATH", unset = file.path("Brandon", "gibbs_prior", "true_tree.newick"))
  fasta_path <- Sys.getenv("FASTA_PATH", unset = file.path("Brandon", "gibbs_prior", "sequences.fasta"))

  if (!file.exists(tree_path)) {
    stop("Could not find TREE_PATH: ", tree_path)
  }
  if (!file.exists(fasta_path)) {
    stop("Could not find FASTA_PATH: ", fasta_path)
  }

  true_tree <- read.tree(tree_path)
  sequences <- as.phyDat(read.FASTA(fasta_path))
} else if (identical(data_mode, "simulate")) {
  init_results <- phylodyn:::generate_true_M_and_data(
    num_tips = num_tips,
    rate = rate,
    seq_len = seq_len
  )
  true_tree <- init_results$M_true_tree
  sequences <- init_results$sequences
} else {
  stop("DATA_MODE must be either 'simulate' or 'saved'.")
}

inter_coal_times <- coalescent.intervals(true_tree)$interval.length
inter_coal_times[inter_coal_times <= 0.001] <- 0.01
coal_times <- cumsum(inter_coal_times)

upgma_tree <- upgma(dist.hamming(sequences))
upgma_timed_tree <- phylodyn:::update_time(upgma_tree, coal_times)

F_true <- round(phylodyn:::gen_Fmat(true_tree, tol = 8), 0)
F_init <- round(phylodyn:::gen_Fmat(upgma_timed_tree, tol = 8), 0)

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
true_idx <- match(fmat_key(F_true), support_keys)
init_idx <- match(fmat_key(F_init), support_keys)

if (is.na(true_idx)) {
  stop("True F matrix was not found in ", flist_name)
}
if (is.na(init_idx)) {
  warning("UPGMA init F matrix was not found in ", flist_name)
}

cat("Run directory:", run_dir, "\n")
cat("Data mode:", data_mode, "\n")
cat("num_tips:", num_tips, "rate:", rate, "seq_len:", seq_len, "seed:", seed, "\n")
cat("support size:", length(all_Fmats), "scan_R:", scan_R,
    "refine_R:", refine_R, "refine_top:", refine_top, "\n")
cat("true support idx:", true_idx, "init support idx:", init_idx, "\n")

true_labeled_ll <- pml_loglik(true_tree, sequences, rate)
upgma_labeled_ll <- pml_loglik(upgma_timed_tree, sequences, rate)

scan_indices <- seq_along(all_Fmats)
scan_df <- score_support(
  indices = scan_indices,
  all_Fmats = all_Fmats,
  F_true = F_true,
  F_init = F_init,
  coal_times = coal_times,
  sequences = sequences,
  rate = rate,
  R = scan_R,
  label = "scan"
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
  label = "refine"
)

refine_df$rank_logmeanexp_within_refined <- rank_desc(refine_df$logmeanexp)
refine_df$rank_max_within_refined <- rank_desc(refine_df$max)
refine_df$rank_meanlog_within_refined <- rank_desc(refine_df$meanlog)

true_random_scan <- scan_df[scan_df$is_true_F, , drop = FALSE]
init_random_scan <- scan_df[scan_df$is_init_F, , drop = FALSE]
true_random_refine <- refine_df[refine_df$is_true_F, , drop = FALSE]
init_random_refine <- refine_df[refine_df$is_init_F, , drop = FALSE]

true_cherries <- fmatrix:::num_cherries(F_true)
single_draw_true_label_prob <- 2^true_cherries / factorial(num_tips)

anchor_df <- data.frame(
  quantity = c(
    "true_labeled_tree",
    "upgma_labeled_tree",
    "true_F_random_labels_scan_logmeanexp",
    "true_F_random_labels_scan_max",
    "init_F_random_labels_scan_logmeanexp",
    "init_F_random_labels_scan_max",
    "true_F_random_labels_refine_logmeanexp",
    "true_F_random_labels_refine_max",
    "init_F_random_labels_refine_logmeanexp",
    "init_F_random_labels_refine_max"
  ),
  log_likelihood = c(
    true_labeled_ll,
    upgma_labeled_ll,
    true_random_scan$logmeanexp,
    true_random_scan$max,
    init_random_scan$logmeanexp,
    init_random_scan$max,
    true_random_refine$logmeanexp,
    true_random_refine$max,
    init_random_refine$logmeanexp,
    init_random_refine$max
  )
)

summary_df <- data.frame(
  data_mode = data_mode,
  num_tips = num_tips,
  rate = rate,
  seq_len = seq_len,
  seed = seed,
  support_size = length(all_Fmats),
  scan_R = scan_R,
  refine_R = refine_R,
  refine_top = refine_top,
  refined_size = nrow(refine_df),
  true_idx = true_idx,
  init_idx = init_idx,
  true_labeled_ll = true_labeled_ll,
  upgma_labeled_ll = upgma_labeled_ll,
  upgma_minus_true_labeled_ll = upgma_labeled_ll - true_labeled_ll,
  true_scan_rank_logmeanexp = true_random_scan$rank_logmeanexp,
  true_scan_rank_max = true_random_scan$rank_max,
  init_scan_rank_logmeanexp = init_random_scan$rank_logmeanexp,
  init_scan_rank_max = init_random_scan$rank_max,
  true_refine_rank_logmeanexp = true_random_refine$rank_logmeanexp_within_refined,
  true_refine_rank_max = true_random_refine$rank_max_within_refined,
  init_refine_rank_logmeanexp = init_random_refine$rank_logmeanexp_within_refined,
  init_refine_rank_max = init_random_refine$rank_max_within_refined,
  true_num_cherries = true_cherries,
  single_draw_true_label_prob = single_draw_true_label_prob,
  scan_hit_true_label_prob = 1 - (1 - single_draw_true_label_prob)^scan_R,
  refine_hit_true_label_prob = 1 - (1 - single_draw_true_label_prob)^refine_R,
  top_scan_logmeanexp_dist_mean = mean(head(scan_df[order(-scan_df$logmeanexp), "dist_to_true"], 20)),
  top_scan_max_dist_mean = mean(head(scan_df[order(-scan_df$max), "dist_to_true"], 20)),
  top_refine_logmeanexp_dist_mean = mean(head(refine_df[order(-refine_df$logmeanexp), "dist_to_true"], 20)),
  top_refine_max_dist_mean = mean(head(refine_df[order(-refine_df$max), "dist_to_true"], 20))
)

write.csv(scan_df, file.path(run_dir, "support_scan.csv"), row.names = FALSE)
write.csv(refine_df, file.path(run_dir, "support_refined.csv"), row.names = FALSE)
write.csv(anchor_df, file.path(run_dir, "anchor_likelihoods.csv"), row.names = FALSE)
write.csv(summary_df, file.path(run_dir, "summary.csv"), row.names = FALSE)
write.csv(head(scan_df[order(-scan_df$logmeanexp), ], 50),
          file.path(run_dir, "top_scan_logmeanexp.csv"), row.names = FALSE)
write.csv(head(scan_df[order(-scan_df$max), ], 50),
          file.path(run_dir, "top_scan_max.csv"), row.names = FALSE)
write.csv(head(refine_df[order(-refine_df$logmeanexp), ], 50),
          file.path(run_dir, "top_refined_logmeanexp.csv"), row.names = FALSE)
write.csv(head(refine_df[order(-refine_df$max), ], 50),
          file.path(run_dir, "top_refined_max.csv"), row.names = FALSE)

saveRDS(
  list(
    config = summary_df,
    anchors = anchor_df,
    scan = scan_df,
    refined = refine_df,
    F_true = F_true,
    F_init = F_init,
    true_tree = true_tree,
    upgma_timed_tree = upgma_timed_tree,
    coal_times = coal_times
  ),
  file.path(run_dir, "diagnostic_results.rds")
)

cat("\nAnchor likelihoods:\n")
print(anchor_df, row.names = FALSE)

cat("\nSummary:\n")
print(summary_df, row.names = FALSE)

cat("\nTop refined by logmeanexp:\n")
print(head(refine_df[order(-refine_df$logmeanexp), ], 10), row.names = FALSE)

cat("\nTop refined by max:\n")
print(head(refine_df[order(-refine_df$max), ], 10), row.names = FALSE)
