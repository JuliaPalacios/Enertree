suppressPackageStartupMessages({
  library(stats)
})

script_dir <- local({
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", command_args[grepl(file_arg, command_args)])

  if (length(script_path) > 0L && file.exists(script_path[[1L]])) {
    dirname(normalizePath(script_path[[1L]]))
  } else {
    candidate <- "Variational Ranked Methods/fmatrix_blocks/heterochronous_fmatrix_utils.R"
    if (file.exists(candidate)) {
      dirname(normalizePath(candidate))
    } else {
      getwd()
    }
  }
})

source(file.path(script_dir, "heterochronous_fmatrix_utils.R"))

simulate_sampler <- function(diagonal,
                             sampler,
                             n_replicates,
                             seed,
                             all_keys,
                             show_progress = TRUE,
                             progress_label = "sampling") {
  set.seed(seed)
  counts <- setNames(numeric(length(all_keys)), all_keys)

  progress_bar <- NULL
  if (show_progress) {
    progress_bar <- txtProgressBar(
      min = 0,
      max = n_replicates,
      initial = 0,
      style = 3
    )
    cat("\n", progress_label, "\n", sep = "")
  }

  progress_every <- max(1L, floor(n_replicates / 100L))
  start_time <- proc.time()[["elapsed"]]

  for (i in seq_len(n_replicates)) {
    key <- matrix_key(sampler(diagonal))
    counts[[key]] <- counts[[key]] + 1

    if (show_progress && (i %% progress_every == 0L || i == n_replicates)) {
      setTxtProgressBar(progress_bar, i)
    }
  }

  if (show_progress) {
    close(progress_bar)
  }

  elapsed <- proc.time()[["elapsed"]] - start_time
  cat(progress_label, "elapsed seconds:", round(elapsed, 2), "\n")

  counts[counts > 0]
}

summarize_sample_counts <- function(observed_counts, all_keys) {
  counts <- setNames(integer(length(all_keys)), all_keys)
  counts[names(observed_counts)] <- as.integer(observed_counts)

  expected <- rep(sum(counts) / length(all_keys), length(all_keys))
  chisq <- suppressWarnings(
    chisq.test(x = counts, p = rep(1 / length(all_keys), length(all_keys)))
  )

  data.frame(
    n_states = length(all_keys),
    n_replicates = sum(counts),
    expected_per_state = expected[[1L]],
    min_count = min(counts),
    max_count = max(counts),
    mean_count = mean(counts),
    sd_count = sd(counts),
    chisq_statistic = unname(chisq$statistic),
    chisq_df = unname(chisq$parameter),
    chisq_p_value = chisq$p.value
  )
}

n_tips <- 5
n_replicates <- 10000000
show_progress <- TRUE
diagonal <- fixed_isochronous_completion_diagonal(n_tips)

enumerated <- enumerate_fixed_diagonal_fmats(diagonal)
all_keys <- vapply(enumerated, function(state) matrix_key(state$F), character(1))
exact_probabilities <- vapply(enumerated, function(state) state$probability, numeric(1))

cat("Fixed diagonal:\n")
print(diagonal)

cat("\nNumber of enumerated matrices:\n")
print(length(enumerated))

cat("\nClosed-form count for this diagonal, prod(diagonal) / 2^(n - 1):\n")
print(prod(diagonal) / 2^(n_tips - 1L))

cat("\nExact probabilities from fixed-diagonal sampler over enumerated matrices:\n")
print(summary(exact_probabilities))
cat("Number of unique exact probabilities:\n")
print(length(unique(round(exact_probabilities, 14))))
cat("Sum of exact probabilities:\n")
print(sum(exact_probabilities))

sample_counts <- simulate_sampler(
  diagonal = diagonal,
  sampler = sample_fixed_diagonal_f_matrix,
  n_replicates = n_replicates,
  seed = 1,
  all_keys = all_keys,
  show_progress = show_progress,
  progress_label = "Fixed-diagonal sampler"
)

cat("\nSimulation summary for fixed-diagonal sampler:\n")
print(summarize_sample_counts(sample_counts, all_keys))
