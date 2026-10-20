script_dir <- local({
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", command_args[grepl(file_arg, command_args)])

  if (length(script_path) > 0L && file.exists(script_path[[1L]])) {
    dirname(normalizePath(script_path[[1L]]))
  } else {
    candidate <- "Variational Ranked Methods/fmatrix_blocks/label_swap_likelihood_experiment.R"
    if (file.exists(candidate)) dirname(normalizePath(candidate)) else getwd()
  }
})

source(file.path(script_dir, "heterochronous_fmatrix_utils.R"))
source(file.path(script_dir, "swap_utils.R"))

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(phyclust)
  library(phylotools)
  library(ggplot2)
})

distance_columns <- c(D = "D_distance", "Phylodyn A" = "Amat_distance")

generate_tree_and_sequences <- function(n_tips,
                                        mutation_rate,
                                        sequence_length,
                                        tree_file,
                                        fasta_file) {
  simulated_tree <- phylodyn:::generate_newick(
    phylodyn:::coalsim(
      samp_times = 0,
      n_sampled = n_tips,
      traj = phylodyn:::exp_traj,
      method = "tt",
      val_upper = 11
    )
  )
  true_tree <- simulated_tree$newick
  true_tree$tip.label <- sub("_0", "", true_tree$tip.label)

  write.tree(true_tree, file = tree_file)

  temporary_phylip <- tempfile("seqgen_", fileext = ".phylip")
  on.exit(unlink(temporary_phylip), add = TRUE)

  seqgen_options <- paste0(
    "-mHKY -t2.0 -f0.25,0.25,0.25,0.25 -l",
    sequence_length,
    " -s",
    mutation_rate
  )
  seqgen(
    opts = seqgen_options,
    newick.tree = write.tree(true_tree),
    temp.file = temporary_phylip
  )

  dat2fasta(read.phylip(temporary_phylip), outfile = fasta_file)

  list(
    tree = true_tree,
    sequences = read.phyDat(fasta_file, format = "fasta", type = "DNA")
  )
}

tree_log_likelihood <- function(tree, sequences, mutation_rate) {
  fit <- pml(
    tree,
    sequences,
    bf = rep(0.25, 4),
    Q = c(1, 2, 1, 1, 2, 1),
    rate = mutation_rate,
    model = "USER"
  )
  as.numeric(logLik(fit))
}

scaled_D_and_Amat_distances <- function(true_tree, candidate_tree, leaf_order) {
  raw_distances <- tree_matrix_component_distances(
    true_tree,
    candidate_tree,
    matrices = c("D", "Amat"),
    norm = "l1",
    leaf_order = leaf_order
  )
  n_tips <- length(leaf_order)

  c(
    raw_D_distance = raw_distances[["D"]],
    raw_Amat_distance = raw_distances[["Amat"]],
    D_distance = raw_distances[["D"]] / n_tips^1.5,
    Amat_distance = raw_distances[["Amat"]] / n_tips^1.5
  )
}

summarize_metric_group <- function(group,
                                   metric,
                                   distance_column,
                                   group_name,
                                   settings) {
  distances <- group[[distance_column]]
  likelihood_changes <- group$delta_log_likelihood

  data.frame(
    n_tips = settings$n_tips,
    sequence_length = settings$sequence_length,
    mutation_rate = settings$mutation_rate,
    true_tree_log_likelihood = settings$true_tree_log_likelihood,
    metric = metric,
    n_swaps = group_name,
    n = nrow(group),
    mean_distance = mean(distances),
    median_distance = median(distances),
    sd_distance = sd(distances),
    mean_delta_log_likelihood = mean(likelihood_changes),
    median_delta_log_likelihood = median(likelihood_changes),
    sd_delta_log_likelihood = sd(likelihood_changes),
    pearson_correlation = suppressWarnings(cor(
      distances,
      likelihood_changes,
      method = "pearson"
    )),
    spearman_correlation = suppressWarnings(cor(
      distances,
      likelihood_changes,
      method = "spearman"
    )),
    row.names = NULL
  )
}

summarize_distance_likelihood_results <- function(candidate_results, settings) {
  groups <- c(
    list(all = candidate_results),
    split(candidate_results, candidate_results$n_swaps)
  )

  summary_rows <- lapply(names(groups), function(group_name) {
    group <- groups[[group_name]]
    do.call(rbind, lapply(names(distance_columns), function(metric) {
      summarize_metric_group(
        group,
        metric,
        distance_columns[[metric]],
        group_name,
        settings
      )
    }))
  })

  do.call(rbind, summary_rows)
}

plot_distance_vs_likelihood <- function(results, n_tips, output_dir) {
  long_results <- do.call(rbind, lapply(names(distance_columns), function(metric) {
    data.frame(
      metric = metric,
      distance = results[[distance_columns[[metric]]]],
      delta_log_likelihood = results$delta_log_likelihood,
      n_swaps = results$n_swaps,
      is_reference = results$is_reference
    )
  }))
  long_results$metric <- factor(
    long_results$metric,
    levels = c("D", "Phylodyn A")
  )

  candidate_results <- long_results[!long_results$is_reference, ]
  reference_results <- long_results[long_results$is_reference, ]
  candidate_results$n_swaps <- factor(candidate_results$n_swaps)

  scatter_plot <- ggplot(
    candidate_results,
    aes(x = distance, y = delta_log_likelihood, colour = n_swaps)
  ) +
    geom_point(alpha = 0.45) +
    geom_point(
      data = reference_results,
      aes(x = distance, y = delta_log_likelihood),
      inherit.aes = FALSE,
      colour = "black",
      shape = 4,
      size = 3
    ) +
    facet_wrap(~metric, nrow = 1, scales = "free_x") +
    scale_colour_brewer(palette = "Set2") +
    labs(
      title = paste0("Distance from the true tree versus likelihood (", n_tips, " tips)"),
      x = "Scaled L1 distance",
      y = "Log likelihood relative to the true tree",
      colour = "Number of swaps"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(
    file.path(output_dir, "distance_vs_likelihood.png"),
    scatter_plot,
    width = 12,
    height = 6,
    dpi = 300
  )
  ggsave(
    file.path(output_dir, "distance_vs_likelihood.pdf"),
    scatter_plot,
    width = 12,
    height = 6
  )

  scatter_plot
}

run_label_swap_likelihood_experiment <- function(
    n_tips = 20L,
    n_iterations = 1000L,
    n_swaps_values = c(1L, 3L, 5L, 10L, 20L),
    mutation_rate = 0.01,
    sequence_length = 1000L,
    seed = NULL,
    results_dir = file.path(script_dir, "label_swap_likelihood_results"),
    show_progress = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  n_swaps_values <- sort(unique(n_swaps_values))

  tip_results_dir <- file.path(results_dir, paste0("n_tips_", n_tips))
  if (dir.exists(tip_results_dir)) unlink(tip_results_dir, recursive = TRUE)
  dir.create(tip_results_dir, recursive = TRUE)

  generated_data <- generate_tree_and_sequences(
    n_tips = n_tips,
    mutation_rate = mutation_rate,
    sequence_length = sequence_length,
    tree_file = file.path(tip_results_dir, "true_tree.newick"),
    fasta_file = file.path(tip_results_dir, "sequences.fasta")
  )
  true_tree <- generated_data$tree
  sequences <- generated_data$sequences
  leaf_order <- paste0("t", seq_len(n_tips))
  true_tree_log_likelihood <- tree_log_likelihood(true_tree, sequences, mutation_rate)
  experiment_settings <- list(
    n_tips = n_tips,
    n_iterations = n_iterations,
    n_swaps_values = n_swaps_values,
    mutation_rate = mutation_rate,
    sequence_length = sequence_length,
    seed = seed,
    true_tree_log_likelihood = true_tree_log_likelihood,
    leaf_order = leaf_order
  )

  total_candidates <- n_iterations * length(n_swaps_values)
  progress_bar <- if (show_progress) {
    txtProgressBar(min = 0, max = total_candidates, style = 3)
  }
  if (!is.null(progress_bar)) on.exit(close(progress_bar), add = TRUE)

  result_rows <- vector("list", total_candidates)
  row_index <- 1L
  for (iteration in seq_len(n_iterations)) {
    for (n_swaps in n_swaps_values) {
      candidate_tree <- apply_swap_n_times(true_tree, random_tip_label_swap, n_swaps)
      distances <- scaled_D_and_Amat_distances(true_tree, candidate_tree, leaf_order)
      candidate_log_likelihood <- tree_log_likelihood(candidate_tree, sequences, mutation_rate)

      result_rows[[row_index]] <- data.frame(
        iteration = iteration,
        n_swaps = n_swaps,
        as.list(distances),
        log_likelihood = candidate_log_likelihood,
        delta_log_likelihood = candidate_log_likelihood - true_tree_log_likelihood,
        is_reference = FALSE,
        row.names = NULL
      )
      row_index <- row_index + 1L
      if (!is.null(progress_bar)) setTxtProgressBar(progress_bar, row_index - 1L)
    }
  }

  reference_row <- data.frame(
    iteration = 0L,
    n_swaps = 0L,
    raw_D_distance = 0,
    raw_Amat_distance = 0,
    D_distance = 0,
    Amat_distance = 0,
    log_likelihood = true_tree_log_likelihood,
    delta_log_likelihood = 0,
    is_reference = TRUE
  )
  results <- rbind(reference_row, do.call(rbind, result_rows))
  candidate_results <- results[!results$is_reference, ]
  summary <- summarize_distance_likelihood_results(candidate_results, experiment_settings)

  write.csv(
    results,
    file.path(tip_results_dir, "distance_likelihood_results.csv"),
    row.names = FALSE
  )
  write.csv(
    summary,
    file.path(tip_results_dir, "summary.csv"),
    row.names = FALSE
  )
  plot <- plot_distance_vs_likelihood(results, n_tips, tip_results_dir)
  saveRDS(experiment_settings, file.path(tip_results_dir, "config.rds"))

  list(
    results = results,
    summary = summary,
    true_tree = true_tree,
    sequences = sequences,
    plot = plot,
    results_dir = tip_results_dir
  )
}

initial_label_swap_likelihood_results <- run_label_swap_likelihood_experiment(n_tips = 80,
  n_swaps_values = c(1L, 3L, 5L, 10L, 20L),
  seed = 42L
)
