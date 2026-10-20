script_dir <- local({
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", command_args[grepl(file_arg, command_args)])

  if (length(script_path) > 0L && file.exists(script_path[[1L]])) {
    dirname(normalizePath(script_path[[1L]]))
  } else {
    candidate <- "Variational Ranked Methods/fmatrix_blocks/swap_distance_experiment.R"
    if (file.exists(candidate)) dirname(normalizePath(candidate)) else getwd()
  }
})

source(file.path(script_dir, "swap_distance_experiment.R"))

default_component_scaling <- c(D = 1.5, Amat = 1.5, F_A = 2)

add_joint_metrics <- function(component_distances, component_weight = 4) {
  D <- component_distances[["D"]]
  Amat <- component_distances[["Amat"]]
  F_A <- component_distances[["F_A"]]

  data.frame(
    scaled_D = D,
    scaled_Amat = Amat,
    D_plus_F_A = D + F_A,
    Amat_plus_F_A = Amat + F_A,
    D_heavy_plus_F_A = component_weight * D + F_A,
    F_A_heavy_plus_D = D + component_weight * F_A,
    Amat_heavy_plus_F_A = component_weight * Amat + F_A,
    F_A_heavy_plus_Amat = Amat + component_weight * F_A,
    check.names = FALSE
  )
}

scale_component_distances <- function(component_distances,
                                      n_tips,
                                      component_scaling = default_component_scaling) {
  data.frame(
    D = component_distances[["D"]] / n_tips^component_scaling[["D"]],
    Amat = component_distances[["Amat"]] / n_tips^component_scaling[["Amat"]],
    F_A = component_distances[["F_A"]] / n_tips^component_scaling[["F_A"]]
  )
}

metric_names <- function() {
  c(
    "scaled_D",
    "scaled_Amat",
    "D_plus_F_A",
    "Amat_plus_F_A",
    "D_heavy_plus_F_A",
    "F_A_heavy_plus_D",
    "Amat_heavy_plus_F_A",
    "F_A_heavy_plus_Amat"
  )
}

metric_labels <- function(component_weight, component_scaling) {
  D <- paste0("D / N^", component_scaling[["D"]])
  Amat <- paste0("Phylodyn A / N^", component_scaling[["Amat"]])
  F_A <- paste0("F_A / N^", component_scaling[["F_A"]])

  c(
    scaled_D = D,
    scaled_Amat = Amat,
    D_plus_F_A = paste(D, "+", F_A),
    Amat_plus_F_A = paste(Amat, "+", F_A),
    D_heavy_plus_F_A = paste0("D-heavy: ", component_weight, "(", D, ") + ", F_A),
    F_A_heavy_plus_D = paste0("F_A-heavy: ", D, " + ", component_weight, "(", F_A, ")"),
    Amat_heavy_plus_F_A = paste0("A-heavy: ", component_weight, "(", Amat, ") + ", F_A),
    F_A_heavy_plus_Amat = paste0("F_A-heavy: ", Amat, " + ", component_weight, "(", F_A, ")")
  )
}

as_long_metric_distances <- function(distances) {
  metrics <- metric_names()
  rows <- lapply(metrics, function(metric) {
    data.frame(
      n_tips = distances$n_tips,
      iteration = distances$iteration,
      n_swaps = distances$n_swaps,
      swap_type = distances$swap_type,
      metric = metric,
      distance = distances[[metric]],
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

summarize_metric_distances <- function(long_distances) {
  groups <- split(
    long_distances,
    list(long_distances$swap_type, long_distances$n_swaps, long_distances$metric),
    drop = TRUE
  )

  summary_rows <- lapply(groups, function(group) {
    values <- group$distance
    data.frame(
      n_tips = group$n_tips[[1L]],
      n_swaps = group$n_swaps[[1L]],
      swap_type = group$swap_type[[1L]],
      metric = group$metric[[1L]],
      n = length(values),
      mean = mean(values),
      median = median(values),
      sd = sd(values),
      iqr = IQR(values),
      min = min(values),
      p05 = unname(quantile(values, 0.05)),
      p25 = unname(quantile(values, 0.25)),
      p75 = unname(quantile(values, 0.75)),
      p95 = unname(quantile(values, 0.95)),
      max = max(values),
      prop_zero = mean(values == 0),
      row.names = NULL
    )
  })

  do.call(rbind, summary_rows)
}

write_combined_summary <- function(results_dir, summary_data = NULL) {
  if (is.null(summary_data)) {
    result_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
    summary_files <- file.path(result_dirs, "summary.csv")
    summary_files <- summary_files[file.exists(summary_files)]
    if (length(summary_files) == 0L) return(invisible(NULL))
    summaries <- lapply(summary_files, utils::read.csv)
    all_columns <- unique(unlist(lapply(summaries, names)))
    summaries <- lapply(summaries, function(summary) {
      missing_columns <- setdiff(all_columns, names(summary))
      summary[missing_columns] <- NA
      summary[all_columns]
    })
    summary_data <- do.call(rbind, summaries)
  }

  utils::write.csv(
    summary_data,
    file.path(results_dir, "summary_all_n_tips.csv"),
    row.names = FALSE
  )
  summary_data
}

plot_metric_histograms <- function(long_distances,
                                   n_tips,
                                   output_dir,
                                   component_weight,
                                   component_scaling) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to make histogram plots.")
  }

  labels <- metric_labels(component_weight, component_scaling)
  long_distances$metric <- factor(
    long_distances$metric,
    levels = names(labels),
    labels = labels
  )

  move_labels <- c(label_swap = "Tip-label swap", nni_swap = "NNI swap")
  swap_names <- unique(as.character(long_distances$swap_type))

  lapply(swap_names, function(swap_name) {
    move_distances <- long_distances[long_distances$swap_type == swap_name, ]
    move_distances$n_swaps <- factor(
      move_distances$n_swaps,
      levels = sort(unique(move_distances$n_swaps))
    )
    histogram <- ggplot2::ggplot(
      move_distances,
      ggplot2::aes(x = distance, fill = n_swaps)
    ) +
      ggplot2::geom_histogram(bins = 30, position = "identity", alpha = 0.45) +
      ggplot2::facet_wrap(~metric, scales = "free", ncol = 2) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::labs(
        title = paste0(
          "Matrix-distance distributions (", n_tips, " tips): ",
          if (swap_name %in% names(move_labels)) {
            move_labels[[swap_name]]
          } else {
            gsub("_", " ", swap_name)
          }
        ),
        x = "Scaled L1 distance",
        y = "Count",
        fill = "Number of swaps"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(legend.position = "bottom")

    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(swap_name, "_distance_histograms.png")),
      plot = histogram,
      width = 12,
      height = 14,
      dpi = 300
    )
    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(swap_name, "_distance_histograms.pdf")),
      plot = histogram,
      width = 12,
      height = 14
    )

    histogram
  })
}

run_distance_comparison_experiment <- function(
    n_tips,
    n_iterations = 1000L,
    n_swaps_values = c(1L, 3L, 5L, 10L, 20L),
    tree_mode = c("new_tree_each_iteration", "same_tree"),
    swap_fns = list(label_swap = random_tip_label_swap, nni_swap = random_nni_swap),
    leaf_order = paste0("t", seq_len(n_tips)),
    component_weight = 4,
    component_scaling = default_component_scaling,
    seed = NULL,
    results_dir = file.path(script_dir, "swap_distance_results"),
    progress_fn = NULL,
    update_combined_summary = TRUE) {
  tree_mode <- match.arg(tree_mode)

  if (!is.null(seed)) set.seed(seed)

  tip_results_dir <- file.path(results_dir, paste0("n_tips_", n_tips))
  if (dir.exists(tip_results_dir)) unlink(tip_results_dir, recursive = TRUE)
  dir.create(file.path(tip_results_dir, "plots"), recursive = TRUE)

  base_trees <- generate_base_trees(
    n_tips = n_tips,
    n_iterations = n_iterations,
    tree_mode = tree_mode,
    leaf_order = leaf_order
  )
  component_distance_fn <- function(tree1, tree2) {
    tree_matrix_component_distances(
      tree1,
      tree2,
      matrices = c("D", "Amat", "F_A"),
      leaf_order = leaf_order
    )
  }
  simulations <- lapply(swap_fns, function(swap_fn) {
    simulate_swap_distances(
      n_tips = n_tips,
      n_iterations = n_iterations,
      n_swaps = n_swaps_values,
      tree_mode = tree_mode,
      swap_fn = swap_fn,
      distance_fns = list(component_distance_fn),
      base_trees = base_trees,
      progress_fn = progress_fn
    )
  })

  distance_rows <- lapply(names(simulations), function(swap_name) {
    components <- simulations[[swap_name]]$distances
    scaled_components <- scale_component_distances(
      components[, c("D", "Amat", "F_A")],
      n_tips = n_tips,
      component_scaling = component_scaling
    )
    metrics <- add_joint_metrics(
      scaled_components,
      component_weight = component_weight
    )

    data.frame(
      n_tips = n_tips,
      iteration = components$iteration,
      n_swaps = components$n_swaps,
      swap_type = swap_name,
      raw_D = components$D,
      raw_Amat = components$Amat,
      raw_F_A = components$F_A,
      metrics,
      scaled_F_A = scaled_components$F_A,
      check.names = FALSE
    )
  })
  distances <- do.call(rbind, distance_rows)
  long_distances <- as_long_metric_distances(distances)
  summary <- summarize_metric_distances(long_distances)
  saveRDS(
    list(
      n_tips = n_tips,
      n_iterations = n_iterations,
      n_swaps_values = n_swaps_values,
      tree_mode = tree_mode,
      swap_types = names(swap_fns),
      component_weight = component_weight,
      component_scaling = component_scaling,
      seed = seed
    ),
    file.path(tip_results_dir, "config.rds")
  )
  utils::write.csv(distances, file.path(tip_results_dir, "distances.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(tip_results_dir, "summary.csv"), row.names = FALSE)
  plot_metric_histograms(
    long_distances = long_distances,
    n_tips = n_tips,
    output_dir = file.path(tip_results_dir, "plots"),
    component_weight = component_weight,
    component_scaling = component_scaling
  )
  if (update_combined_summary) write_combined_summary(results_dir)

  list(
    distances = distances,
    summary = summary,
    results_dir = tip_results_dir
  )
}

run_all_distance_comparison_experiments <- function(
    n_tips_values = c(5L, 10L, 20L, 40L, 80L),
    n_iterations = 1000L,
    n_swaps_values = c(1L, 3L, 5L, 10L, 20L),
    tree_mode = "new_tree_each_iteration",
    component_weight = 4,
    component_scaling = default_component_scaling,
    seed = NULL,
    results_dir = file.path(script_dir, "swap_distance_results"),
    swap_fns = list(label_swap = random_tip_label_swap, nni_swap = random_nni_swap)) {
  total_steps <- length(n_tips_values) * n_iterations * length(n_swaps_values) *
    length(swap_fns)
  progress_bar <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)
  completed_steps <- 0L
  advance_progress <- function() {
    completed_steps <<- completed_steps + 1L
    utils::setTxtProgressBar(progress_bar, completed_steps)
  }
  on.exit(close(progress_bar), add = TRUE)

  results <- lapply(seq_along(n_tips_values), function(index) {
    run_distance_comparison_experiment(
      n_tips = n_tips_values[[index]],
      n_iterations = n_iterations,
      n_swaps_values = n_swaps_values,
      tree_mode = tree_mode,
      component_weight = component_weight,
      component_scaling = component_scaling,
      seed = if (is.null(seed)) NULL else seed + index - 1L,
      results_dir = results_dir,
      swap_fns = swap_fns,
      progress_fn = advance_progress,
      update_combined_summary = FALSE
    )
  })

  combined_summary <- do.call(rbind, lapply(results, `[[`, "summary"))
  write_combined_summary(results_dir, summary_data = combined_summary)

  list(results = results, summary = combined_summary)
}

initial_distance_comparison_results <- run_all_distance_comparison_experiments(
  n_tips_values = c(5L, 10L, 20L, 40L, 80L),
  n_iterations = 1000L,
  n_swaps_values = c(1L, 3L, 5L, 10L, 20L),
  seed = 20261018L
)
