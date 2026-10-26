# Compare final A-matrix tree draws with BEAST, following the rooted RF and
# cherry-count comparisons in Variational Ranked Methods/rf_boltzmann.
#
# Rscript compare_tree_samples.R [EXPERIMENT_ID_OR_OUTPUT_DIR] [BEAST_TREES] [TRUE_TREE] [PLOTS_DIR]
# In RStudio, set experiment_id below to the output folder name you want.
# With no selection, a single exported output folder is used automatically.
# By default, BEAST trees must be inside the selected experiment's output
# directory. Optional paths can override the BEAST file and true tree.

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggplot2)
  library(TreeTools)
  library(coda)
})

# Edit this when running the whole file with RStudio's Run button.
# Example: experiment_id <- "20261024_184601"
# A command-line experiment argument takes precedence over this setting.
experiment_id <- "20261025_231748"

source_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
compare_console_run <- FALSE
if (is.null(source_file)) {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg)) source_file <- sub("^--file=", "", file_arg[1L])
}
if (!is.null(source_file) && !file.exists(source_file)) {
  source_file <- gsub("~\\+~", " ", source_file)
}
if (is.null(source_file) || !file.exists(source_file)) {
  # RStudio's Run button sends lines to the Console without an `ofile` value.
  candidates <- file.path(getwd(), c(
    "A matrix analysis/compare_tree_samples.R", "compare_tree_samples.R",
    "Brandon/gibbs_prior/A matrix analysis/compare_tree_samples.R"
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    stop("Cannot find compare_tree_samples.R from working directory: ", getwd(),
         ". Set the working directory to gibbs_prior or run with Rscript.")
  }
  source_file <- candidates[1L]
  compare_console_run <- TRUE
}
analysis_dir <- dirname(normalizePath(source_file, mustWork = TRUE))
repo_dir <- normalizePath(file.path(analysis_dir, "..", "..", ".."), mustWork = TRUE)

default_experiment <- function() {
  root <- file.path(analysis_dir, "outputs")
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "sequences.fasta")) &
                 (file.exists(file.path(dirs, "sampled_trees.newick")) |
                    file.exists(file.path(dirs, "sampled_trees.trees")))]
  if (!length(dirs)) stop("No exported experiment found in: ", root)
  if (length(dirs) > 1L) {
    stop("Multiple exported experiments found. Set experiment_id near the top of ",
         "compare_tree_samples.R or pass an experiment ID on the command line. ",
         "Available IDs: ", paste(sort(basename(dirs)), collapse = ", "))
  }
  dirs[1L]
}

read_tree_sample <- function(path) {
  if (!file.exists(path)) stop("Tree file does not exist: ", path)
  header <- readLines(path, n = 20L, warn = FALSE)
  trees <- if (any(grepl("#NEXUS", header, ignore.case = TRUE))) {
    ape::read.nexus(path)
  } else {
    ape::read.tree(path)
  }
  if (inherits(trees, "phylo")) trees <- list(trees)
  if (!length(trees)) stop("No trees found in: ", path)
  class(trees) <- "multiPhylo"
  trees
}

same_alignment <- function(first, second) {
  x <- tryCatch(suppressWarnings(ape::read.FASTA(first)), error = function(e) NULL)
  y <- tryCatch(suppressWarnings(ape::read.FASTA(second)), error = function(e) NULL)
  if (is.null(x) || is.null(y) || !setequal(names(x), names(y)) ||
      anyDuplicated(names(x)) || anyDuplicated(names(y))) return(FALSE)
  labels <- sort(names(x))
  identical(unname(as.matrix(x[labels])), unname(as.matrix(y[labels])))
}

find_experiment_beast <- function(output_dir) {
  files <- list.files(output_dir, pattern = "\\.trees$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  files <- files[basename(files) != "sampled_trees.trees"]
  named_beast <- files[grepl("beast", basename(files), ignore.case = TRUE)]
  if (length(named_beast)) files <- named_beast
  if (length(files) != 1L) {
    stop("Expected one BEAST .trees file inside ", output_dir, ", found ",
         length(files), ". Add it there or pass its path as argument 2.",
         if (length(files)) paste0(" Matches: ", paste(files, collapse = "; ")))
  }
  normalizePath(files, mustWork = TRUE)
}

find_true_tree <- function(output_dir, export_fasta) {
  files <- list.files(output_dir, pattern = "\\.newick$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  files <- files[grepl("true", basename(files), ignore.case = TRUE)]
  if (!length(files)) {
    # Simulated fits keep the true tree in result.rds. Recover it only when
    # that fit's alignment exactly matches the exported FASTA.
    summary_path <- file.path(output_dir, "export_summary.rds")
    source_path <- if (file.exists(summary_path)) {
      readRDS(summary_path)$source_result
    } else NULL
    result_paths <- unique(c(source_path, file.path(
      dirname(analysis_dir), "logs", "VB_master_A",
      basename(output_dir), "result.rds"
    )))
    for (result_path in result_paths[file.exists(result_paths)]) {
      fit <- readRDS(result_path)
      if (!inherits(fit$true_tree, "phylo") ||
          !inherits(fit$sequences, "phyDat")) next
      exported <- ape::read.FASTA(export_fasta)
      saved <- ape::as.DNAbin(fit$sequences)
      labels <- sort(names(exported))
      if (!setequal(labels, rownames(saved)) ||
          !setequal(labels, fit$true_tree$tip.label) ||
          !identical(unname(as.matrix(exported[labels])),
                     unname(as.matrix(saved[labels, , drop = FALSE])))) next
      true_path <- file.path(output_dir, "true_tree.newick")
      ape::write.tree(fit$true_tree, file = true_path)
      return(normalizePath(true_path, mustWork = TRUE))
    }
  }
  if (!length(files)) {
    # A simulation's source directory may contain the true tree even when the
    # comparison samples are kept together in a separate experiment folder.
    fasta <- list.files(repo_dir, pattern = "\\.(fasta|fa|fas)$",
                        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    matching_fasta <- fasta[vapply(fasta, same_alignment, logical(1L),
                                   second = export_fasta)]
    files <- unlist(lapply(dirname(matching_fasta), function(directory) {
      candidates <- list.files(directory, pattern = "\\.newick$",
                               full.names = TRUE, ignore.case = TRUE)
      candidates[grepl("true", basename(candidates), ignore.case = TRUE)]
    }), use.names = FALSE)
    files <- unique(files)
  }
  if (length(files) != 1L) {
    stop("Expected one true-tree .newick for this experiment, found ",
         length(files), ". Pass the true-tree file as argument 3.")
  }
  normalizePath(files, mustWork = TRUE)
}

summarize_sample <- function(trees, true_tree, name) {
  if (!ape::is.rooted(true_tree)) stop("The true tree is not rooted.")
  rf <- numeric(length(trees))
  cherries <- integer(length(trees))
  for (i in seq_along(trees)) {
    tree <- trees[[i]]
    if (!setequal(tree$tip.label, true_tree$tip.label) ||
        anyDuplicated(tree$tip.label)) {
      stop(name, " tree ", i, " has tip labels different from the true tree.")
    }
    if (!ape::is.rooted(tree)) stop(name, " tree ", i, " is not rooted.")
    rf[i] <- as.numeric(phangorn::RF.dist(
      tree, true_tree, normalize = FALSE, check.labels = TRUE, rooted = TRUE
    ))
    cherries[i] <- TreeTools::Cherries(tree)
    if (i %% 10000L == 0L) cat(name, ": summarized ", i, " trees\n", sep = "")
  }
  data.frame(sample = name, tree_index = seq_along(trees),
             rf_distance = rf, n_cherries = cherries)
}

tree_clades <- function(tree) {
  n_tips <- length(tree$tip.label)
  children <- split(tree$edge[, 2L], tree$edge[, 1L])
  descendant_labels <- function(node) {
    if (node <= n_tips) return(tree$tip.label[node])
    unlist(lapply(children[[as.character(node)]], descendant_labels),
           use.names = FALSE)
  }
  keys <- lapply(seq.int(n_tips + 1L, n_tips + tree$Nnode), function(node) {
    labels <- descendant_labels(node)
    if (length(labels) < 2L || length(labels) == n_tips) return(NULL)
    paste(sort(labels), collapse = "|")
  })
  unlist(keys, use.names = FALSE)
}

ess_of <- function(values) {
  if (length(values) < 3L || anyNA(values) ||
      all(values == values[1L])) return(NA_real_)
  tryCatch(as.numeric(coda::effectiveSize(as.numeric(values))),
           error = function(e) NA_real_)
}

sample_chains <- function(a_data, beast_data, manifest_path) {
  if (!file.exists(manifest_path)) {
    stop("ESS needs sample_manifest.csv to identify the A-matrix chains: ",
         manifest_path)
  }
  manifest <- utils::read.csv(manifest_path)
  if (!all(c("chain_id", "draw_in_chain") %in% names(manifest)) ||
      nrow(manifest) != nrow(a_data) || anyNA(manifest$chain_id) ||
      anyNA(manifest$draw_in_chain)) {
    stop("sample_manifest.csv does not match the exported A-matrix trees.")
  }
  groups <- split(seq_len(nrow(manifest)), manifest$chain_id)
  if (any(vapply(groups, function(index) {
    !identical(as.integer(manifest$draw_in_chain[index]),
               seq_along(index))
  }, logical(1L)))) {
    stop("Within-chain A-matrix draws must be in sample order.")
  }
  chains <- lapply(names(groups), function(id) {
    list(sample = "A-matrix VI", chain = id, index = groups[[id]],
         data = a_data[groups[[id]], , drop = FALSE])
  })
  c(chains, list(list(sample = "BEAST", chain = "1",
                      index = seq_len(nrow(beast_data)), data = beast_data)))
}

calculate_ess <- function(a_trees, beast_trees, a_data, beast_data,
                          manifest_path) {
  chains <- sample_chains(a_data, beast_data, manifest_path)
  clades <- c(lapply(a_trees, tree_clades), lapply(beast_trees, tree_clades))
  all_keys <- sort(unique(unlist(clades, use.names = FALSE)))
  incidence <- matrix(0L, nrow = length(clades), ncol = length(all_keys),
                      dimnames = list(NULL, all_keys))
  for (i in seq_along(clades)) {
    if (length(clades[[i]])) incidence[i, match(clades[[i]], all_keys)] <- 1L
  }
  n_a <- length(a_trees)
  summaries <- vector("list", length(chains))
  clade_rows <- vector("list", length(chains))
  for (k in seq_along(chains)) {
    chain <- chains[[k]]
    row_index <- if (chain$sample == "BEAST") n_a + chain$index else chain$index
    indicators <- incidence[row_index, , drop = FALSE]
    clade_ess <- apply(indicators, 2L, ess_of)
    varying <- is.finite(clade_ess)
    n_draws <- length(row_index)
    summaries[[k]] <- data.frame(
      sample = chain$sample, chain = chain$chain, n_draws = n_draws,
      rf_ess = ess_of(chain$data$rf_distance),
      cherries_ess = ess_of(chain$data$n_cherries),
      min_clade_ess = if (any(varying)) min(clade_ess[varying]) else NA_real_,
      median_clade_ess = if (any(varying)) stats::median(clade_ess[varying]) else NA_real_,
      n_varying_clades = sum(varying)
    )
    clade_rows[[k]] <- data.frame(
      sample = chain$sample, chain = chain$chain, n_draws = n_draws,
      clade = all_keys, frequency = colMeans(indicators), ess = clade_ess
    )
  }
  list(summary = do.call(rbind, summaries),
       clades = do.call(rbind, clade_rows))
}

save_ess_plot <- function(ess_summary, path) {
  metrics <- c(rf_ess = "Rooted RF", cherries_ess = "Cherries",
               min_clade_ess = "Minimum varying clade")
  data <- do.call(rbind, lapply(names(metrics), function(metric) {
    data.frame(chain = paste(ess_summary$sample, ess_summary$chain),
               metric = metrics[[metric]],
               ess_per_draw = ess_summary[[metric]] / ess_summary$n_draws)
  }))
  data$metric <- factor(data$metric, levels = unname(metrics))
  plot <- ggplot(data, aes(x = chain, y = ess_per_draw, fill = chain)) +
    geom_col(na.rm = TRUE) + facet_wrap(~metric, nrow = 1L) +
    labs(title = "Tree-summary effective sample size",
         x = "Sample and chain", y = "ESS / retained draws") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  ggplot2::ggsave(path, plot = plot, width = 11, height = 5, dpi = 300)
}

save_histogram <- function(data, variable, x_label, title, path) {
  # Proportions make samples of different sizes directly comparable.
  counts <- as.data.frame(table(sample = data$sample, value = data[[variable]]),
                          stringsAsFactors = FALSE)
  counts <- counts[counts$Freq > 0L, , drop = FALSE]
  counts$value <- as.numeric(as.character(counts$value))
  counts$proportion <- counts$Freq / as.numeric(table(data$sample)[counts$sample])
  plot <- ggplot(counts, aes(x = value, y = proportion, fill = sample)) +
    geom_col(width = 0.9, color = "white") +
    facet_wrap(~sample, ncol = 1) +
    scale_x_continuous(breaks = sort(unique(counts$value))) +
    labs(title = title, x = x_label, y = "Proportion of trees") +
    theme_minimal(base_size = 13) + theme(legend.position = "none")
  ggplot2::ggsave(path, plot = plot, width = 8, height = 6, dpi = 300)
}

compare_tree_samples <- function(experiment, beast_path = NULL, true_path = NULL,
                                 plots_dir = NULL) {
  if (length(experiment) != 1L || !nzchar(experiment)) {
    stop("Pass an experiment ID or its exported output directory.")
  }
  output_dir <- if (dir.exists(experiment)) experiment else
    file.path(analysis_dir, "outputs", experiment)
  if (!dir.exists(output_dir)) stop("Cannot find exported experiment: ", output_dir)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  a_tree_path <- file.path(output_dir, "sampled_trees.newick")
  if (!file.exists(a_tree_path)) a_tree_path <- file.path(output_dir, "sampled_trees.trees")
  if (!file.exists(a_tree_path)) stop("No sampled_trees.newick or sampled_trees.trees in: ", output_dir)
  export_fasta <- file.path(output_dir, "sequences.fasta")
  if (!file.exists(export_fasta)) stop("No sequences.fasta in: ", output_dir)

  if (is.null(beast_path)) beast_path <- find_experiment_beast(output_dir)
  beast_path <- normalizePath(beast_path, mustWork = TRUE)
  if (is.null(true_path)) true_path <- find_true_tree(output_dir, export_fasta)
  true_path <- normalizePath(true_path, mustWork = TRUE)
  if (is.null(plots_dir)) {
    plots_dir <- file.path(analysis_dir, "plots", basename(output_dir))
  }

  true_tree <- ape::read.tree(true_path)
  if (!inherits(true_tree, "phylo")) stop("True-tree file must contain one tree.")
  labels <- names(ape::read.FASTA(export_fasta))
  if (!setequal(labels, true_tree$tip.label)) {
    stop("The true tree's taxon labels do not match the exported sequences.")
  }
  cat("A-matrix trees: ", a_tree_path, "\nBEAST trees: ", beast_path,
      "\nTrue tree: ", true_path, "\n", sep = "")
  a_trees <- read_tree_sample(a_tree_path)
  beast_trees <- read_tree_sample(beast_path)
  a_data <- summarize_sample(a_trees, true_tree, "A-matrix VI")
  beast_data <- summarize_sample(beast_trees, true_tree, "BEAST")
  ess <- calculate_ess(a_trees, beast_trees, a_data, beast_data,
                       file.path(output_dir, "sample_manifest.csv"))
  data <- rbind(a_data, beast_data)

  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  plots_dir <- normalizePath(plots_dir, mustWork = TRUE)
  utils::write.csv(data, file.path(plots_dir, "tree_summary.csv"), row.names = FALSE)
  utils::write.csv(ess$summary, file.path(plots_dir, "ess_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(ess$clades, file.path(plots_dir, "clade_ess.csv"),
                   row.names = FALSE)
  save_histogram(data, "rf_distance", "Rooted RF distance",
                 "Rooted RF distance from true tree",
                 file.path(plots_dir, "rooted_rf_vs_beast_rf_distance_histogram.png"))
  save_histogram(data, "n_cherries", "Number of cherries",
                 "Number of cherries",
                 file.path(plots_dir, "rooted_rf_vs_beast_cherry_histogram.png"))
  save_ess_plot(ess$summary, file.path(plots_dir, "ess_per_draw.png"))
  cat("Wrote plots, tree_summary.csv, ess_summary.csv, and clade_ess.csv to: ",
      plots_dir, "\n", sep = "")
  invisible(list(trees = data, ess = ess))
}

if (sys.nframe() == 0L || compare_console_run) {
  args <- if (compare_console_run) character() else commandArgs(trailingOnly = TRUE)
  if (length(args) > 4L) {
    stop("Usage: Rscript compare_tree_samples.R [EXPERIMENT_ID_OR_OUTPUT_DIR] [BEAST_TREES] [TRUE_TREE] [PLOTS_DIR]")
  }
  selected_experiment <- if (length(args)) args[1L] else
    if (!is.null(experiment_id)) experiment_id else default_experiment()
  compare_tree_samples(selected_experiment,
    beast_path = if (length(args) >= 2L) args[2L] else NULL,
    true_path = if (length(args) >= 3L) args[3L] else NULL,
    plots_dir = if (length(args) >= 4L) args[4L] else NULL)
}
