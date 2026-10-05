suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggplot2)
  library(TreeTools)
})

# Compare two tree samples against the true tree:
#   1. RF Boltzmann trees from the estimated beta
#   2. BEAST trees
#
# The script makes two histograms:
#   1. rooted RF distance from the true tree
#   2. number of cherries

true_tree_file <- "rf_true_tree_10tip.newick"
rf_beta_tree_file <- "rooted_rf_boltzmann_estimated_beta_100000_thinned_trees.newick"
beast_tree_file <- "rf_boltzman_beast_trees-rf_sequences_10tip.trees"

output_summary_file <- "rooted_rf_vs_beast_tree_summary.csv"
output_rf_histogram <- "rooted_rf_vs_beast_rf_distance_histogram.png"
output_cherry_histogram <- "rooted_rf_vs_beast_cherry_histogram.png"

read_tree_sample <- function(file) {
  if (!file.exists(file)) {
    stop("Tree file does not exist: ", file)
  }

  first_lines <- readLines(file, n = 20, warn = FALSE)
  is_nexus <- any(grepl("#NEXUS", first_lines, ignore.case = TRUE))

  trees <- if (is_nexus) {
    read.nexus(file)
  } else {
    read.tree(file)
  }

  if (inherits(trees, "phylo")) {
    trees <- list(trees)
  }

  class(trees) <- "multiPhylo"
  trees
}

count_cherries <- function(tree) {
  Cherries(tree)
}

summarize_tree_sample <- function(trees, true_tree, sample_name) {
  rows <- vector("list", length(trees))

  if (!is.rooted(true_tree)) {
    stop("The true tree must be rooted for rooted RF comparison.")
  }

  for (i in seq_along(trees)) {
    tree <- trees[[i]]

    if (!setequal(tree$tip.label, true_tree$tip.label)) {
      stop(sample_name, " tree ", i, " has tip labels that do not match the true tree.")
    }

    if (!is.rooted(tree)) {
      stop(sample_name, " tree ", i, " is not rooted.")
    }

    rows[[i]] <- data.frame(
      sample = sample_name,
      tree_index = i,
      rf_distance = as.numeric(RF.dist(
        tree,
        true_tree,
        normalize = FALSE,
        check.labels = TRUE,
        rooted = TRUE
      )),
      n_cherries = count_cherries(tree)
    )

    if (i %% 1000 == 0) {
      cat(sample_name, ": summarized", i, "trees\n")
    }
  }

  do.call(rbind, rows)
}

save_rf_histogram <- function(summary_data, file) {
  plot <- ggplot(summary_data, aes(x = rf_distance, fill = sample)) +
    geom_histogram(binwidth = 1, boundary = -0.5, color = "white") +
    facet_wrap(~sample, ncol = 1) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
      title = "Rooted RF distance from true tree",
      x = "Rooted RF distance",
      y = "Number of trees"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")

  ggsave(file, plot = plot, width = 8, height = 6, dpi = 300)
}

save_cherry_histogram <- function(summary_data, file) {
  plot <- ggplot(summary_data, aes(x = n_cherries, fill = sample)) +
    geom_histogram(binwidth = 1, boundary = -0.5, color = "white") +
    facet_wrap(~sample, ncol = 1) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
      title = "Number of cherries",
      x = "Number of cherries",
      y = "Number of trees"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")

  ggsave(file, plot = plot, width = 8, height = 6, dpi = 300)
}

true_tree <- read.tree(true_tree_file)
rf_beta_trees <- read_tree_sample(rf_beta_tree_file)
beast_trees <- read_tree_sample(beast_tree_file)

rf_beta_summary <- summarize_tree_sample(
  trees = rf_beta_trees,
  true_tree = true_tree,
  sample_name = "RF Boltzmann"
)

beast_summary <- summarize_tree_sample(
  trees = beast_trees,
  true_tree = true_tree,
  sample_name = "BEAST"
)

summary_data <- rbind(rf_beta_summary, beast_summary)

write.csv(summary_data, output_summary_file, row.names = FALSE)
save_rf_histogram(summary_data, output_rf_histogram)
save_cherry_histogram(summary_data, output_cherry_histogram)

cat("Wrote summary to", output_summary_file, "\n")
cat("Wrote RF histogram to", output_rf_histogram, "\n")
cat("Wrote cherry histogram to", output_cherry_histogram, "\n")
