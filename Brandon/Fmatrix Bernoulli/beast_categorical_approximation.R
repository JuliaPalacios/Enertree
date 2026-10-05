library(ape)
library(fmatrix)
library(phylodyn)
setwd("~/Desktop/Research/Palacios Lab/Untitled/Enertree")
source("Brandon/Fmatrix Bernoulli/utils.R")
source("Brandon/gibbs_prior/utils.R")
source("Brandon/gibbs_prior/gibbs_prior_joint_standardized_helpers.R", chdir = TRUE)

gen_Fmat <- phylodyn:::gen_Fmat

experiment_id <- "beast_vs_categorical"
output_root <- file.path(
  "Brandon", "Fmatrix Bernoulli", "plots", "beast_categorical_approximation", experiment_id
)
figures_dir <- file.path(output_root, "figures")
summaries_dir <- file.path(output_root, "summaries")

dir.create("plots", showWarnings = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summaries_dir, recursive = TRUE, showWarnings = FALSE)

save_path_for_plot_function <- function(filename) {
  file.path(
    "..", "Brandon", "Fmatrix Bernoulli", "plots",
    "beast_categorical_approximation", experiment_id, "figures", filename
  )
}

tree_file_candidates <- file.path(
  c(
    file.path(
      "Brandon", "BEAST Gibbs Approximation", "inputs",
      "beast_trees_summarized-sequences.trees"
    ),
    file.path(
      "Brandon", "Fmatrix Bernoulli",
      c("beast_trees_summarized-sequences.trees", "HKY1000Gibbs-sequences.trees")
    )
  )
)
trees_path <- tree_file_candidates[file.exists(tree_file_candidates)][1]
if (is.na(trees_path)) {
  stop("Could not find a BEAST .trees file.")
}

true_tree_candidates <- c(
  file.path("Brandon", "BEAST Gibbs Approximation", "inputs", "true_tree.newick"),
  file.path("Brandon", "Fmatrix Bernoulli", "true_tree.newick")
)
true_tree_path <- true_tree_candidates[file.exists(true_tree_candidates)][1]
if (is.na(true_tree_path)) {
  stop("Could not find true tree file.")
}

true_tree <- read.tree(true_tree_path)
true_fmat <- gen_Fmat(true_tree, tol = 8)

trees_beast <- read.nexus(trees_path)
beast_fmats <- lapply(trees_beast, gen_Fmat, tol = 8)

cat("Read", length(beast_fmats), "BEAST trees from", trees_path, "\n")
cat("Read true tree from", true_tree_path, "\n")

fit_result <- fit_and_generate_from_fmats(
  F_mats = beast_fmats,
  generated_sample_size = length(beast_fmats),
  seed = 1,
  method = "BFGS"
)

beast_average_fmat <- fit_result$beast_mean$projected
generated_average_fmat <- fit_result$generated_mean$projected
generated_fmats <- fit_result$generated_fmats

fmat_res_true <- fmat_posterior_comparison_plot(
  chain_fmats = generated_fmats,
  beast_fmats = beast_fmats,
  true_fmat = true_fmat,
  threshold = 0.01,
  title = "Categorical Approximation: F-matrix posterior vs BEAST (true tree reference)",
  save_path = save_path_for_plot_function("fmat_posterior_vs_beast_true_reference.png")
)

fmat_res_beast_average <- fmat_posterior_comparison_plot(
  chain_fmats = generated_fmats,
  beast_fmats = beast_fmats,
  true_fmat = beast_average_fmat,
  threshold = 0.01,
  title = "Categorical Approximation: F-matrix posterior vs BEAST (BEAST average reference)",
  save_path = save_path_for_plot_function("fmat_posterior_vs_beast_average_reference.png")
)

plot_subsample_size <- min(500, length(beast_fmats), length(generated_fmats))
set.seed(11)
beast_idx <- sample(length(beast_fmats), plot_subsample_size)
set.seed(12)
generated_idx <- sample(length(generated_fmats), plot_subsample_size)

mds_res <- tree_MDS_comparison_plot(
  trees_gen = generated_fmats[generated_idx],
  trees_data = beast_fmats[beast_idx],
  true_tree = true_fmat,
  M_tree = beast_average_fmat,
  reference_trees = list(generated_average_fmat),
  reference_labels = c("Generated Avg"),
  true_tree_label = "True Tree",
  M_tree_label = "BEAST Avg",
  title = "Categorical Approximation: tree MDS vs BEAST",
  save_path = save_path_for_plot_function("tree_mds_vs_beast.png")
)

hist_res_true <- tree_histogram_comparison_plot(
  chain_gen = generated_fmats,
  chain_data = beast_fmats,
  M_true_fmat = true_fmat,
  title = "Categorical Approximation: distance to true tree",
  save_path = save_path_for_plot_function("tree_histogram_vs_true.png")
)

hist_res_beast_average <- tree_histogram_comparison_plot(
  chain_gen = generated_fmats,
  chain_data = beast_fmats,
  M_true_fmat = beast_average_fmat,
  title = "Categorical Approximation: distance to BEAST average",
  save_path = save_path_for_plot_function("tree_histogram_vs_beast_average.png")
)

summary_df <- cbind(
  data.frame(
    experiment = experiment_id,
    beast_tree_file = basename(trees_path),
    true_tree_file = basename(true_tree_path),
    plot_subsample_size = plot_subsample_size
  ),
  fit_result$summary,
  data.frame(
    fmat_l1_distance_vs_beast = fmat_res_true$l1_distance,
    beast_average_l2_to_true = distance_Fmat(beast_average_fmat, true_fmat, dist = "l2"),
    generated_average_l2_to_true = distance_Fmat(generated_average_fmat, true_fmat, dist = "l2"),
    generated_average_l2_to_beast_average = distance_Fmat(
      generated_average_fmat,
      beast_average_fmat,
      dist = "l2"
    )
  )
)

fitted_vectors_df <- flatten_categorical_vectors(fit_result$fitted_vectors)

write.csv(
  summary_df,
  file = file.path(summaries_dir, "experiment_summary.csv"),
  row.names = FALSE
)
saveRDS(
  summary_df,
  file = file.path(summaries_dir, "experiment_summary.rds")
)
write.csv(
  fitted_vectors_df,
  file = file.path(summaries_dir, "fitted_categorical_vectors.csv"),
  row.names = FALSE
)
saveRDS(
  fit_result$fitted_vectors,
  file = file.path(summaries_dir, "fitted_categorical_vectors.rds")
)
saveRDS(
  fmat_res_true,
  file = file.path(summaries_dir, "fmat_posterior_comparison_true_reference.rds")
)
saveRDS(
  fmat_res_beast_average,
  file = file.path(summaries_dir, "fmat_posterior_comparison_beast_average_reference.rds")
)
saveRDS(
  mds_res,
  file = file.path(summaries_dir, "tree_mds_vs_beast.rds")
)
saveRDS(
  hist_res_true,
  file = file.path(summaries_dir, "tree_histogram_vs_true.rds")
)
saveRDS(
  hist_res_beast_average,
  file = file.path(summaries_dir, "tree_histogram_vs_beast_average.rds")
)
saveRDS(
  list(
    true_fmat = true_fmat,
    beast_average_fmat = beast_average_fmat,
    generated_average_fmat = generated_average_fmat
  ),
  file = file.path(summaries_dir, "reference_fmats.rds")
)
saveRDS(
  fmat_res_true$df,
  file = file.path(summaries_dir, "fmat_posterior_comparison_true_reference_df.rds")
)
saveRDS(
  fmat_res_beast_average$df,
  file = file.path(summaries_dir, "fmat_posterior_comparison_beast_average_reference_df.rds")
)

print(summary_df)
