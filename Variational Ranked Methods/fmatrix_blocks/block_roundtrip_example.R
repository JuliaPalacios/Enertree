suppressPackageStartupMessages({
  library(ape)
})

if (file.exists("heterochronous_fmatrix_utils.R")) {
  source("heterochronous_fmatrix_utils.R")
} else {
  source("Variational Ranked Methods/fmatrix_blocks/heterochronous_fmatrix_utils.R")
}

set.seed(1)
n_tips <- 6
leaf_order <- paste0("tip_", seq_len(n_tips))

tree <- rcoal(n_tips, tip.label = leaf_order)
coal_times <- cumsum(coalescent.intervals(tree)$interval.length)

plot(tree, main = "Original Tree")

labeled_fmat <- gen_labeled_fmat(
  tree = tree,
  leaf_order = leaf_order
)

cat("Labeled F-matrix:\n")
print(labeled_fmat)

blocks <- split_f_matrix_blocks(labeled_fmat, include_C = TRUE)

cat("\nA block:\n")
print(blocks$A)

cat("\nB block:\n")
print(blocks$B)

cat("\nC block:\n")
print(blocks$C)

cat("\nSaved coalescent times:\n")
print(coal_times)

reconstructed_without_times <- tree_from_labeled_fmat(
  f_mat = labeled_fmat,
  tip_labels = leaf_order
)

plot(reconstructed_without_times, main = "Reconstructed, Default Times")

reconstructed_with_times <- tree_from_labeled_fmat(
  f_mat = labeled_fmat,
  tip_labels = leaf_order,
  coal_times = coal_times
)

plot(reconstructed_with_times, main = "Reconstructed, Saved Times")


A_mat <- labeled_fmat_to_Amat(labeled_fmat, leaf_order = leaf_order)
cat("\nA matrix:\n")
print(A_mat)

