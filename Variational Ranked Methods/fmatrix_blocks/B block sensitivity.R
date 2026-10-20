source("heterochronous_fmatrix_utils.R")

print_distance_summary <- function(original_tree,
                                   comparison_tree,
                                   comparison_name,
                                   original_name = "Original",
                                   leaf_order) {
  cat("\n", original_name, " vs ", comparison_name, " full labeled F-matrix L1 distance:\n", sep = "")
  print(dist_labeled_fmat(original_tree, comparison_tree, dist = "l1", leaf_order = leaf_order))

  cat(original_name, " vs ", comparison_name, " full labeled F-matrix dist = \"l2\" value:\n", sep = "")
  print(dist_labeled_fmat(original_tree, comparison_tree, dist = "l2", leaf_order = leaf_order))

  cat(original_name, " vs ", comparison_name, " B-block L1 distance:\n", sep = "")
  print(dist_labeled_fmat(original_tree, comparison_tree, dist = "l1", onlyB = TRUE, leaf_order = leaf_order))

  cat(original_name, " vs ", comparison_name, " B-block dist = \"l2\" value:\n", sep = "")
  print(dist_labeled_fmat(original_tree, comparison_tree, dist = "l2", onlyB = TRUE, leaf_order = leaf_order))

  original_B <- split_f_matrix_blocks(
    gen_labeled_fmat(original_tree, leaf_order = leaf_order)
  )$B
  comparison_B <- split_f_matrix_blocks(
    gen_labeled_fmat(comparison_tree, leaf_order = leaf_order)
  )$B
  cat(original_name, " vs ", comparison_name, " B-block difference:\n", sep = "")
  print(comparison_B - original_B)
}

swap_tip_labels <- function(tree, tip_a, tip_b) {
  swapped_tree <- tree
  tip_a_index <- match(tip_a, swapped_tree$tip.label)
  tip_b_index <- match(tip_b, swapped_tree$tip.label)

  if (is.na(tip_a_index) || is.na(tip_b_index)) {
    stop("Both tip labels must be present in tree.")
  }

  swapped_tree$tip.label[tip_a_index] <- tip_b
  swapped_tree$tip.label[tip_b_index] <- tip_a
  swapped_tree
}

move_mrca_between <- function(tree, target_tips, first_tips, second_tips) {
  moved_tree <- tree
  node_depths <- node.depth.edgelength(moved_tree)

  target_node <- getMRCA(moved_tree, target_tips)
  first_node <- getMRCA(moved_tree, first_tips)
  second_node <- getMRCA(moved_tree, second_tips)

  if (is.null(target_node) || is.null(first_node) || is.null(second_node)) {
    stop("Could not find one of the requested MRCA nodes.")
  }

  node_depths[target_node] <- mean(node_depths[c(first_node, second_node)])
  moved_tree$edge.length <- node_depths[moved_tree$edge[, 2]] -
    node_depths[moved_tree$edge[, 1]]

  if (any(moved_tree$edge.length <= 0)) {
    stop("Moving the MRCA created a non-positive edge length.")
  }

  moved_tree
}



n_tips <-10
leaf_order <- paste0("t", seq_len(n_tips))
set.seed(42)
tree <- rcoal(n_tips)

depths <- node.depth.edgelength(tree)
internal_nodes <- (n_tips + 1L):(n_tips + tree$Nnode)
tip_depth <- max(depths[seq_len(n_tips)])
internal_ages <- tip_depth - depths[internal_nodes]
ranked_internal_nodes <- internal_nodes[order(internal_ages)]
new_node_depths <- numeric(length(depths))
new_node_depths[seq_len(n_tips)] <- n_tips - 1L
new_node_depths[ranked_internal_nodes] <- rev(seq_len(tree$Nnode) - 1L)
tree$edge.length <- new_node_depths[tree$edge[, 2]] -
  new_node_depths[tree$edge[, 1]]

plot(tree, main = "Original tree: evenly spaced coalescent events")


cat("Original tree tip labels:\n")
print(tree$tip.label)

F <- gen_labeled_fmat(tree, leaf_order = leaf_order)
blocks <- split_f_matrix_blocks(F)


cherry_augment_tree <- swap_tip_labels(tree, "t4", "t7")
cat("Cherry-swap tree tip labels:\n")
print(cherry_augment_tree$tip.label)
plot(cherry_augment_tree, main = "Cherry tip-label swap: t4 <-> t7")


print_distance_summary(tree, cherry_augment_tree, "cherry-swap tree", leaf_order = leaf_order)



small_augment_tree <- swap_tip_labels(tree, "t6", "t4")
cat("Small-swap tree tip labels:\n")
print(small_augment_tree$tip.label)
plot(small_augment_tree, main = "Small tip-label swap: t6 <-> t4")

print_distance_summary(tree, small_augment_tree, "small-swap tree", leaf_order = leaf_order)


same_small_swap_problem_tree <- swap_tip_labels(tree, "t6", "t7")
cat("Small-swap tree tip labels:\n")
print(same_small_swap_problem_tree$tip.label)
plot(same_small_swap_problem_tree, main = "Small tip-label swap: t6 <-> t7")

print_distance_summary(tree, small_augment_tree, "small-swap tree", leaf_order = leaf_order)

large_augment_tree <- swap_tip_labels(tree, "t9", "t6")
cat("Large-swap tree tip labels:\n")
print(large_augment_tree$tip.label)
plot(large_augment_tree, main = "Large tip-label swap: t9 <-> t6")

print_distance_summary(tree, large_augment_tree, "large-swap tree", leaf_order = leaf_order)



farthest_augment_tree <- swap_tip_labels(tree, "t8", "t7")
cat("Largest-swap tree tip labels:\n")
print(farthest_augment_tree$tip.label)
plot(farthest_augment_tree, main = "Largest tip-label swap: t8 <-> t7")

print_distance_summary(tree, farthest_augment_tree, "largest-swap tree", leaf_order = leaf_order)



cat("\nSecond example: first tree with t3/t8 coalescence between t7/t4 and t6/t5\n")

second_tree <- move_mrca_between(
  tree,
  target_tips = c("t3", "t8"),
  first_tips = c("t7", "t4"),
  second_tips = c("t6", "t5")
)

plot(second_tree, main = "Second Tree")

cat("Second original tree tip labels:\n")
print(second_tree$tip.label)


second_small_augment_tree <- swap_tip_labels(second_tree, "t6", "t7")
cat("Second small-swap tree tip labels:\n")
print(second_small_augment_tree$tip.label)
plot(second_small_augment_tree, main = "Second small tip-label swap: t6 <-> t7")

print_distance_summary(
  second_tree,
  second_small_augment_tree,
  "second small-swap tree",
  original_name = "Second original tree",
  leaf_order = leaf_order
)

second_largest_augment_tree <- swap_tip_labels(second_tree, "t8", "t7")
cat("Second largest-swap tree tip labels:\n")
print(second_largest_augment_tree$tip.label)
plot(second_largest_augment_tree, main = "Second largest tip-label swap: t8 <-> t7")

print_distance_summary(
  second_tree,
  second_largest_augment_tree,
  "second largest-swap tree",
  original_name = "Second original tree",
  leaf_order = leaf_order
)
