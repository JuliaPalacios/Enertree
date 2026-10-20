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

random_tip_label_swap <- function(tree) {
  selected_tips <- sample(tree$tip.label, size = 2L, replace = FALSE)
  swap_tip_labels(tree, selected_tips[[1L]], selected_tips[[2L]])
}

rank_tree_topology <- function(tree) {
  n_tips <- length(tree$tip.label)
  internal_nodes <- (n_tips + 1L):(n_tips + tree$Nnode)
  topology_tree <- tree
  topology_tree$edge.length <- rep(1, nrow(topology_tree$edge))
  depths <- ape::node.depth.edgelength(topology_tree)
  ordered_internal_nodes <- internal_nodes[order(depths[internal_nodes], internal_nodes)]

  node_heights <- numeric(n_tips + tree$Nnode)
  node_heights[ordered_internal_nodes] <- seq.int(0, n_tips - 2L)
  node_heights[seq_len(n_tips)] <- n_tips - 1L

  ranked_tree <- tree
  ranked_tree$edge.length <- node_heights[ranked_tree$edge[, 2]] -
    node_heights[ranked_tree$edge[, 1]]
  ranked_tree
}

random_nni_swap <- function(tree) {
  proposal <- phangorn::rNNI(tree, moves = 1L)
  rank_tree_topology(proposal)
}

apply_swap_n_times <- function(tree, swap_fn, n_swaps) {
  for (swap_index in seq_len(n_swaps)) {
    tree <- swap_fn(tree)
  }
  tree
}
