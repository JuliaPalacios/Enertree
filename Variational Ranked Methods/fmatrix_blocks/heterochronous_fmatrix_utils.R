suppressPackageStartupMessages({
  library(ape)
})

gen_labeled_fmat <- function(tree,
                             leaf_order = tree$tip.label,
                             epsilon = NULL,
                             tol = 13) {
  if (!requireNamespace("phylodyn", quietly = TRUE)) {
    stop("The phylodyn package is required.")
  }
  if (!is.rooted(tree)) {
    stop("tree must be rooted.")
  }
  if (is.null(tree$edge.length)) {
    stop("tree must have branch lengths.")
  }

  n_tips <- length(tree$tip.label)
  if (length(leaf_order) != n_tips || !setequal(leaf_order, tree$tip.label)) {
    stop("leaf_order must contain exactly the tree tip labels.")
  }

  depths <- node.depth.edgelength(tree)
  leaf_depths <- depths[seq_len(n_tips)]
  internal_depths <- depths[(n_tips + 1L):(n_tips + tree$Nnode)]

  if (is.null(epsilon)) {
    internal_times_before_present <- max(leaf_depths) - internal_depths
    positive_internal_times <- internal_times_before_present[
      internal_times_before_present > 0
    ]
    if (length(positive_internal_times) == 0L) {
      stop("Could not choose epsilon from tree depths.")
    }
    epsilon <- min(positive_internal_times) / (10 * n_tips)
  }

  leaf_rank <- setNames(seq_along(leaf_order), leaf_order)
  leaf_edges <- match(seq_len(n_tips), tree$edge[, 2])
  if (any(is.na(leaf_edges))) {
    stop("Could not find every leaf edge.")
  }

  augmented_tree <- tree
  augmented_tree$edge.length[leaf_edges] <-
    augmented_tree$edge.length[leaf_edges] +
    epsilon * leaf_rank[tree$tip.label]

  phylodyn:::gen_Fmat(augmented_tree, tol = tol)
}

tree_from_labeled_fmat <- function(f_mat,
                                   tip_labels = NULL,
                                   coal_times = NULL,
                                   coalescent_spacing = 1) {
  if (!requireNamespace("fmatrix", quietly = TRUE)) {
    stop("The fmatrix package is required.")
  }

  blocks <- split_f_matrix_blocks(f_mat)
  A <- blocks$A
  B <- blocks$B
  n_tips <- nrow(A) + 1L

  if (is.null(tip_labels)) {
    tip_labels <- paste0("t", seq_len(n_tips))
  }
  if (length(tip_labels) != n_tips) {
    stop("tip_labels must have one label for each tip.")
  }

  leaf_parent_ranks <- leaf_parent_ranks_from_B(A, B)
  names(leaf_parent_ranks) <- tip_labels

  if (is.null(coal_times)) {
    coal_times <- seq_len(n_tips - 1L) * coalescent_spacing
  } else {
    if (length(coal_times) != n_tips - 1L) {
      stop("coal_times must have length n_tips - 1.")
    }
    if (is.unsorted(coal_times, strictly = TRUE)) {
      stop("coal_times must be strictly increasing.")
    }
  }

  tree <- fmatrix::mytree_from_F(A, coal.times = coal_times)
  tree <- reorder(tree, "postorder")

  n_tree_tips <- length(tree$tip.label)
  internal_nodes <- (n_tree_tips + 1L):(n_tree_tips + tree$Nnode)
  depths <- node.depth.edgelength(tree)
  ranked_nodes <- internal_nodes[order(depths[internal_nodes], internal_nodes)]
  internal_ranks <- setNames(seq_along(ranked_nodes) - 1L, ranked_nodes)

  tip_parent_ranks <- integer(n_tree_tips)
  for (edge_index in seq_len(nrow(tree$edge))) {
    child <- tree$edge[edge_index, 2]
    if (child <= n_tree_tips) {
      parent <- tree$edge[edge_index, 1]
      tip_parent_ranks[child] <- internal_ranks[[as.character(parent)]]
    }
  }

  available_tips <- split(seq_len(n_tree_tips), tip_parent_ranks)
  new_tip_labels <- rep(NA_character_, n_tree_tips)

  for (label in tip_labels) {
    parent_rank <- as.character(leaf_parent_ranks[[label]])
    candidate_tips <- available_tips[[parent_rank]]
    if (length(candidate_tips) == 0L) {
      stop("No reconstructed tip is available for parent rank ", parent_rank, ".")
    }

    tip_node <- candidate_tips[[1L]]
    new_tip_labels[[tip_node]] <- label
    available_tips[[parent_rank]] <- candidate_tips[-1L]
  }

  if (anyNA(new_tip_labels)) {
    stop("Could not assign every tip label.")
  }

  tree$tip.label <- new_tip_labels
  reorder(tree, "postorder")
}

labeled_fmat_to_Amat <- function(f_mat,
                                  leaf_order,
                                  coal_times = NULL,
                                  coalescent_spacing = 1) {
  if (!requireNamespace("phylodyn", quietly = TRUE)) {
    stop("The phylodyn package is required.")
  }

  tree <- tree_from_labeled_fmat(
    f_mat = f_mat,
    tip_labels = leaf_order,
    coal_times = coal_times,
    coalescent_spacing = coalescent_spacing
  )
  A_phy <- phylodyn:::gen_Amat(tree)
  A_phy[, match(leaf_order, tree$tip.label), drop = FALSE]
}

tree_matrix_component_distances <- function(tree1,
                                            tree2,
                                            matrices = c("D", "F_A", "Amat"),
                                            norm = c("l1", "l2"),
                                            leaf_order = NULL) {
  norm <- match.arg(norm)
  tree_leaf_order <- if (is.null(leaf_order)) sort(tree1$tip.label) else leaf_order

  matrices1 <- list()
  matrices2 <- list()

  if (any(matrices %in% c("D", "F_A"))) {
    f_mat1 <- gen_labeled_fmat(tree1, leaf_order = tree_leaf_order)
    f_mat2 <- gen_labeled_fmat(tree2, leaf_order = tree_leaf_order)
    blocks1 <- split_f_matrix_blocks(f_mat1)
    blocks2 <- split_f_matrix_blocks(f_mat2)
    if ("D" %in% matrices) {
      matrices1$D <- make_D_matrix(blocks1$A, blocks1$B)
      matrices2$D <- make_D_matrix(blocks2$A, blocks2$B)
    }
    if ("F_A" %in% matrices) {
      matrices1$F_A <- blocks1$A
      matrices2$F_A <- blocks2$A
    }
  }

  if ("Amat" %in% matrices) {
    A1 <- phylodyn:::gen_Amat(tree1)
    A2 <- phylodyn:::gen_Amat(tree2)
    matrices1$Amat <- A1[, match(tree_leaf_order, tree1$tip.label), drop = FALSE]
    matrices2$Amat <- A2[, match(tree_leaf_order, tree2$tip.label), drop = FALSE]
  }

  vapply(matrices, function(matrix_name) {
    difference <- matrices2[[matrix_name]] - matrices1[[matrix_name]]
    if (norm == "l1") sum(abs(difference)) else sqrt(sum(difference^2))
  }, numeric(1))
}

tree_matrix_distance <- function(matrices = "D",
                                 weights = rep(1, length(matrices)),
                                 norm = c("l1", "l2"),
                                 leaf_order = NULL) {
  norm <- match.arg(norm)

  function(tree1, tree2) {
    distances <- tree_matrix_component_distances(
      tree1,
      tree2,
      matrices = matrices,
      norm = norm,
      leaf_order = leaf_order
    )
    sum(weights * distances)
  }
}

make_deterministic_C_block <- function(n_tips) {
  c_block <- matrix(0L, nrow = n_tips - 1L, ncol = n_tips - 1L)

  for (row in seq_len(n_tips - 1L)) {
    c_block[row, seq_len(row)] <- n_tips - row
  }

  c_block
}

split_f_matrix_blocks <- function(f_mat, include_C = FALSE) {
  m <- nrow(f_mat)
  if (!all(dim(f_mat) == c(m, m))) {
    stop("f_mat must be square.")
  }
  if (m %% 2L != 0L) {
    stop("f_mat must have even dimension 2 * n_tips - 2.")
  }

  n_tips <- as.integer(m / 2L + 1L)
  internal <- seq_len(n_tips - 1L)
  sampled <- n_tips:m

  blocks <- list(
    A = f_mat[internal, internal, drop = FALSE],
    B = f_mat[sampled, internal, drop = FALSE]
  )

  if (include_C) {
    blocks$C <- f_mat[sampled, sampled, drop = FALSE]
  }

  blocks
}

combine_A_B_blocks <- function(A, B) {
  n_tips <- nrow(A) + 1L
  if (!all(dim(A) == c(n_tips - 1L, n_tips - 1L))) {
    stop("A must be square.")
  }
  if (!all(dim(B) == c(n_tips - 1L, n_tips - 1L))) {
    stop("B must have the same dimensions as A.")
  }

  top_right <- matrix(0L, nrow = n_tips - 1L, ncol = n_tips - 1L)
  C <- make_deterministic_C_block(n_tips)

  rbind(
    cbind(A, top_right),
    cbind(B, C)
  )
}

dist_labeled_fmat <- function(tree1,
                              tree2,
                              dist = c("l2", "l1"),
                              onlyB = FALSE,
                              leaf_order = NULL,
                              ...) {
  dist <- match.arg(dist)
  
  if (is.null(leaf_order)) {
    if (!setequal(tree1$tip.label, tree2$tip.label)) {
      stop("tree1 and tree2 must have the same tip labels.")
    }
    leaf_order <- sort(tree1$tip.label)
  }

  F1 <- gen_labeled_fmat(tree1, leaf_order = leaf_order, ...)
  F2 <- gen_labeled_fmat(tree2, leaf_order = leaf_order, ...)

  if (onlyB) {
    F1 <- split_f_matrix_blocks(F1)$B
    F2 <- split_f_matrix_blocks(F2)$B
  }

  delta <- F2 - F1
  if (dist == "l1") {
    return(sum(abs(delta)))
  }

  d2 <- sum(delta^2)
  return(d2)
}

previous_d_counts <- function(previous_row) {
  if (length(previous_row) < 1L) {
    stop("previous_row must contain at least one entry.")
  }

  as.integer(c(previous_row[1L], diff(previous_row)))
}

matrix_key <- function(f_mat) {
  paste(f_mat[lower.tri(f_mat, diag = TRUE)], collapse = ",")
}

# Finds the parent rank from a leaf-change string; e.g. c(0, 0, 1, 1) -> 2.
find_parent_rank_from_leaf <- function(leaf_change) {
  leaf_change <- as.integer(leaf_change)
  first_one <- match(1L, leaf_change)

  if (is.na(first_one)) {
    stop("Expected at least one 1.")
  }

  expected <- as.integer(seq_along(leaf_change) >= first_one)
  if (!identical(leaf_change, expected)) {
    stop("Expected zeros followed by ones.")
  }

  first_one - 1L
}

leaf_parent_ranks_from_B <- function(A, B) {
  apply(make_D_matrix(A, B), 1, find_parent_rank_from_leaf)
}

make_D_matrix <- function(A, B) {
  stacked_rows <- rbind(
    A[nrow(A), , drop = FALSE],
    B,
    matrix(0L, nrow = 1L, ncol = ncol(B))
  )

  stacked_rows[-nrow(stacked_rows), , drop = FALSE] -
    stacked_rows[-1L, , drop = FALSE]
}

fixed_isochronous_completion_diagonal <- function(n_tips) {
  if (length(n_tips) != 1L || n_tips < 2L || n_tips != as.integer(n_tips)) {
    stop("n_tips must be a single integer at least 2.")
  }

  c(seq.int(2L, n_tips), seq.int(n_tips - 1L, 1L))
}

sample_fixed_diagonal_f_matrix <- function(diagonal) {
  if (length(diagonal) < 2L || any(diagonal <= 0L)) {
    stop("diagonal must contain at least two positive entries.")
  }

  m <- length(diagonal)
  fmat <- matrix(0L, nrow = m, ncol = m)
  diag(fmat) <- diagonal

  for (i in 2:m) {
    cols <- seq_len(i - 1L)
    previous_row <- fmat[i - 1L, cols]
    D <- previous_d_counts(previous_row)
    L <- sample.int(i - 1L, size = 1L, prob = D)

    fmat[i, cols] <- previous_row - as.integer(cols >= L)
  }

  fmat
}

uniform_f_matrix <- function(n_tips) {
  sample_fixed_diagonal_f_matrix(
    fixed_isochronous_completion_diagonal(n_tips)
  )
}

enumerate_fixed_diagonal_fmats <- function(diagonal) {
  m <- length(diagonal)
  
  fmat <- matrix(0L, nrow = m, ncol = m)
  diag(fmat) <- diagonal
  
  states <- list(list(F = fmat, probability = 1))
  
  for (i in 2:m) {
    new_states <- list()
    
    for (state in states) {
      fmat <- state$F
      prob <- state$probability
      
      cols <- seq_len(i - 1L)
      previous_row <- fmat[i - 1L, cols]
      D <- previous_d_counts(previous_row)
      
      for (L in which(D > 0L)) {
        new_fmat <- fmat
        new_fmat[i, cols] <- previous_row - as.integer(cols >= L)
        
        new_prob <- prob * D[L] / sum(D)
        
        new_states[[length(new_states) + 1L]] <- list(
          F = new_fmat,
          probability = new_prob
        )
      }
    }
    
    states <- new_states
  }
  
  keys <- vapply(states, function(state) matrix_key(state$F), character(1))
  if (any(duplicated(keys))) {
    stop("Enumeration produced duplicate F-matrices.")
  }
  
  states
}
