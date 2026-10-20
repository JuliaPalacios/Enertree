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
source(file.path(script_dir, "swap_utils.R"))

generate_base_trees <- function(n_tips,
                                n_iterations,
                                tree_mode = c("same_tree", "new_tree_each_iteration"),
                                leaf_order = paste0("t", seq_len(n_tips))) {
  tree_mode <- match.arg(tree_mode)
  generate_tree <- function() ape::rcoal(n_tips, tip.label = leaf_order)

  if (tree_mode == "same_tree") {
    return(rep(list(generate_tree()), n_iterations))
  }

  lapply(seq_len(n_iterations), function(iteration) generate_tree())
}

simulate_swap_distances <- function(n_tips,
                                    n_iterations,
                                    n_swaps,
                                    tree_mode = c("same_tree", "new_tree_each_iteration"),
                                    swap_fn,
                                    distance_fns = list(D = tree_matrix_distance()),
                                    leaf_order = paste0("t", seq_len(n_tips)),
                                    seed = NULL,
                                    keep_trees = FALSE,
                                    base_trees = NULL,
                                    progress_fn = NULL) {
  tree_mode <- match.arg(tree_mode)
  n_swaps <- sort(unique(n_swaps))

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(base_trees)) {
    base_trees <- generate_base_trees(
      n_tips = n_tips,
      n_iterations = n_iterations,
      tree_mode = tree_mode,
      leaf_order = leaf_order
    )
  }

  shared_tree <- if (tree_mode == "same_tree") base_trees[[1L]] else NULL

  distance_rows <- vector("list", n_iterations * length(n_swaps))
  kept_base_trees <- if (keep_trees) vector("list", n_iterations) else NULL
  swapped_trees <- if (keep_trees) vector("list", n_iterations) else NULL
  row_index <- 1L

  for (iteration in seq_len(n_iterations)) {
    base_tree <- base_trees[[iteration]]
    swapped_tree <- base_tree
    applied_swaps <- 0L

    for (target_n_swaps in n_swaps) {
      swapped_tree <- apply_swap_n_times(
        swapped_tree,
        swap_fn,
        target_n_swaps - applied_swaps
      )
      applied_swaps <- target_n_swaps

      distances <- unlist(lapply(
        distance_fns,
        function(distance_fn) distance_fn(base_tree, swapped_tree)
      ), use.names = TRUE)
      distance_rows[[row_index]] <- data.frame(
        iteration = iteration,
        n_swaps = target_n_swaps,
        as.list(distances),
        check.names = FALSE
      )
      row_index <- row_index + 1L
      if (!is.null(progress_fn)) progress_fn()
    }

    if (keep_trees) {
      kept_base_trees[[iteration]] <- base_tree
      swapped_trees[[iteration]] <- swapped_tree
    }
  }

  list(
    distances = do.call(rbind, distance_rows),
    tree_mode = tree_mode,
    n_tips = n_tips,
    n_iterations = n_iterations,
    n_swaps = n_swaps,
    leaf_order = leaf_order,
    shared_tree = shared_tree,
    base_trees = kept_base_trees,
    swapped_trees = swapped_trees
  )
}
