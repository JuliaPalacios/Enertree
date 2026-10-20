
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

block_key <- function(M) {
  paste(as.integer(M), collapse = ",")
}

check_n_tips <- function(n_tips) {
  diagonal <- fixed_isochronous_completion_diagonal(n_tips)
  states <- enumerate_fixed_diagonal_fmats(diagonal)
  F_list <- lapply(states, `[[`, "F")

  AB_keys <- character(length(F_list))
  A_parent_rank_keys <- character(length(F_list))
  roundtrip_ok <- logical(length(F_list))

  for (i in seq_along(F_list)) {
    F_mat <- F_list[[i]]
    blocks <- split_f_matrix_blocks(F_mat)

    roundtrip_ok[[i]] <- identical(
      F_mat,
      combine_A_B_blocks(blocks$A, blocks$B)
    )

    AB_keys[[i]] <- paste(
      block_key(blocks$A),
      block_key(blocks$B),
      sep = "||"
    )

    parent_ranks <- leaf_parent_ranks_from_B(blocks$A, blocks$B)
    A_parent_rank_keys[[i]] <- paste(
      block_key(blocks$A),
      paste(parent_ranks, collapse = ","),
      sep = "||"
    )
  }

  data.frame(
    n_tips = n_tips,
    full_F_count = length(F_list),
    AB_roundtrip_all = all(roundtrip_ok),
    unique_AB_count = length(unique(AB_keys)),
    duplicate_AB = any(duplicated(AB_keys)),
    unique_A_parent_rank_count = length(unique(A_parent_rank_keys)),
    duplicate_A_parent_rank = any(duplicated(A_parent_rank_keys))
  )
}

results <- do.call(rbind, lapply(3:7, check_n_tips))
print(results)
