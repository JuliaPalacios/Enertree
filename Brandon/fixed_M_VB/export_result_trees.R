## Export final-beta tree samples from saved fixed-M VB results.
##
## This does not rerun the optimizer. It reads result_*.rds files that contain
## final_beta, M_fixed, and coal_times, then samples trees from the exact prior.

suppressPackageStartupMessages({
  library(ape)
  library(phylodyn)
  library(fmatrix)
})

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  getwd()
})

source(file.path(script_dir, "VB_Utils.R"))

env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (value == "") {
    return(default)
  }
  as.integer(value)
}

env_character_vector <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (value == "") {
    return(default)
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

default_result_paths <- Sys.glob(file.path(
  script_dir,
  "logs",
  "VB_fixed_M_exact_prior",
  "*",
  "result_*.rds"
))

result_paths <- env_character_vector("VB_RESULT_PATHS", default_result_paths)
num_trees <- env_integer("VB_EXPORT_NUM_TREES", 100000L)
chunk_size <- env_integer("VB_EXPORT_CHUNK_SIZE", 5000)
export_seed <- env_integer("VB_EXPORT_SEED", 42)

if (length(result_paths) == 0) {
  stop("No result_*.rds files found. Set VB_RESULT_PATHS to one or more result files.")
}

write_reconstructed_tree <- function(Fmat, coal_times, file, seed) {
  set.seed(seed)
  tree <- mytree_from_F(Fmat, coal_times)
  write.tree(tree, file = file)
  invisible(tree)
}

format_count <- function(x) {
  format(x, scientific = FALSE, trim = TRUE)
}

write_tree_or_reconstruction <- function(result, tree_name, Fmat, coal_times, file, seed) {
  if (!is.null(result[[tree_name]])) {
    write.tree(result[[tree_name]], file = file)
    return("saved")
  }

  write_reconstructed_tree(Fmat, coal_times, file, seed)
  "reconstructed"
}

export_one_result <- function(result_path, result_index) {
  result <- readRDS(result_path)
  required_fields <- c("fixed_M_source", "M_fixed", "M_true", "M_upgma", "coal_times", "final_beta")
  missing_fields <- setdiff(required_fields, names(result))
  if (length(missing_fields) > 0) {
    stop(
      "Result file is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      "\nFile: ",
      result_path
    )
  }

  num_tips <- nrow(result$M_fixed) + 1L
  all_Fmats <- get(paste0("F.list", num_tips))
  exact_prior_cache <- phylodyn:::precompute_tree_chain_distance_cache(all_Fmats)

  result_id <- sub("\\.rds$", "", sub("^result_", "", basename(result_path)))
  export_dir <- file.path(dirname(result_path), paste0("final_beta_tree_export_", result_id))
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

  fixed_M_source <- result$fixed_M_source
  prefix <- paste0(fixed_M_source, "_", result_id)

  sequences_fasta_path <- file.path(export_dir, paste0("sequences_", prefix, ".fasta"))
  fixed_tree_path <- file.path(export_dir, paste0("fixed_tree_", prefix, ".newick"))
  true_tree_path <- file.path(export_dir, paste0("true_tree_", prefix, ".newick"))
  upgma_tree_path <- file.path(export_dir, paste0("upgma_tree_", prefix, ".newick"))
  num_trees_label <- format_count(num_trees)
  generated_trees_path <- file.path(
    export_dir,
    paste0("generated_", num_trees_label, "_trees_final_beta_fixed_M_", prefix, ".trees")
  )
  generated_indices_path <- file.path(
    export_dir,
    paste0("generated_", num_trees_label, "_tree_indices_final_beta_fixed_M_", prefix, ".rds")
  )

  write.csv(result$M_fixed, file.path(export_dir, paste0("M_fixed_", prefix, ".csv")), row.names = FALSE)
  write.csv(result$M_true, file.path(export_dir, paste0("M_true_", prefix, ".csv")), row.names = FALSE)
  write.csv(result$M_upgma, file.path(export_dir, paste0("M_upgma_", prefix, ".csv")), row.names = FALSE)

  sequences_status <- "missing"
  if (!is.null(result$sequences)) {
    phangorn::write.phyDat(result$sequences, file = sequences_fasta_path, format = "fasta")
    sequences_status <- "saved"
  }

  fixed_tree_status <- write_tree_or_reconstruction(
    result, "fixed_tree", result$M_fixed, result$coal_times,
    fixed_tree_path, export_seed + 1000L * result_index + 1L
  )
  true_tree_status <- write_tree_or_reconstruction(
    result, "M_true_tree", result$M_true, result$coal_times,
    true_tree_path, export_seed + 1000L * result_index + 2L
  )
  upgma_tree_status <- write_tree_or_reconstruction(
    result, "upgma_tree", result$M_upgma, result$coal_times,
    upgma_tree_path, export_seed + 1000L * result_index + 3L
  )

  cat("Generating", num_trees, "trees for", result_path, "\n")
  generated_info <- write_prior_trees_to_file(
    M = result$M_fixed,
    b = result$final_beta,
    num_trees = num_trees,
    all_Fmats = all_Fmats,
    coal_times = result$coal_times,
    file = generated_trees_path,
    cache = exact_prior_cache,
    chunk_size = chunk_size,
    seed = export_seed + 1000L * result_index + 10L,
    indices_file = generated_indices_path
  )

  note_path <- file.path(export_dir, paste0("README_", prefix, ".txt"))
  writeLines(
    c(
      "This export was generated from a saved fixed-M VB result RDS.",
      "The optimizer was not rerun.",
      "",
      paste("Sequences FASTA status:", sequences_status),
      paste("Fixed tree status:", fixed_tree_status),
      paste("True tree status:", true_tree_status),
      paste("UPGMA tree status:", upgma_tree_status),
      "Status 'reconstructed' means the file was rebuilt from saved F-matrices and coal_times.",
      "",
      paste("Source result:", normalizePath(result_path, mustWork = FALSE)),
      paste("Fixed M source:", fixed_M_source),
      paste("Final beta:", signif(result$final_beta, 16)),
      paste("Generated tree count:", num_trees),
      paste("Export seed:", export_seed)
    ),
    con = note_path
  )

  summary <- data.frame(
    result_path = normalizePath(result_path, mustWork = FALSE),
    export_dir = normalizePath(export_dir, mustWork = FALSE),
    fixed_M_source = fixed_M_source,
    num_tips = num_tips,
    final_beta = result$final_beta,
    final_log_beta = result$final_log_beta,
    num_generated_trees = num_trees,
    chunk_size = chunk_size,
    export_seed = export_seed,
    exact_prior_log_Z = generated_info$log_Z,
    exact_prior_prob_sum = generated_info$prob_sum,
    exact_prior_ess = generated_info$ess,
    exact_prior_support_size = generated_info$support_size,
    unique_generated_fmats = generated_info$unique_fmats,
    generated_trees_path = normalizePath(generated_trees_path, mustWork = FALSE),
    generated_indices_path = normalizePath(generated_indices_path, mustWork = FALSE),
    sequences_fasta_path = if (file.exists(sequences_fasta_path)) normalizePath(sequences_fasta_path, mustWork = FALSE) else NA_character_,
    sequences_status = sequences_status,
    fixed_tree_path = normalizePath(fixed_tree_path, mustWork = FALSE),
    fixed_tree_status = fixed_tree_status,
    true_tree_path = normalizePath(true_tree_path, mustWork = FALSE),
    true_tree_status = true_tree_status,
    upgma_tree_path = normalizePath(upgma_tree_path, mustWork = FALSE),
    upgma_tree_status = upgma_tree_status,
    note_path = normalizePath(note_path, mustWork = FALSE)
  )

  summary_path <- file.path(export_dir, paste0("export_summary_", prefix, ".csv"))
  write.csv(summary, summary_path, row.names = FALSE)
  saveRDS(
    list(summary = summary, generated_info = generated_info, source_result = result),
    file = file.path(export_dir, paste0("export_summary_", prefix, ".rds"))
  )

  summary
}

summaries <- do.call(
  rbind,
  Map(export_one_result, result_paths, seq_along(result_paths))
)

combined_summary_path <- file.path(
  script_dir,
  "logs",
  "VB_fixed_M_exact_prior",
  paste0("final_beta_tree_export_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
)
write.csv(summaries, combined_summary_path, row.names = FALSE)
cat("Wrote combined export summary:", combined_summary_path, "\n")
