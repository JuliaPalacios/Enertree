# Export sequences and fresh trees sampled at the fitted final A-matrix weights.
# Run with:
# Rscript "Brandon/gibbs_prior/A matrix analysis/export_final_A_samples.R"
# An optional first argument selects an experiment ID, result.rds file, or its containing folder.
# An optional second argument selects the output directory.

# Edit this when running the whole file with RStudio's Run button.
# Example: experiment_id <- "20261025_221200"
# A command-line experiment argument takes precedence over this setting.
experiment_id <- "20261024_184601"

a_export_source_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
a_export_console_run <- FALSE
if (is.null(a_export_source_file)) {
  a_export_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(a_export_file_arg)) {
    a_export_source_file <- sub("^--file=", "", a_export_file_arg[1L])
  }
}
# Rscript encodes spaces as ~+~ in the --file argument on this platform.
if (!is.null(a_export_source_file) && !file.exists(a_export_source_file)) {
  a_export_source_file <- gsub("~\\+~", " ", a_export_source_file)
}
if (is.null(a_export_source_file) || !file.exists(a_export_source_file)) {
  # RStudio's Run button sends lines to the Console, where neither `ofile`
  # nor `--file` is available. Resolve the script from the working directory.
  a_export_candidates <- file.path(
    getwd(), c("A matrix analysis/export_final_A_samples.R",
               "export_final_A_samples.R",
               "Brandon/gibbs_prior/A matrix analysis/export_final_A_samples.R")
  )
  a_export_candidates <- a_export_candidates[file.exists(a_export_candidates)]
  if (!length(a_export_candidates)) {
    stop("Cannot find export_final_A_samples.R from working directory: ", getwd(),
         ". Set the working directory to gibbs_prior or run with Rscript.")
  }
  a_export_source_file <- a_export_candidates[1L]
  a_export_console_run <- TRUE
}
a_export_dir <- dirname(normalizePath(a_export_source_file, mustWork = TRUE))
a_gibbs_dir <- dirname(a_export_dir)
find_A_result <- function(path = NULL) {
  if (!is.null(path)) {
    if (length(path) != 1L || is.na(path) || !nzchar(path)) {
      stop("Pass one nonempty experiment ID, result.rds path, or experiment folder.")
    }
    # A bare ID names a folder under the saved VB_master_A runs.
    if (!file.exists(path) && !dir.exists(path) && basename(path) == path) {
      candidate <- file.path(a_gibbs_dir, "logs", "VB_master_A", path)
      if (dir.exists(candidate)) path <- candidate
    }
    if (dir.exists(path)) path <- file.path(path, "result.rds")
    if (!file.exists(path)) {
      stop("Cannot find result.rds for: ", path,
           ". Pass an experiment ID from logs/VB_master_A, its folder, or the result.rds path.")
    }
    return(normalizePath(path, mustWork = TRUE))
  }

  log_root <- file.path(a_gibbs_dir, "logs", "VB_master_A")
  candidates <- if (dir.exists(log_root)) {
    list.files(log_root, pattern = "^result\\.rds$", recursive = TRUE,
               full.names = TRUE)
  } else character()
  if (!length(candidates)) {
    stop("No completed VB_master_A result.rds found in ", log_root,
         ". Wait for the run to finish or pass a result path explicitly.")
  }
  candidates[which.max(file.info(candidates)$mtime)]
}

export_final_A_samples <- function(result_path = NULL, output_dir = NULL,
                                   n_samples = 1000L, seed = 20261024L) {
  if (length(n_samples) != 1L || !is.finite(n_samples) ||
      n_samples != as.integer(n_samples) || n_samples < 1L) {
    stop("n_samples must be a positive integer.")
  }
  if (length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed)) {
    stop("seed must be one finite integer.")
  }
  result_path <- find_A_result(result_path)
  fit <- readRDS(result_path)
  needed <- c("theta", "chain_states", "coal_times", "sequences", "config")
  missing <- setdiff(needed, names(fit))
  if (length(missing)) stop("The saved fit is missing: ", paste(missing, collapse = ", "))

  suppressPackageStartupMessages(library(phangorn))
  if (!inherits(fit$sequences, "phyDat")) stop("fit$sequences must be a phyDat alignment.")
  labels <- names(fit$sequences)
  theta <- fit$theta
  starts <- fit$chain_states
  if (!is.matrix(theta) || any(!is.finite(theta)) ||
      ncol(theta) != length(labels) || nrow(theta) != length(labels) - 1L ||
      anyNA(labels) || anyDuplicated(labels) ||
      !identical(colnames(theta), labels)) {
    stop("The fitted theta and sequence taxon labels do not align.")
  }
  if (!length(starts) || any(!vapply(starts, function(A) {
    is.matrix(A) && identical(dim(A), dim(theta)) &&
      identical(colnames(A), labels) && phylodyn:::a_matrix_valid(A)
  }, logical(1L)))) {
    stop("The saved chain states must be valid A matrices with matching taxa.")
  }
  if (length(fit$coal_times) != length(labels) - 1L ||
      any(!is.finite(fit$coal_times)) || any(fit$coal_times <= 0) ||
      is.unsorted(fit$coal_times, strictly = TRUE)) {
    stop("Invalid fixed coalescent times in the saved fit.")
  }
  burnin <- as.integer(fit$config$burnin)
  thin <- as.integer(fit$config$take_every)
  if (length(burnin) != 1L || is.na(burnin) || burnin < 0L ||
      length(thin) != 1L || is.na(thin) || thin < 1L) {
    stop("Invalid burn-in or thinning in the saved fit configuration.")
  }

  if (is.null(output_dir)) {
    output_dir <- file.path(a_export_dir, "outputs", basename(dirname(result_path)))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  # Share draws as evenly as possible among the saved final chain states.
  draws_per_chain <- rep(n_samples %/% length(starts), length(starts))
  remainder <- n_samples %% length(starts)
  if (remainder) draws_per_chain[seq_len(remainder)] <-
    draws_per_chain[seq_len(remainder)] + 1L

  set.seed(seed)
  chains <- lapply(seq_along(starts), function(k) {
    if (draws_per_chain[k] == 0L) return(NULL)
    phylodyn:::a_rnni_sample(
      theta = theta, initial = starts[[k]],
      n_samples = draws_per_chain[k], burnin = burnin, thin = thin
    )
  })
  A_samples <- unlist(lapply(chains, function(chain) {
    if (is.null(chain)) list() else chain$samples
  }), recursive = FALSE)
  stopifnot(length(A_samples) == n_samples)
  chain_id <- rep(seq_along(chains), draws_per_chain)
  draw_in_chain <- unlist(lapply(draws_per_chain, seq_len), use.names = FALSE)

  # Reuse the conversion in the fitting script so ranks, times and labels agree.
  helpers <- new.env(parent = globalenv())
  sys.source(file.path(a_gibbs_dir, "VB_master_A.R"), envir = helpers)
  trees <- lapply(A_samples, helpers$vb_A_to_tree,
                  coal_times = fit$coal_times, leaf_order = labels)
  names(trees) <- sprintf("sample_%04d", seq_len(n_samples))
  if (!all(vapply(trees, function(tree) {
    setequal(tree$tip.label, labels) &&
      all(is.finite(tree$edge.length)) && all(tree$edge.length > 0)
  }, logical(1L)))) {
    stop("At least one sampled tree has invalid labels or branch lengths.")
  }
  class(trees) <- "multiPhylo"

  fasta_path <- file.path(output_dir, "sequences.fasta")
  newick_path <- file.path(output_dir, "sampled_trees.newick")
  nexus_path <- file.path(output_dir, "sampled_trees.trees")
  sample_path <- file.path(output_dir, "sampled_A_matrices.rds")
  manifest_path <- file.path(output_dir, "sample_manifest.csv")
  dna <- ape::as.DNAbin(fit$sequences)
  ape::write.FASTA(dna, file = fasta_path)
  ape::write.tree(trees, file = newick_path)
  ape::write.nexus(trees, file = nexus_path, translate = FALSE)
  if (!is.null(fit$true_tree)) {
    if (!inherits(fit$true_tree, "phylo") ||
        !setequal(fit$true_tree$tip.label, labels)) {
      stop("The saved true tree does not match the exported sequence taxa.")
    }
    ape::write.tree(fit$true_tree,
                    file = file.path(output_dir, "true_tree.newick"))
  }
  saveRDS(A_samples, sample_path)
  manifest <- data.frame(
    sample_id = names(trees), chain_id = chain_id,
    draw_in_chain = draw_in_chain,
    log_unnormalized_q = vapply(A_samples, function(A) sum(theta * A), numeric(1L))
  )
  utils::write.csv(manifest, manifest_path, row.names = FALSE)

  # Check that the shared exports can be read back with all taxa and sites.
  fasta_check <- ape::read.FASTA(fasta_path)
  tree_check <- ape::read.tree(newick_path)
  nexus_check <- ape::read.nexus(nexus_path)
  if (length(fasta_check) != length(labels) ||
      !setequal(names(fasta_check), labels) ||
      !identical(unname(as.matrix(fasta_check)), unname(dna)) ||
      length(tree_check) != n_samples || length(nexus_check) != n_samples ||
      !all(vapply(tree_check, function(tree) setequal(tree$tip.label, labels), logical(1L))) ||
      !all(vapply(nexus_check, function(tree) setequal(tree$tip.label, labels), logical(1L)))) {
    stop("An exported FASTA or tree file failed its read-back check.")
  }

  summary <- list(source_result = result_path, output_dir = output_dir,
                  n_samples = n_samples, seed = seed,
                  draws_per_chain = draws_per_chain,
                  burnin_per_chain = burnin, thin = thin,
                  acceptance_rate = vapply(chains, function(chain) {
                    if (is.null(chain)) NA_real_ else chain$acceptance_rate
                  }, numeric(1L)),
                  unique_A_matrices = length(unique(vapply(
                    A_samples, phylodyn:::a_matrix_key, character(1L)
                  ))))
  saveRDS(summary, file.path(output_dir, "export_summary.rds"))
  cat("Source fit:", result_path, "\n")
  cat("Exported", n_samples, "tree draws to:", output_dir, "\n")
  cat("Distinct sampled A matrices:", summary$unique_A_matrices, "\n")
  invisible(summary)
}

if (sys.nframe() == 0L || a_export_console_run) {
  args <- if (a_export_console_run) character() else commandArgs(trailingOnly = TRUE)
  export_final_A_samples(
    result_path = if (length(args) >= 1L) args[1L] else experiment_id,
    output_dir = if (length(args) >= 2L) args[2L] else NULL
  )
}

