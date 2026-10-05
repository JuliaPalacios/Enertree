suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(phylodyn)
  library(phyclust)
  library(phylotools)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
setwd(script_dir)

set.seed(42)

num_tips <- 10
mutation_rate <- 0.01
sequence_length <- 1000

true_tree_file <- "rf_true_tree_10tip.newick"
sequence_file <- "rf_sequences_10tip.fasta"

generate_true_tree_and_sequences <- function(num_tips,
                                             traj = phylodyn:::exp_traj,
                                             rate = 1,
                                             seq_len = 1000,
                                             name_fasta = "sequences.fasta",
                                             name_tree = "true_tree.newick") {
  samp_times <- 0
  true_tree <- phylodyn:::generate_newick(
    phylodyn:::coalsim(
      samp_times = samp_times,
      n_sampled = num_tips,
      traj = traj,
      method = "tt",
      val_upper = 11
    )
  )

  true_tree$newick$tip.label <- vapply(
    true_tree$newick$tip.label,
    function(x) sub("_0", "", x),
    character(1)
  )
  true_tree$labels <- vapply(
    true_tree$labels,
    function(x) sub("_0", "", x),
    character(1)
  )

  true_tree_newick <- write.tree(true_tree$newick)
  seqgen_opts <- paste0(
    "-mHKY -t2.0 -f0.25,0.25,0.25,0.25 -l",
    seq_len,
    " -s",
    rate
  )

  seqgen(
    opts = seqgen_opts,
    newick.tree = true_tree_newick,
    temp.file = "seqgen_temp.phylip"
  )

  seqgen_output <- read.phylip("seqgen_temp.phylip")
  dat2fasta(seqgen_output, outfile = name_fasta)
  sequences <- as.phyDat(read.FASTA(name_fasta))
  write.tree(true_tree$newick, file = name_tree)

  if (file.exists("seqgen_temp.phylip")) {
    unlink("seqgen_temp.phylip")
  }

  list(
    true_tree = true_tree$newick,
    sequences = sequences,
    settings = list(
      num_tips = num_tips,
      mutation_rate = rate,
      sequence_length = seq_len,
      true_tree_file = name_tree,
      sequence_file = name_fasta,
      seed = 20261005
    )
  )
}

sim <- generate_true_tree_and_sequences(
  num_tips = num_tips,
  rate = mutation_rate,
  seq_len = sequence_length,
  name_fasta = sequence_file,
  name_tree = true_tree_file
)

writeLines(
  c(
    "RF Boltzmann data set",
    paste("num_tips:", sim$settings$num_tips),
    paste("mutation_rate:", sim$settings$mutation_rate),
    paste("sequence_length:", sim$settings$sequence_length),
    paste("seed:", sim$settings$seed),
    paste("true_tree_file:", sim$settings$true_tree_file),
    paste("sequence_file:", sim$settings$sequence_file)
  ),
  con = "README_data.txt"
)

cat("Wrote", true_tree_file, "and", sequence_file, "in", script_dir, "\n")

