library(phylodyn)
library(ape)

set.seed(42)

sim <- phylodyn:::generate_true_M_and_data(
  num_tips = 10,
  rate = 0.01,
  seq_len = 10000,
  write_files = TRUE,
  name_fasta = "sequences.fasta",
  name_tree = "true_tree.newick"
)

true_tree <- sim$M_true_tree
sequences <- sim$sequences
