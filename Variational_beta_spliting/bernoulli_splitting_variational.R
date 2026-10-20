seed                <- 42L
num_tips            <- 10
rate                <- 0.01
sequence_length     <- 10000
step_size           <- 0.01

fasta_file <- "sequences.fasta"
tree_file <- "true_tree.newick"


phylodyn:::generate_true_M_and_data(
  seed,
  num_tips,
  rate = rate,
  seq_len = sequence_length,
  write_files = save_input_data,
  name_fasta = fasta_file,
  name_tree = tree_file
)