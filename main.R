library(ape)
library(phylodyn)
library(phyclust)
library(phangorn)
library(phylotools)
library(devtools)
library(fmatrix)
library(CVXR)

source("./estimate_topology_utils.R")
source("./expectations.R")
source("./estimate_coal_utils.R")

num_tips <- 10
step_size <- 0.1

num_samps <- 100
num_grad_desc_steps <- 200
num_tip_label_iters <- 100
num_gibbs_iters <- 2
n <- num_tips - 1

# Creates new data, a new M, etc. then estimates topology
# in order to compute diagnostics such as average distance to the true topology.
#
# @param {numeric} step_size: step size for gradient descent
# @param {integer} num_samps: number of samples to estimate expectations.
# @param {integer} num_gibbs_iters: number of Gibbs iterations to run.
run_gibbs_instance <- function(num_tips, step_size, num_samps, num_gibbs_iters) {
  init_results <- generate_true_M_and_data(num_tips)
  # We will fix the initial guess of M_est to be the balanced tree
  # to evaluate the ELBO better.
  balanced <- matrix(c(2,0,0,0,0,0,0,0,0,
                       1,3,0,0,0,0,0,0,0,
                       0,2,4,0,0,0,0,0,0,
                       0,1,3,5,0,0,0,0,0,
                       0,0,2,4,6,0,0,0,0,
                       0,0,1,3,5,7,0,0,0,
                       0,0,0,2,4,6,8,0,0,
                       0,0,0,1,3,5,7,9,0,
                       0,0,0,0,2,4,6,8,10), nrow=9, ncol=9, byrow=TRUE)
  # set branch lengths for the initial guess of M.
  inter_coal_times <- coalescent.intervals(init_results$M_true_tree)$interval.length
  inter_coal_times[which(inter_coal_times<=0.001)] <- 0.01
  coal_times <- cumsum(inter_coal_times)
  M_est_tree <- mytree_from_F(balanced, coal_times)
  M<-init_results$M_true_fmat
  M_full_tree <- init_results$M_true_tree
  output_prefix <- "sample"
  b_est <- 1
  distance_to_truth <- c()
  coal_times_est <- coalsim(samp_times = 0, 
                            n_sampled = num_tips, 
                            traj = exp_traj,
                            method="tt",
                            val_upper=11)$coal_times
  for (i in 1:num_gibbs_iters) {
    # Update topology
    M_est <- round(gen_Fmat(M_est_tree, tol=8),0)
    estimated_params <- estimate_topology(M=M_est, 
                                          b=b_est, 
                                          num_tips=num_tips, 
                                          num_samps=num_samps,
                                          burnin = round(num_samps*0.1),
                                          data=init_results$sequences,
                                          coal_times=coal_times_est,
                                          step_size=step_size,
                                          num_grad_desc_steps = num_grad_desc_steps,
                                          num_tip_label_iters = num_tip_label_iters)
    M_est <- estimated_params$M_est
    b_est <- 1
    M_est_tree <- mytree_from_F(nearby_Fmat(M_est), coal_times_est)
    write.tree(M_est_tree, 
               file=paste0(output_prefix, '.newick'))
    # Update coalescent times
    infer_coal_times(data_filename="out.fasta",
                     tree_filename=paste0(output_prefix, '.newick'),
                     output_prefix=output_prefix)
    coal_times_est <- read_coal_times(filename=output_prefix,
                                      num_tips=num_tips)
    M_est_tree <- update_time(M_est_tree, coal_times_est)
    
  }
  return(list(M_est = M_est,
              F_est=nearby_Fmat(M_est),
              F_true=M,
              dist_to_truth=norm(M_est - M, type="F"),
              coal_times_est=coal_times_est,
              true_coal_times=coal_times))
}

avg_dist_to_true_topology <- function(num_tips, 
                                      step_size, 
                                      num_samps, 
                                      num_gibbs_iters) {
  dists <- c()
  F_est <- NULL
  M_est <- NULL
  F_true <- NULL
  est_coal_times <- NULL
  true_coal_times <- NULL
  for (i in 1:num_gibbs_iters) {
    result <- run_gibbs_instance(num_tips=num_tips,
                                 step_size=step_size,
                                 num_samps=num_samps,
                                 num_gibbs_iters=num_gibbs_iters)
    dists <- c(dists, result$dist_to_truth)
    F_true <- result$F_true
    F_est <- result$F_est
    M_est <- result$M_est
    est_coal_times <- result$coal_times_est
    true_coal_times <- result$true_coal_times
  }
  print(dists)
  hist(dists)
  return(list(F_est=F_est,
              M_est=M_est,
              F_true=F_true,
              true_coal_times=true_coal_times,
              est_coal_times=est_coal_times))
}

tree_results <- avg_dist_to_true_topology(num_tips, step_size, num_samps, num_gibbs_iters)

