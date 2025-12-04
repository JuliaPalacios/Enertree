#setwd("~/My Drive/Statistics/Phylo-VB")
library(ape)
library(phylodyn)
library(phyclust)
library(phangorn)
library(phylotools)
library(devtools)
#library(gurobi)
library(fmatrix)
library(CVXR)

source("./estimate_topology_utils.R")
source("./expectations.R")
source("./estimate_coal_utils.R")

gen_caterpillar <- function(n){
  
  #"Return caterpillar tree"
  F_mat <- matrix(rep(seq(1,n-1),n-1),nrow=n-1,byrow = T)
  F_mat[upper.tri(F_mat)] <- 0
  diag(F_mat) <- seq(2,n)
  return(F_mat)
}

num_tips <- 10
step_size <- 0.1

num_samps <- 200
num_grad_desc_steps <- 50
num_tip_label_iters <- 30
num_gibbs_iters <- 50
n <- num_tips - 1


###We run directly "estimate topology"
###############################
"Generate data and initialize"
###############################
init_results <- generate_true_M_and_data(num_tips)

# set branch lengths for the initial guess of M.
inter_coal_times <- coalescent.intervals(init_results$M_true_tree)$interval.length
inter_coal_times[which(inter_coal_times<=0.001)] <- 0.01
coal_times <- cumsum(inter_coal_times)
#Initialization
init_Fmat <- gen_caterpillar(num_tips)
M_est_tree <- mytree_from_F(init_Fmat, coal_times)
#init_Fmat <-  sampleF(init_results$M_true_fmat, 1,1000 , diam=1)$chainF[100]
#M_est_tree <- mytree_from_F(init_Fmat[[1]], coal_times)
#M<-init_results$M_true_fmat
#M_full_tree <- init_results$M_true_tree
#F <- Variable(rows=n, cols=n)
#constraints <- generate_constraints(F, num_tips)
M_est <- round(gen_Fmat(M_est_tree, tol=8),0)
output_prefix <- "sample"

g_est <- log(0.01)
distance_to_truth <- c()
elbos <- c()
elbos2 <- c()
burnin <- round(num_samps*0.1)
########################################
"Essentially, it runs estimate topology"
########################################
  
  ##add a termination criteria for difference in M, whatever happens firs
  M_grads <- list(matrix(step_size, nrow=nrow(M_est), ncol=nrow(M_est)))
  g_grads <- list(step_size)
  joint <- TRUE
  start_time <- Sys.time()
  for (i in 1:num_grad_desc_steps) {
  print(paste0("Gradient Descent Iteration ", i))
  print(M_est)
  print(paste0("Estimated log beta: ", g_est))
  samples <- generate_sample_markov_chain(M=M_est,
                                          b=exp(g_est),
                                          num_samps=num_samps,
                                          burnin=burnin)
  samples_uniform <- generate_sample_uniform(num_tips=num_tips,
                                             num_samps=num_samps)
  log_Z <- estimate_log_Z(samples=samples_uniform,
                          M_est=M_est,
                          b_est=exp(g_est),
                          num_tips=num_tips)
  if (joint==FALSE){
    g_up <- estimate_grad_g2(data=init_results$sequences, 
                              coal_times=coal_times,
                              samples_gibbs=samples, 
                              M_est=M_est, 
                              b_est=exp(g_est), 
                              num_tips=num_tips,
                              log_Z=log_Z,
                              num_tip_label_iters=num_tip_label_iters)
    g_grad <- g_up$grad
    g_grads[[length(g_grads)+1]] <- g_grad ** 2
    g_est <- g_est + step_size * 1/sqrt(Reduce('+', g_grads)) * g_grad
    samples <- generate_sample_markov_chain(M=M_est,
                                            b=exp(g_est),
                                            num_samps=num_samps,
                                            burnin=burnin)
    samples_uniform <- generate_sample_uniform(num_tips=num_tips,
                                               num_samps=num_samps)
    log_Z <- estimate_log_Z(samples=samples_uniform,
                            M_est=M_est,
                            b_est=exp(g_est),
                            num_tips=num_tips)
    M_up <- estimate_grad_M2(data=init_results$sequences, 
                              coal_times=coal_times,
                              samples_gibbs=samples, 
                              M_est=M_est, 
                              b_est=exp(g_est), 
                              num_tips=num_tips,
                              log_Z=log_Z,
                              num_tip_label_iters=num_tip_label_iters)
    M_grad <- M_up$grad
    M_grads[[length(M_grads)+1]] <- M_grad ** 2
    M_est <- M_est + step_size * 1/(sqrt(Reduce('+', M_grads))) * M_grad
    elbos2 <- c(elbos2,M_up$elbo)
  } else {
    J_up <- estimate_grad_M_g(data=init_results$sequences, 
                             coal_times=coal_times,
                             samples_gibbs=samples, 
                             M_est=M_est, 
                             b_est=exp(g_est), 
                             num_tips=num_tips,
                             log_Z=log_Z,
                             num_tip_label_iters=num_tip_label_iters)
    M_grad <- J_up$grad_M
    M_grads[[length(M_grads)+1]] <- M_grad ** 2
    M_est <- M_est + step_size * 1/(sqrt(Reduce('+', M_grads))) * M_grad
    g_grad <- J_up$grad_g
    g_grads[[length(g_grads)+1]] <- g_grad ** 2
    g_est <- g_est + step_size * 1/sqrt(Reduce('+', g_grads)) * g_grad
    elbos2 <- c(elbos2,J_up$elbo)
  }
  #Potentially, this can be removed to be more efficient
  # samples <- generate_sample_markov_chain(M=M_est,
  #                                         b=b_est,
  #                                         num_samps=num_samps,
  #                                         burnin=burnin)
  # samples_uniform <- generate_sample_uniform(num_tips=num_tips,
  #                                            num_samps=num_samps)
  # log_Z <- estimate_log_Z(samples=samples_uniform,
  #                         M_est=M_est,
  #                         b_est=b_est,
  #                         num_tips=num_tips)
  # elbo_est <- estimate_elbo(data=init_results$sequences, 
  #                           coal_times=coal_times,
  #                           samples_gibbs=samples,
  #                           M_est=M_est, 
  #                           b_est=b_est, 
  #                           num_tips=num_tips,
  #                           log_Z=log_Z,
  #                           num_tip_label_iters=num_tip_label_iters)
  #print(elbo_est)
  #elbos <- c(elbos, elbo_est)
  }
  end_time <- Sys.time()
  end_time - start_time
  plot(elbos2)
  
  M_estimated_tree <- mytree_from_F(nearby_Fmat_sampling(M_est), coal_times)
plot(M_estimated_tree)
plot(init_results$M_true_tree)


# @param {numeric} step_size: step size for gradient descent
# @param {integer} num_samps: number of samples to estimate expectations.
# @param {integer} num_gibbs_iters: number of Gibbs iterations to run.
#run_gibbs_instance(num_tips, step_size, num_samps, num_gibbs_iters) 




run_gibbs_instance <- function(num_tips, step_size, num_samps, num_gibbs_iters) {
  init_results <- generate_true_M_and_data(num_tips)
  # We will fix the initial guess of M_est to be the balanced tree
  # to evaluate the ELBO better.

  # set branch lengths for the initial guess of M.
  inter_coal_times <- coalescent.intervals(init_results$M_true_tree)$interval.length
  inter_coal_times[which(inter_coal_times<=0.001)] <- 0.01
  coal_times <- cumsum(inter_coal_times)
  #Initialization
  init_Fmat <- gen_caterpillar(num_tips)
  M_est_tree <- mytree_from_F(init_Fmat, coal_times)
  M<-init_results$M_true_fmat
  M_full_tree <- init_results$M_true_tree
  #F <- Variable(rows=n, cols=n)
  #constraints <- generate_constraints(F, num_tips)
  M_est <- round(gen_Fmat(M_est_tree, tol=8),0)
  output_prefix <- "sample"
  b_est <- 1
  distance_to_truth <- c()
  elbos <- c()
  for (i in 1:num_gibbs_iters) {
    estimated_params <- estimate_topology(M=M_est, 
                                          b=b_est, 
                                          num_tips=num_tips, 
                                          num_samps=num_samps,
                                          burnin = round(num_samps*0.1),
                                          data=init_results$sequences,
                                          coal_times=coal_times,
                                          step_size=step_size,
                                          num_grad_desc_steps = num_grad_desc_steps,
                                          num_tip_label_iters = num_tip_label_iters)
    M_est <- estimated_params$M_est
    b_est <- estimated_params$b_est
    elbos <- c(elbos,estimated_params$elbos)
    M_est_tree <- mytree_from_F(nearby_Fmat_sampling(M_est), coal_times)
    write.tree(M_est_tree, 
               file=paste0(output_prefix, '.newick'))
    
    
  }
  return(list(dist_to_truth=norm(M_est - M, type="F")))
}

avg_dist_to_true_topology <- function(num_iters, 
                                      num_tips, 
                                      step_size, 
                                      num_samps, 
                                      num_gibbs_iters) {
  dists <- c()
  for (i in 1:num_iters) {
    dists <- c(dists, run_gibbs_instance(num_tips=num_tips,
                                         step_size=step_size,
                                         num_samps=num_samps,
                                         num_gibbs_iters=num_gibbs_iters)$dist_to_truth)
  }
  print(dists)
  hist(dists)
}

#avg_dist_to_true_topology(num_iters=1, num_tips, step_size, num_samps, num_gibbs_iters)









