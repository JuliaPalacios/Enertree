# CONTAINS ALL THE HELPER CODE FOR GENERATING SAMPLES USING
# THE ARGMAX TRICK, GENERATING TRUE DATA, AND ESTIMATING
# TOPOLOGY VIA GRADIENT ASCENT.

library("ape")
library(phylodyn)
library("phyclust",quiet=TRUE)
library("phylotools")
library("phangorn")
library(phylotools)
library(devtools)
library(fmatrix)
library(CVXR)

source("./expectations.R")
source("./compute_exactly_utils.R")
source("./estimate_coal_utils.R")
source("./all_functions.R")
generate_true_M_and_data<- function(num_tips,traj="exp_traj") {
  samp_times <- 0
  M_true_tree <- generate_newick(
    coalsim(
      samp_times = samp_times, 
      n_sampled = num_tips, 
      exp_traj ,
      method="tt",
      val_upper=11
    )
  )
  M_true_tree$newick$tip.label <- sapply(M_true_tree$newick$tip.label, 
                                         function(x) sub("_0", "", x))
  M_true_tree$labels <- sapply(M_true_tree$labels, 
                               function(x) sub("_0", "", x))
  M_true_newick_string <-write.tree(M_true_tree$newick)
  # generate data using our true M 
  seqgen(opts="-mHKY -t0.5 -f0.25,0.25,0.25,0.25 -l100",
         newick.tree=M_true_newick_string,
         temp.file="temp.fasta")
  data2<-read.phylip("temp.fasta")
  
  data<-dat2fasta(data2,outfile="out.fasta")
  mydna<-as.phyDat(read.FASTA("out.fasta"))
  return(list(M_true_tree = M_true_tree$newick,
              M_true_fmat = gen_Fmat(M_true_tree$newick,tol=8),
              sequences = mydna))
}

unnormalized_gibbs_pmf <- function(M, b, mat) {
  return(exp(-b * norm(as.numeric(mat)-M, type="F") ** 2))
}


generate_sample_markov_chain <- function(M, b,num_samps,burnin) {
  # space out the indices by some number of iterations.
  indices <- c()
  for (i in 1:num_samps) {
    indices <- c(indices, 20 * i)
  }
  samples <- sampleF(M, b, burnin+20*num_samps, diam=1)$chainF[indices]
  return(samples)
}

#
# Generate an approximate sample of F matrices from Uniform distribution.
#
# @param {integer} num_tips: number of tips
# @param {integer} num_samps: number of F-matrices to sample from distribution
#
# @returns {array}: (num_tips-1, num_tips-1, num_samps) dimensional array.

generate_sample_uniform <- function(num_tips, num_samps) {
  samples <- c()
  for (i in 1:num_samps) {
    samples <- c(samples, rFmat(num_tips))
  }
  return(array(samples, c(num_tips-1, num_tips-1, num_samps)))
}

estimate_topology <- function(M, 
                              b, 
                              num_tips, 
                              num_samps,
                              burnin,
                              data,
                              coal_times,
                              step_size,
                              num_grad_desc_steps,
                              num_tip_label_iters) {
  ##add a termination criteria for difference in M, whatever happens first
  M_est <- M
  g_est <- exp(1)
  elbos <- c()
  M_grads <- list(matrix(step_size, nrow=nrow(M_est), ncol=nrow(M_est)))
  g_grads <- list(step_size)
  #for (i in 1:num_grad_desc_steps) {
    print(paste0("Gradient Descent Iteration ", i))
    print(M_est)
    print(g_est)
    samples <- generate_sample_markov_chain(M=M_est,
                                            b=exp(g_est),
                                            num_samps=num_samps,
                                            burnin=burnin)
    samples_uniform <- generate_sample_uniform(num_tips=num_tips,
                                               num_samps=num_samps)
    log_Z <- estimate_log_Z(samples=samples_uniform,
                            M_est=M_est,
                            g_est=exp(b_est),
                            num_tips=num_tips)
    g_grad <- estimate_grad_g(data=data, 
                              coal_times=coal_times,
                              samples_gibbs=samples, 
                              M_est=M_est, 
                              g_est=exp(g_est), 
                              num_tips=num_tips,
                              log_Z=log_Z,
                              num_tip_label_iters=num_tip_label_iters)
    g_grads[[length(g_grads)+1]] <- g_grad ** 2
    g_est <- g_est + step_size * 1/sqrt(Reduce('+', g_grads)) * g_grad
    samples <- generate_sample_markov_chain(M=M_est,
                                            g=g_est,
                                            num_samps=num_samps,
                                            burnin=burnin)
    samples_uniform <- generate_sample_uniform(num_tips=num_tips,
                                               num_samps=num_samps)
    log_Z <- estimate_log_Z(samples=samples_uniform,
                            M_est=M_est,
                            g_est=g_est,
                            num_tips=num_tips)
    M_grad <- estimate_grad_M(data=data, 
                              coal_times=coal_times,
                              samples_gibbs=samples, 
                              M_est=M_est, 
                              g_est=b_est, 
                              num_tips=num_tips,
                              log_Z=log_Z,
                              num_tip_label_iters=num_tip_label_iters)
    M_grads[[length(M_grads)+1]] <- M_grad ** 2
    M_est <- M_est + step_size * 1/(sqrt(Reduce('+', M_grads))) * M_grad
    
    #Potentially, this can be removed to be more efficient
    samples <- generate_sample_markov_chain(M=M_est,
                                            g=g_est,
                                            num_samps=num_samps,
                                            burnin=burnin)
    samples_uniform <- generate_sample_uniform(num_tips=num_tips,
                                               num_samps=num_samps)
    log_Z <- estimate_log_Z(samples=samples_uniform,
                            M_est=M_est,
                            g_est=g_est,
                            num_tips=num_tips)
    elbo_est <- estimate_elbo(data=data, 
                              coal_times=coal_times,
                              samples_gibbs=samples,
                              M_est=M_est, 
                              g_est=g_est, 
                              num_tips=num_tips,
                              log_Z=log_Z,
                              num_tip_label_iters=num_tip_label_iters)
    print(elbo_est)
    elbos <- c(elbos, elbo_est)
  #}
  
  #plot(elbos)
  
  return(list(M_est=M_est, b_est=exp(g_est),elbos=elbos))
}




estimate_topology2 <-  function(M_est, 
                                g_est, 
                                num_tips, 
                                num_samps,
                                burnin,
                                data,
                                coal_times,
                                step_size,
                                num_grad_desc_steps,
                                num_tip_label_iters,
                                joint=TRUE) {

"Mostly taken from main2.R, I see g_est modification. To double check"
   
  elbos2 <- c()
  ##add a termination criteria for difference in M, whatever happens firs
  M_grads <- list(matrix(step_size, nrow=nrow(M_est), ncol=nrow(M_est)))
  g_grads <- list(step_size)
  joint <- TRUE
  #start_time <- Sys.time()
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
  }
return(list(M_est=M_est, g_est=g_est,elbos=elbos2))
}
