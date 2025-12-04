# CONTAINS WRAPPERS TO ESTIMATE EXPECTATIONS 
library(fmatrix)
library(phangorn)


# Computes nth zig-zag number.
zigzag <- function(n) {
  fact<-list()
  zig<-list()
  fact[[1]] = 1;
  for (i in (1:n)){
    fact[[i+1]] = fact[[i]] *i
  }
  zig[[1]] = 1;
  zig[[2]] = 1;
  for (i in (2:(n-1))){
    sum = 0
    for (k in (0:(i-1))){
      sum = sum+(fact[[i]]/(fact[[i- k]]*fact[[k+1]]))*zig[[k+1]]*zig[[i-k]]
    }
    zig[[i+1]] = sum / 2
  }
  return(zig[[length(zig)]])
}

# Computes the expected square distance w.r.t Gibbs.
#
# @param {array} samples: (num_tips-1, num_tips-1, num_samps) sized array
#                         of samples from the Gibbs or Uniform distribution.
# @param {matrix} M_est: estimate of mean parameter.
#
# @returns estimated expected distance squared to M.

expect_dst_sq <- function(samples, M_est) {
  dist_sq <- function(mat) {
    return(norm(as.numeric(mat)-M_est, type="F") ** 2)
  }
  dists <- sapply(samples, FUN=dist_sq, simplify=FALSE)
  return(Reduce('+', dists)/length(dists))
}

# Computes the expected F-matrix wrt Gibbs.
#
# @param {array} samples: (num_tips-1, num_tips-1, num_samps) sized array
#                         of samples from the Gibbs or Uniform distribution.
# @param {matrix} M_est: estimate of mean parameter.
#
# @returns estimated F-matrix expectation.

expect_F <- function(samples) {
  return(Reduce('+', samples)/length(samples))
}

# Computes the expected unnormalized Gibbs mass function. 
#
# @param {array} samples: (num_tips-1, num_tips-1, num_samps) sized array
#                         of samples from Uniform distribution.
# @param {matrix} M_est: estimate of mean parameter.
# @param {numeric} b_est: estimate of shape parameter
#
# @returns estimated expected unnormalized Gibbs mass function.

expect_gibbs_pmf <- function(samples, M_est, g_est) {
  gibbs_pmf <- function(mat) {
    return(exp(-exp(g_est) * norm(as.numeric(mat)-M_est, type="F") ** 2))
  }
  gibbs_probs <- apply(X=samples, MARGIN=c(3), FUN=gibbs_pmf)
  return(mean(gibbs_probs))
}

# Estimates log-normalizer of the Gibbs distribution.
#
# @param {array} samples: (num_tips-1, num_tips-1, num_samps) sized array
#                         of samples from Uniform distribution.
# @param {matrix} M_est: estimate of mean parameter.
# @param {numeric} b_est: estimate of shape parameter.
# @param {numeric} num_tips: number of tips

estimate_log_Z <- function(samples, M_est, b_est, num_tips) {
  return(log(zigzag(num_tips-1)) + log(expect_gibbs_pmf(samples, M_est, b_est)))
}

avg_log_lik_tip_labels <- function(data, 
                                   mat, 
                                   coal_times, 
                                   num_tip_label_iters) {
  
  lls <- c()
  for (i in 1:num_tip_label_iters) {
    tree <- mytree_from_F(mat, coal_times)
    curr_ll <- pml(tree=tree, data=data, model="JC69")$logLik
    lls <- c(lls, curr_ll)
  }
  return(mean(lls))
}

estimate_grad_g2 <- function(data, 
                              coal_times,
                              samples_gibbs, 
                              M_est, 
                              b_est, 
                              num_tips,
                              log_Z,
                              num_tip_label_iters) {
  expectation <- 0
  elbo <- 0
  for (mat in samples_gibbs) {
    tree <- mytree_from_F(mat, coal.times=coal_times)
    l_lik <- avg_log_lik_tip_labels(data=data, 
                                    mat=mat,
                                    coal_times = coal_times,
                                    num_tip_label_iters = num_tip_label_iters)
    l_var <- -log(2 ** (num_tips - 1 -num_cherries(mat))) + 
      log(factorial(num_tips-1)) -log_Z - b_est * norm(M_est-mat, type="F") ** 2
    grad_g_log_q <- b_est * (expect_dst_sq(samples_gibbs, M_est) - norm(mat - M_est, type="F") ** 2)
    expectation <- expectation + grad_b_log_q * (l_lik- l_var)
    elbo <- elbo + (l_lik- l_var)
    
  }
  return(list(grad=expectation,elbo=elbo))
}

estimate_grad_M2 <- function(data, 
                            coal_times,
                            samples_gibbs, 
                            M_est, 
                            b_est, 
                            num_tips,
                            log_Z,
                            num_tip_label_iters) {
  expectation <- 0
  elbo <- 0
  
  for (mat in samples_gibbs) {
    tree <- mytree_from_F(mat, coal.times=coal_times)
    l_lik <- avg_log_lik_tip_labels(data=data, 
                                    mat=mat,
                                    coal_times = coal_times,
                                    num_tip_label_iters = num_tip_label_iters)
    grad_M_log_q <- 2 * b_est * (mat - expect_F(samples_gibbs))
    l_var <- -log(2 ** (num_tips - 1 -num_cherries(mat))) + 
      log(factorial(num_tips-1)) -log_Z - b_est * norm(M_est-mat, type="F") ** 2
    expectation <- expectation + grad_M_log_q  * (l_lik - l_var)
    elbo <- elbo + (l_lik - l_var)
  }
  return(list(grad=expectation,elbo=elbo))
}


estimate_grad_M_g <- function(data, 
                             coal_times,
                             samples_gibbs, 
                             M_est, 
                             b_est, 
                             num_tips,
                             log_Z,
                             num_tip_label_iters) {
  expectation_M <- 0
  expectation_g <- 0
  elbo <- 0
  
  for (mat in samples_gibbs) {
    tree <- mytree_from_F(mat, coal.times=coal_times)
    l_lik <- avg_log_lik_tip_labels(data=data, 
                                    mat=mat,
                                    coal_times = coal_times,
                                    num_tip_label_iters = num_tip_label_iters)
    grad_M_log_q <- 2 * b_est * (mat - expect_F(samples_gibbs))
    l_var <- -log(2 ** (num_tips - 1 -num_cherries(mat))) + 
      log(factorial(num_tips-1)) -log_Z - b_est * norm(M_est-mat, type="F") ** 2
    expectation_M <- expectation_M + grad_M_log_q  * (l_lik - l_var)
    grad_g_log_q <- b_est * (expect_dst_sq(samples_gibbs, M_est) - norm(mat - M_est, type="F") ** 2)
    expectation_g <- expectation_g + grad_g_log_q * (l_lik- l_var)
    elbo <- elbo + (l_lik - l_var)
  }
  return(list(grad_g=expectation_g,grad_M=expectation_M,elbo=elbo))
}


estimate_elbo <- function(data, 
                          coal_times,
                          samples_gibbs,
                          M_est, 
                          b_est, 
                          num_tips,
                          log_Z,
                          num_tip_label_iters) {
  
  expectation <- 0
  for (mat in samples_gibbs) {
    tree <- mytree_from_F(mat, coal.times=coal_times)
    l_lik <- avg_log_lik_tip_labels(data=data, 
                                    mat=mat,
                                    coal_times = coal_times,
                                    num_tip_label_iters = num_tip_label_iters)
    l_var <- -log(2 ** (num_tips - 1 - num_cherries(mat)))+ 
      log(factorial(num_tips-1))-log_Z - b_est * norm(M_est-mat, type="F") ** 2
    expectation <- expectation +  (l_lik - l_var)
  }
  return(expectation)
}
