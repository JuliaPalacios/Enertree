library(fmatrix)
library(phangorn)
# Hard-coded for num_tips=7, as a result do not change number of tips.
# Checked
compute_log_Z <- function(M_est, b_est) {
  sum_pmf <- 0
  for (mat in F.list7) {
    sum_pmf <- sum_pmf + 
      exp(-b_est * norm(as.numeric(mat)-M_est, type="F") ** 2)
  }
  return(log(sum_pmf))
}

compute_EF_gibbs <- function(M, b) {
  multiply_by_gibbs <- function(mat) {
    gibbs_pmf <- exp(-b * norm(as.numeric(mat)-M, 
                               type="F") ** 2)/exp(compute_log_Z(M,b))
    return(mat * gibbs_pmf)
  }
  dotted <- lapply(F.list7, FUN=multiply_by_gibbs)
  return(Reduce('+', dotted))
}


compute_Edstsq_gibbs <- function(M, b) {
  multiply_by_gibbs <- function(mat) {
    gibbs_pmf <- exp(-b * norm(as.numeric(mat)-M, 
                               type="F") ** 2)/exp(compute_log_Z(M,b))
    return(norm(as.numeric(mat)-M, type="F") ** 2 * gibbs_pmf)
  }
  dotted <- sapply(F.list7, FUN=multiply_by_gibbs)
  return(Reduce('+', dotted))
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

compute_grad_b <- function(data, 
                           M_est, 
                           b, 
                           num_tips, 
                           coal_times,
                           log_Z,
                           num_tip_label_iters) {
  expectation <- 0
  for (mat in F.list7) {
    l_lik <- avg_log_lik_tip_labels(data=data, 
                                    mat=mat,
                                    coal_times = coal_times,
                                    num_tip_label_iters = num_tip_label_iters)
    l_var <- -log(2 ** (num_tips - 1 -num_cherries(mat))) + 
      log(factorial(num_tips-1)) -log_Z - b * norm(M_est-mat, type="F") ** 2
    gibbs_pmf <- exp(-b * norm(mat-M_est, 
                               type="F") ** 2)/exp(log_Z)
    grad_b_log_q <- compute_Edstsq_gibbs(M_est, b) - norm(mat - M_est, type="F") ** 2
    expectation <- expectation + grad_b_log_q * gibbs_pmf * (l_lik - l_var)
    
  }
  return(expectation)
}

compute_grad_M <- function(data, 
                           M_est, 
                           b, 
                           num_tips, 
                           coal_times,
                           log_Z,
                           num_tip_label_iters) {
  expectation <- 0
  for (mat in F.list7) {
    gibbs_pmf <- exp(-b * norm(mat-M_est, 
                               type="F") ** 2)/exp(log_Z)
    l_lik <- avg_log_lik_tip_labels(data=data, 
                                    mat=mat,
                                    coal_times = coal_times,
                                    num_tip_label_iters = num_tip_label_iters)
    grad_M_log_q <- 2 * b * (mat - compute_EF_gibbs(M_est, b))
    l_var <- -log(2 ** (num_tips - 1 -num_cherries(mat))) + 
      log(factorial(num_tips-1)) -log_Z - b * norm(M_est-mat, type="F") ** 2
    expectation <- expectation + grad_M_log_q * gibbs_pmf * (l_lik - l_var)
  }
  return(expectation)
}

compute_elbo <- function(data, 
                         M_est, 
                         b, 
                         num_tips, 
                         coal_times,
                         log_Z,
                         num_tip_label_iters) {
  expectation <- 0
  for (mat in F.list7) {
    gibbs_pmf <- exp(-b * norm(mat-M_est, 
                               type="F") ** 2)/exp(log_Z)
    l_lik <- avg_log_lik_tip_labels(data=data, 
                                    mat=mat,
                                    coal_times = coal_times,
                                    num_tip_label_iters = num_tip_label_iters)
    l_var <- -log(2 ** (num_tips - 1 - num_cherries(mat)))+ 
      log(factorial(num_tips-1))-log_Z - b * norm(M_est-mat, type="F") ** 2
    expectation <- expectation + gibbs_pmf * (l_lik - l_var)
  }
  return(expectation)
}
