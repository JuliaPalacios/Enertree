# Small smoke test for gradient_ascent_M_beta to diagnose "nothing updated"
suppressPackageStartupMessages({
  library(ape); library(phylodyn); library(phangorn); library(phyclust)
  library(phylotools); library(fmatrix); library(coda); library(ggplot2)
})

source("gibbs_prior_joint_standardized_helpers.R")
source("utils.R")
source("gradient_ascent_M_beta.R", local = FALSE, echo = FALSE,
       # Don't re-run the example block at the bottom of that file:
       chdir = TRUE)
# Re-define just our functions to be safe (the source above already re-runs them).

n_tips <- 25
set.seed(123)
sequences <- as.phyDat(read.FASTA("sequences.fasta"))
M_true_tree <- read.tree("true_tree.newick")
inter_coal_times <- coalescent.intervals(M_true_tree)$interval.length
inter_coal_times[which(inter_coal_times <= 1e-6)] <- 1e-7
true_coalescent_times <- cumsum(inter_coal_times)
M_true_fmat <- gen_Fmat(M_true_tree, tol = 8)
dm <- dist.hamming(sequences)
upgma_tree <- upgma(dm)
upgma_tree <- update_time(upgma_tree, true_coalescent_times)
M_upgma_fmat <- gen_Fmat(upgma_tree, tol = 8)

cat("\n=== Initial state ===\n")
cat("M_upgma_fmat dim:", paste(dim(M_upgma_fmat), collapse = "x"),
    "  range:", range(M_upgma_fmat), "\n")
cat("M_true_fmat  dim:", paste(dim(M_true_fmat),  collapse = "x"),
    "  range:", range(M_true_fmat),  "\n")
cat("upper-tri zeros in M_upgma? ",
    all(M_upgma_fmat[upper.tri(M_upgma_fmat)] == 0), "\n")

cat("\n=== Building small cache (n=500, R=3) ===\n")
t0 <- Sys.time()
cache <- build_gradient_cache(
  n_tips      = n_tips,
  coal_times  = true_coalescent_times,
  sequences   = sequences,
  n_samples   = 500,
  R_labelings = 3,
  seed        = 42
)
cat("cache build time:", format(Sys.time() - t0), "\n")
cat("loglik_F summary:\n"); print(summary(cache$loglik_F))
cat("loglik_F sd:", sd(cache$loglik_F), "\n")

cat("\n=== Probe expectations at M_init, beta=1 ===\n")
res <- expectations_under_M_beta(cache, M_upgma_fmat, beta = 1)
cat("ESS prior:", res$ess_prior, "  ESS post:", res$ess_post, "\n")
cat("log_marginal:", res$log_marginal, "\n")
cat("||E_F_post - E_F_prior||_F:",
    sqrt(sum((res$E_F_post - res$E_F_prior)^2)), "\n")
cat("E_dist2_prior:", res$E_dist2_prior,
    "  E_dist2_post:", res$E_dist2_post, "\n")
cat("grad_beta (prior - post):", res$E_dist2_prior - res$E_dist2_post, "\n")

cat("\n=== Try a single update step ===\n")
M_before <- M_upgma_fmat
grad_M <- 2 * 1 * (res$E_F_post - res$E_F_prior)
cat("||grad_M|| (lr=1, beta=1):", sqrt(sum(grad_M^2)), "\n")
cat("max |grad_M entry|:", max(abs(grad_M)), "\n")

cat("\n=== Run 20 gradient steps ===\n")
fit <- gradient_ascent_M_beta(
  cache         = cache,
  M_init        = M_upgma_fmat,
  beta_init     = 1,
  n_iter        = 20,
  lr_M          = 5e-4,
  lr_g          = 1e-2,
  verbose_every = 1
)

cat("\n=== Final diagnostics ===\n")
cat("beta start:", 1,         " end:", fit$beta, "\n")
cat("||M_continuous - M_init||_F:",
    sqrt(sum((fit$M_continuous - M_upgma_fmat)^2)), "\n")
cat("||M_projected  - M_init||_F:",
    sqrt(sum((fit$M_projected  - M_upgma_fmat)^2)), "\n")
cat("log_marginal trajectory:\n")
print(fit$history$log_marginal)
cat("ess_prior trajectory:\n"); print(fit$history$ess_prior)
cat("ess_post trajectory:\n");  print(fit$history$ess_post)
cat("grad_M_norm trajectory:\n"); print(fit$history$grad_M_norm)
