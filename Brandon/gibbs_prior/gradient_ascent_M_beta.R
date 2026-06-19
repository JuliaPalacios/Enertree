source("gibbs_prior_joint_standardized_helpers.R")

build_gradient_cache <- function(n_tips,
                                 coal_times,
                                 sequences,
                                 n_samples    = 50000,
                                 R_labelings  = 5,
                                 mode         = "average",
                                 seed         = NULL) {
  if (!is.null(seed)) set.seed(seed)

  encod_chain <- rEncod(m = n_samples, n = n_tips, distr = "uniform")
  F_list      <- lapply(encod_chain, Fmat_from_myencod)

  base_cache  <- precompute_tree_chain_distance_cache(F_list)

  cat("Pre-computing log P(S|F) for", n_samples, "F matrices...\n")
  pb <- progress_bar$new(
    format = "  loglik [:bar] :percent eta: :eta",
    total  = n_samples, clear = FALSE, width = 60, show_after = 0
  )
  loglik_F <- numeric(n_samples)
  for (i in seq_len(n_samples)) {
    loglik_F[i] <- log_likelihood_given_tree(
      tree_fmat  = F_list[[i]],
      coal_times = coal_times,
      sequences  = sequences,
      mode       = mode,
      R          = R_labelings
    )$log_likelihood
    pb$tick()
  }
  pb$terminate()

  c(base_cache, list(
    F_list      = F_list,
    loglik_F    = loglik_F,
    R_labelings = R_labelings,
    label_mode  = mode
  ))
}

# -----------------------------------------------------------------------------
# Numeric helpers.
# -----------------------------------------------------------------------------
softmax_weights <- function(log_w) {
  m <- max(log_w)
  w <- exp(log_w - m)
  w / sum(w)
}

ess_of_normalized <- function(w) 1 / sum(w * w)

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Fold a length-feature_count vector (lower-tri + diag) back to a
# lower-triangular F-matrix.  F-matrices in this codebase are lower-tri
# (zeros above the diagonal), so we just drop the values into the lower
# triangle and leave the upper triangle as 0.
expand_lower_tri_to_full <- function(vec, lower_idx, dim) {
  out <- matrix(0, dim[1], dim[2])
  out[lower_idx] <- vec
  out
}

# -----------------------------------------------------------------------------
# expectations_under_M_beta(): both expectations the gradient needs, plus
# diagnostics, from one shared uniform pool.
# -----------------------------------------------------------------------------
expectations_under_M_beta <- function(cache, M, beta, diam = 1) {
  D2     <- tree_chain_squared_distances(cache, M)          # ||F_i - M||^2
  log_pi <- -(beta / diam) * D2                              # Boltzmann factor
  log_W  <- log_pi + cache$loglik_F                          # Boltzmann * likelihood

  pi_w <- softmax_weights(log_pi)
  W    <- softmax_weights(log_W)

  # Weighted means of F over the cache, computed in lower-tri vector form.
  E_F_prior_vec <- as.numeric(crossprod(cache$tree_matrix, pi_w))
  E_F_post_vec  <- as.numeric(crossprod(cache$tree_matrix, W))

  E_F_prior <- expand_lower_tri_to_full(E_F_prior_vec, cache$lower_idx, cache$template_dim)
  E_F_post  <- expand_lower_tri_to_full(E_F_post_vec,  cache$lower_idx, cache$template_dim)

  log_Z_prior <- logsumexp(log_pi) - log(length(log_pi))     # log of (1/N) sum exp(...)
  log_marg    <- logsumexp(log_W)  - logsumexp(log_pi)       # log P(S | M, beta) up to const

  list(
    E_F_prior     = E_F_prior,
    E_F_post      = E_F_post,
    E_dist2_prior = sum(pi_w * D2),
    E_dist2_post  = sum(W    * D2),
    ess_prior     = ess_of_normalized(pi_w),
    ess_post      = ess_of_normalized(W),
    log_Z_est     = log_Z_prior,
    log_marginal  = log_marg
  )
}

# -----------------------------------------------------------------------------
# gradient_ascent_M_beta(): SGA on (M, g = log beta).
#
# Step on g = log beta keeps beta > 0 automatically; chain rule gives
#   dL/dg = beta * dL/dbeta.
#
# Projection onto the valid F-matrix cone is OFF by default — gradient ascent
# moves M through the convex hull of F-space, and projecting every step over-
# constrains. Project once at the end if you want a discrete F.
# -----------------------------------------------------------------------------
gradient_ascent_M_beta <- function(
    cache,
    M_init,
    beta_init      = 1,
    n_iter         = 500,
    lr_M           = 1e-3,
    lr_g           = 1e-2,
    project_M      = FALSE,
    ess_warn_thres = 50,
    verbose_every  = 25) {

  M <- M_init
  g <- log(beta_init)

  history <- list(
    M            = vector("list", n_iter),
    beta         = numeric(n_iter),
    log_marginal = numeric(n_iter),
    ess_prior    = numeric(n_iter),
    ess_post     = numeric(n_iter),
    grad_M_norm  = numeric(n_iter),
    grad_g       = numeric(n_iter)
  )

  for (k in seq_len(n_iter)) {
    beta <- exp(g)

    res <- expectations_under_M_beta(cache, M, beta)

    grad_M    <- 2 * beta * (res$E_F_post - res$E_F_prior)
    grad_beta <- res$E_dist2_prior - res$E_dist2_post
    grad_g    <- beta * grad_beta                            # chain rule

    M <- M + lr_M * grad_M
    g <- g + lr_g * grad_g
    if (project_M) M <- nearby_Fmat(M)

    history$M[[k]]          <- M
    history$beta[k]         <- beta
    history$log_marginal[k] <- res$log_marginal
    history$ess_prior[k]    <- res$ess_prior
    history$ess_post[k]     <- res$ess_post
    history$grad_M_norm[k]  <- sqrt(sum(grad_M * grad_M))
    history$grad_g[k]       <- grad_g

    if (res$ess_post < ess_warn_thres || res$ess_prior < ess_warn_thres) {
      warning(sprintf(
        "iter %d: low ESS (prior=%.1f, post=%.1f) - enlarge cache or temper likelihood",
        k, res$ess_prior, res$ess_post
      ))
    }

    if (k %% verbose_every == 0 || k == 1) {
      cat(sprintf(
        "iter %4d | beta = %7.4f | ESS pr/po = %6.1f / %6.1f | log P(S|M,b) = %9.3f | |gM| = %.3e\n",
        k, beta, res$ess_prior, res$ess_post,
        res$log_marginal, sqrt(sum(grad_M * grad_M))
      ))
    }
  }

  M_final <- if (project_M) M else nearby_Fmat(M)            # also report projected
  list(
    M_continuous = M,
    M_projected  = M_final,
    beta         = exp(g),
    history      = history
  )
}
library(ape)
library(phylodyn)
library(phangorn)
library(phyclust)
library(phylotools)
library(fmatrix)
library(coda)
library(ggplot2)

source("utils.R")


n_tips <- 25
set.seed(123)

sequences <- as.phyDat(read.FASTA("sequences.fasta"))

M_true_tree <- read.tree("true_tree.newick")
inter_coal_times <- coalescent.intervals(M_true_tree)$interval.length
inter_coal_times[which(inter_coal_times <= 0.000001)] <- 0.0000001
true_coalescent_times <- cumsum(inter_coal_times)
M_true_fmat <- gen_Fmat(M_true_tree, tol = 8)

dm <- dist.hamming(sequences)

upgma_tree <- upgma(dm)
upgma_tree <- update_time(upgma_tree, true_coalescent_times)
M_upgma_fmat <- gen_Fmat(upgma_tree, tol = 8)
# -----------------------------------------------------------------------------

 cache <- build_gradient_cache(
   n_tips      = n_tips,
   coal_times  = true_coalescent_times,
   sequences   = sequences,
   n_samples   = 500000,
   R_labelings = 10,
   seed        = 42
 )

 

fit <- gradient_ascent_M_beta(
   cache     = cache,
   M_init    = M_upgma_fmat,
   beta_init = 1,
   n_iter    = 500,
   lr_M      = 5e4,
   lr_g      = 1e2xs
 )

# plot(fit$history$log_marginal, type = "l",
#      xlab = "iter", ylab = "log P(S | M, beta) (up to constant)")
# plot(fit$history$beta,         type = "l", xlab = "iter", ylab = "beta")
# plot(fit$history$ess_post,     type = "l", xlab = "iter", ylab = "ESS posterior")
# -----------------------------------------------------------------------------
