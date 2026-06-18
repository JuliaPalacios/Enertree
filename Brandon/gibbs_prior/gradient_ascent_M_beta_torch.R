suppressPackageStartupMessages({
  library(torch)
})

source("gibbs_prior_joint_standardized_helpers.R")
source("gradient_ascent_M_beta.R")   # for build_gradient_cache, expand_lower_tri_to_full

# -----------------------------------------------------------------------------
# Convert a build_gradient_cache() result into torch tensors that the
# autodiff loop will consume.  Only the cached F's and log-likelihoods need to
# become tensors -- everything (M, beta) on the optimization side is a tensor
# created by fit_M_beta_torch().
# -----------------------------------------------------------------------------
cache_to_torch <- function(cache, dtype = torch_float64()) {
  list(
    F_t        = torch_tensor(cache$tree_matrix, dtype = dtype),   # n x feat
    ll_t       = torch_tensor(cache$loglik_F,    dtype = dtype),   # n
    lower_idx  = cache$lower_idx,
    dim        = cache$template_dim,
    feat_count = sum(cache$lower_idx)
  )
}

# -----------------------------------------------------------------------------
# The differentiable log-marginal estimator.
#
#   log P(y | M, beta)  ~=  logsumexp_i(ell_i - beta ||F_i - M||^2)
#                          - logsumexp_i(    -beta ||F_i - M||^2)
#
# Inputs M_vec_t and log_beta_t must have requires_grad = TRUE for autodiff.
# Everything else is constant w.r.t. the optimization variables.
# -----------------------------------------------------------------------------
torch_log_marginal <- function(M_vec_t, log_beta_t, F_t, ll_t) {
  beta   <- log_beta_t$exp()
  diffs  <- F_t - M_vec_t$unsqueeze(1)              # broadcast: n x feat
  D2     <- diffs$pow(2)$sum(dim = 2)               # length n
  log_pi <- -beta * D2                              # length n
  log_W  <- log_pi + ll_t                           # length n
  torch_logsumexp(log_W, dim = 1) - torch_logsumexp(log_pi, dim = 1)
}

# Helper: effective sample size from a vector of log-weights (no grad needed).
torch_ess <- function(log_w) {
  W <- nnf_softmax(log_w, dim = 1)
  as.numeric(1 / W$pow(2)$sum()$item())
}

# -----------------------------------------------------------------------------
# fit_M_beta_torch(): autodiff gradient ascent on (M_vec, log beta).
#
#   optimizer = "adam"  -> elementwise adaptive step sizes; robust, no line
#                          search; good default.
#   optimizer = "lbfgs" -> quasi-Newton with strong-Wolfe line search; usually
#                          fewer iterations but each step is more expensive.
#
# Parameterisation:
#   - M is parameterised by its lower-triangle vector (length feat_count) so
#     gradient updates can't accidentally fill the upper triangle.
#   - beta lives as g = log(beta) so optimization is unconstrained but
#     beta = exp(g) is always positive.
# -----------------------------------------------------------------------------
fit_M_beta_torch <- function(cache,
                             M_init,
                             beta_init     = 1,
                             optimizer     = c("adam", "lbfgs"),
                             lr            = 1e-2,
                             n_iter        = 500,
                             verbose_every = 25,
                             dtype         = torch_float64()) {
  optimizer <- match.arg(optimizer)
  tc <- cache_to_torch(cache, dtype = dtype)

  # Optimization variables (requires_grad = TRUE so autograd tracks them).
  M_vec    <- torch_tensor(as.numeric(M_init[tc$lower_idx]),
                           dtype = dtype, requires_grad = TRUE)
  log_beta <- torch_tensor(log(beta_init),
                           dtype = dtype, requires_grad = TRUE)

  history <- list(
    log_marginal = numeric(n_iter),
    beta         = numeric(n_iter),
    ess_prior    = numeric(n_iter),
    ess_post     = numeric(n_iter)
  )

  step_and_log <- function(k) {
    # Diagnostics under no-grad to avoid building unneeded graph nodes.
    with_no_grad({
      beta   <- log_beta$exp()
      diffs  <- tc$F_t - M_vec$unsqueeze(1)
      D2     <- diffs$pow(2)$sum(dim = 2)
      log_pi <- -beta * D2
      log_W  <- log_pi + tc$ll_t
      history$ess_prior[k] <<- torch_ess(log_pi)
      history$ess_post[k]  <<- torch_ess(log_W)
      history$beta[k]      <<- log_beta$exp()$item()
    })
    if (k %% verbose_every == 0 || k == 1) {
      cat(sprintf(
        "iter %4d | beta = %7.4f | log P(y|M,b) = %9.3f | ESS pr/po = %6.1f / %6.1f\n",
        k, history$beta[k], history$log_marginal[k],
        history$ess_prior[k], history$ess_post[k]
      ))
    }
  }

  if (optimizer == "adam") {
    opt <- optim_adam(list(M_vec, log_beta), lr = lr)
    for (k in seq_len(n_iter)) {
      opt$zero_grad()
      lm   <- torch_log_marginal(M_vec, log_beta, tc$F_t, tc$ll_t)
      loss <- -lm
      loss$backward()
      opt$step()

      history$log_marginal[k] <- lm$item()
      step_and_log(k)
    }
  } else {
    # L-BFGS with line search; closure pattern required by torch.
    opt <- optim_lbfgs(list(M_vec, log_beta),
                       lr             = lr,
                       max_iter       = 20,
                       history_size   = 10,
                       line_search_fn = "strong_wolfe")
    for (k in seq_len(n_iter)) {
      closure <- function() {
        opt$zero_grad()
        loss <- -torch_log_marginal(M_vec, log_beta, tc$F_t, tc$ll_t)
        loss$backward()
        loss
      }
      loss <- opt$step(closure)
      history$log_marginal[k] <- -loss$item()
      step_and_log(k)
    }
  }

  # Materialize results back to plain R matrices.
  M_vec_final  <- as.numeric(M_vec$detach())
  M_continuous <- expand_lower_tri_to_full(M_vec_final, tc$lower_idx, tc$dim)
  M_projected  <- nearby_Fmat(M_continuous)

  list(
    M_continuous = M_continuous,
    M_projected  = M_projected,
    beta         = log_beta$exp()$item(),
    history      = history
  )
}

# -----------------------------------------------------------------------------


cache <- build_gradient_cache(
  n_tips      = n_tips,
  coal_times  = true_coalescent_times,
  sequences   = sequences,
  n_samples   = 50000,
  R_labelings = 10,
  seed        = 42
)

fit_adam <- fit_M_beta_torch(
  cache         = cache,
  M_init        = M_upgma_fmat,
  beta_init     = 1,
  optimizer     = "adam",
  lr            = 1e-2,
  n_iter        = 500,
  verbose_every = 10
)

# Or quasi-Newton with line search:
fit_lbfgs <- fit_M_beta_torch(
  cache         = cache,
  M_init        = M_upgma_fmat,
  beta_init     = 1,
  optimizer     = "lbfgs",
  lr            = 1,
  n_iter        = 50,
  verbose_every = 1
)

plot(fit_adam$history$log_marginal, type = "l",
     xlab = "iter", ylab = "log P(y | M, beta)")
plot(fit_adam$history$beta,         type = "l", xlab = "iter", ylab = "beta")
plot(fit_adam$history$ess_post,     type = "l", xlab = "iter", ylab = "ESS posterior")
# -----------------------------------------------------------------------------
