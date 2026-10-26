# Exact six-tip fixed-time posterior comparison and four-panel diagnostics.
# Rscript Brandon/gibbs_prior/exact_compare_VB_joint_FA.R [plot.png]
fa_exact_source_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
if (is.null(fa_exact_source_file)) {
  fa_exact_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(fa_exact_arg)) fa_exact_source_file <- sub("^--file=", "", fa_exact_arg[1L])
}
if (is.null(fa_exact_source_file)) {
  fa_exact_source_file <- if (file.exists("Brandon/gibbs_prior/exact_compare_VB_joint_FA.R"))
    "Brandon/gibbs_prior/exact_compare_VB_joint_FA.R" else "exact_compare_VB_joint_FA.R"
}
fa_exact_dir <- dirname(normalizePath(fa_exact_source_file, mustWork = TRUE))
source(file.path(fa_exact_dir, "VB_joint_FA.R"))

fa_exact_probabilities <- function(score) exp(score - fa_logsumexp(score))
fa_exact_kl <- function(q, logp) {
  keep <- q > 0
  sum(q[keep] * (log(q[keep]) - logp[keep]))
}
fa_exact_spread <- function(q, A_values) {
  inclusion <- as.numeric(crossprod(q, A_values))
  sum(2 * inclusion * (1 - inclusion))
}

plot_exact_compare_VB_joint_FA <- function(result, plot_file) {
  stopifnot(is.character(plot_file), length(plot_file) == 1L, nzchar(plot_file))
  dir.create(dirname(plot_file), recursive = TRUE, showWarnings = FALSE)
  p <- result$posterior
  models <- result$models
  labels <- c("Uniform\nprior", "F only\n(uniform labels)",
              "Initial A\nweights", "Learned\nF + A")
  colors <- c("#999999", "#CF6843", "#8063A8", "#1681A1")
  kl <- vapply(models, fa_exact_kl, numeric(1), logp = result$logposterior)
  spread <- c(posterior = fa_exact_spread(p, result$A_values),
              vapply(models, fa_exact_spread, numeric(1), A_values = result$A_values))

  grDevices::png(plot_file, width = 2200, height = 1600, res = 200)
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mfrow = c(2, 2), mar = c(5.5, 5.1, 3.2, 1.2), oma = c(0.3, 0, 0, 0))

  barplot(kl, names.arg = labels, col = colors, border = "#333333",
          ylim = c(0, max(kl) * 1.12), cex.names = 0.78,
          main = "Approximation error", ylab = "KL(q || exact posterior) (nats)")
  box(bty = "l")

  plot(result$exact_clade, result$learned_clade, pch = 19,
       col = grDevices::adjustcolor(colors[4], alpha.f = 0.75),
       xlim = c(0, 1), ylim = c(0, 1), asp = 1,
       xlab = "Exact posterior clade probability",
       ylab = "Joint F + A clade probability",
       main = "Clade support: learned F + A")
  abline(0, 1, lty = 2, col = "#555555")
  box(bty = "l")

  top <- head(order(p, decreasing = TRUE), 10L)
  barplot(rbind("Exact posterior" = p[top],
                "Learned F + A" = models$learned[top]),
          beside = TRUE, names.arg = seq_along(top),
          col = c("#333333", colors[4]), border = NA,
          ylim = c(0, max(c(p[top], models$learned[top])) * 1.1),
          xlab = "Genealogy, ordered by exact posterior probability",
          ylab = "Probability", main = "Top ten labeled ranked genealogies")
  legend("topright", legend = c("Exact posterior", "Learned F + A"),
         fill = c("#333333", colors[4]), bty = "n", cex = 0.8)
  box(bty = "l")

  barplot(spread, names.arg = c("Posterior", labels),
          col = c("#333333", colors), border = "#333333",
          ylim = c(0, max(spread) * 1.12), cex.names = 0.70,
          ylab = "Expected number of differing A entries",
          main = "Structural spread: two independent draws")
  box(bty = "l")
  invisible(plot_file)
}

run_exact_compare_VB_joint_FA <- function(
    steps = 25L, seed = 42L,
    plot_file = file.path(fa_exact_dir, "joint_FA_exact_diagnostics.png")) {
  previous <- readRDS(file.path(fa_source_dir, "logs", "VB_master_A",
                                "20261024_184601", "result.rds"))
  labels <- names(previous$sequences)[1:6]
  fasta <- ape::read.FASTA(file.path(fa_source_dir, "A matrix analysis", "outputs",
                                     "20261024_184601", "sequences.fasta"))
  sequences <- phangorn::as.phyDat(as.matrix(fasta)[labels, , drop = FALSE])
  subtree <- ape::drop.tip(previous$true_tree,
                           setdiff(previous$true_tree$tip.label, labels))
  coal_times <- cumsum(ape::coalescent.intervals(subtree)$interval.length)
  stopifnot(length(coal_times) == 5L, all(diff(coal_times) > 0), coal_times[1] > 0)

  states <- fa_enumerate(labels)
  loglik <- vapply(states, fa_likelihood, numeric(1), sequences = sequences,
                   coal_times = coal_times, rate = previous$config$rate)
  logposterior <- loglik - fa_logsumexp(loglik)
  posterior <- exp(logposterior)
  fit <- fa_fit(sequences, coal_times, rate = previous$config$rate,
                steps = steps, samples_per_chain = 100L, chains = 2L,
                burnin = 300L, thin = 8L, reference_n = 2700L, seed = seed)
  A_values <- t(vapply(states, as.vector, numeric(length(states[[1L]]))))
  d2 <- vapply(states, function(A) sum((fa_F(A, coal_times) - fit$M)^2), numeric(1))
  F_score <- -exp(fit$g) * d2 - vapply(states, fa_log_multiplicity, numeric(1))
  theta_initial <- 0.3 * (2 * vb_A_from_tree(
    phylodyn:::update_time(phangorn::upgma(phangorn::dist.hamming(sequences)),
                            coal_times), names(sequences)) - 1)
  theta_initial[1L, ] <- 0
  score <- F_score + as.numeric(A_values %*% as.vector(fit$theta))
  stopifnot(abs(score[1L] - fa_score(states[[1L]], fit$M, fit$g,
                                    fit$theta, coal_times)) < 1e-10)
  models <- list(prior = rep(1 / length(states), length(states)),
                 F_only = fa_exact_probabilities(F_score),
                 initial_A = fa_exact_probabilities(
                   F_score + as.numeric(A_values %*% as.vector(theta_initial))),
                 learned = fa_exact_probabilities(score))
  variational <- models$learned

  clades <- lapply(states, function(A) vapply(seq.int(2L, nrow(A)),
    function(i) paste(sort(labels[A[i, ] == 1]), collapse = "|"), ""))
  keys <- sort(unique(unlist(clades, use.names = FALSE)))
  incidence <- t(vapply(clades, function(x) as.integer(keys %in% x),
                        integer(length(keys))))
  exact_clade <- as.numeric(crossprod(posterior, incidence))
  variational_clade <- as.numeric(crossprod(variational, incidence))
  summary <- data.frame(n = 6L, states = length(states),
                        total_variation = sum(abs(posterior - variational)) / 2,
                        clade_mae = mean(abs(exact_clade - variational_clade)),
                        kl_to_posterior = fa_exact_kl(variational, logposterior),
                        exact_logZ = fa_logsumexp(score),
                        ais_logZ = fit$ais$logZ,
                        ais_importance_ess = fit$ais$importance_ess)
  result <- list(summary = summary, fit = fit, posterior = posterior,
                 logposterior = logposterior,
                 models = models, A_values = A_values,
                 exact_clade = exact_clade, learned_clade = variational_clade,
                 clade_keys = keys, plot_file = plot_file)
  if (!is.null(plot_file)) plot_exact_compare_VB_joint_FA(result, plot_file)
  result
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  result <- run_exact_compare_VB_joint_FA(
    plot_file = if (length(args)) args[1L] else
      file.path(fa_exact_dir, "joint_FA_exact_diagnostics.png"))
  print(result$summary, row.names = FALSE)
  cat("Plot:", normalizePath(result$plot_file, mustWork = TRUE), "\n")
}
