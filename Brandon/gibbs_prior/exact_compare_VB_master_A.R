# Exact fixed-time posterior comparison for a saved six-tip A-only fit.
# From the repository root:
# Rscript Brandon/gibbs_prior/exact_compare_VB_master_A.R [result.rds] [plot.png]
# The script reads the fit's own sequences, times, rate, and A weights. It does
# not refit the model or use BEAST trees.

exact_A_source <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
if (is.null(exact_A_source)) {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg)) exact_A_source <- sub("^--file=", "", file_arg[1L])
}
if (is.null(exact_A_source)) {
  exact_A_source <- if (file.exists("Brandon/gibbs_prior/exact_compare_VB_master_A.R"))
    "Brandon/gibbs_prior/exact_compare_VB_master_A.R" else "exact_compare_VB_master_A.R"
}
exact_A_dir <- dirname(normalizePath(exact_A_source, mustWork = TRUE))
source(file.path(exact_A_dir, "VB_master_A.R"))

exact_A_logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

exact_A_enumerate <- function(labels) {
  n <- length(labels)
  if (n < 2L || n > 6L || anyNA(labels) || anyDuplicated(labels)) {
    stop("Exact enumeration supports 2 to 6 uniquely named tips.")
  }
  active <- lapply(seq_len(n), function(j) j)
  A <- matrix(0L, n - 1L, n,
              dimnames = list(paste0("rank_", seq_len(n - 1L)), labels))
  states <- list()
  visit <- function(active, A) {
    k <- length(active)
    if (k == 1L) {
      states[[length(states) + 1L]] <<- A
      return(invisible(NULL))
    }
    for (pair in combn(k, 2L, simplify = FALSE)) {
      merged <- sort(c(active[[pair[1L]]], active[[pair[2L]]]))
      B <- A
      B[k - 1L, merged] <- 1L
      visit(c(active[-pair], list(merged)), B)
    }
    invisible(NULL)
  }
  visit(active, A)
  states
}

exact_A_kl <- function(q, logp) {
  keep <- q > 0
  sum(q[keep] * (log(q[keep]) - logp[keep]))
}

exact_A_spread <- function(q, values) {
  inclusion <- as.numeric(crossprod(q, values))
  sum(2 * inclusion * (1 - inclusion))
}

exact_A_find_fit <- function() {
  files <- Sys.glob(file.path(exact_A_dir, "logs", "VB_master_A", "*", "result.rds"))
  if (!length(files)) stop("No saved A-only result.rds found; supply one as the first argument.")
  files <- files[order(file.info(files)$mtime, decreasing = TRUE)]
  for (file in files) {
    fit <- readRDS(file)
    if (length(fit$sequences) == 6L) return(file)
  }
  stop("No saved six-tip A-only fit found; run VB_master_A.R with num_tips=6 first.")
}

plot_exact_compare_VB_master_A <- function(result, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  posterior <- result$posterior
  models <- result$models
  labels <- c("Uniform\nprior", "Initial A\nweights", "Learned A\nweights")
  colors <- c("#999999", "#8063A8", "#1681A1")
  kl <- vapply(models, exact_A_kl, numeric(1), logp = result$logposterior)
  spread <- c(posterior = exact_A_spread(posterior, result$A_values),
              vapply(models, exact_A_spread, numeric(1), values = result$A_values))

  grDevices::png(file, width = 2200, height = 1600, res = 200)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- par(mfrow = c(2, 2), mar = c(5.5, 5.1, 3.2, 1.2))
  on.exit(par(old_par), add = TRUE)

  barplot(kl, names.arg = labels, col = colors, border = "#333333",
          ylim = c(0, max(kl) * 1.12), cex.names = 0.82,
          main = "Approximation error", ylab = "KL(q || exact posterior) (nats)")
  box(bty = "l")

  plot(result$exact_clade, result$learned_clade, pch = 19,
       col = grDevices::adjustcolor(colors[3L], alpha.f = 0.75),
       xlim = c(0, 1), ylim = c(0, 1), asp = 1,
       xlab = "Exact posterior clade probability",
       ylab = "A-only variational clade probability",
       main = "Clade support: learned A weights")
  abline(0, 1, lty = 2, col = "#555555")
  box(bty = "l")

  top <- head(order(posterior, decreasing = TRUE), 10L)
  barplot(rbind("Exact posterior" = posterior[top],
                "Learned A weights" = models$learned[top]),
          beside = TRUE, names.arg = seq_along(top),
          col = c("#333333", colors[3L]), border = NA,
          ylim = c(0, max(c(posterior[top], models$learned[top])) * 1.1),
          xlab = "Genealogy, ordered by exact posterior probability",
          ylab = "Probability", main = "Top ten labeled ranked genealogies")
  legend("topright", legend = c("Exact posterior", "Learned A weights"),
         fill = c("#333333", colors[3L]), bty = "n", cex = 0.8)
  box(bty = "l")

  barplot(spread, names.arg = c("Posterior", labels),
          col = c("#333333", colors), border = "#333333",
          ylim = c(0, max(spread) * 1.12), cex.names = 0.75,
          ylab = "Expected number of differing A entries",
          main = "Structural spread: two independent draws")
  box(bty = "l")
  invisible(file)
}

run_exact_compare_VB_master_A <- function(
    result_file = exact_A_find_fit(),
    plot_file = file.path(dirname(result_file), "exact_posterior_comparison.png")) {
  fit <- readRDS(result_file)
  labels <- names(fit$sequences)
  n <- length(labels)
  if (n != 6L || !identical(colnames(fit$theta), labels) ||
      !identical(dim(fit$theta), c(n - 1L, n)) ||
      !identical(dim(fit$theta_initial), dim(fit$theta)) ||
      length(fit$coal_times) != n - 1L || !is.finite(fit$config$rate)) {
    stop("The saved result must contain a consistent six-tip A-only fit.")
  }
  states <- exact_A_enumerate(labels)
  stopifnot(length(states) == 2700L,
            all(vapply(states, phylodyn:::a_matrix_valid, logical(1L))),
            !anyDuplicated(vapply(states, phylodyn:::a_matrix_key, character(1L))))
  A_values <- t(vapply(states, as.vector, numeric(length(fit$theta))))
  loglik <- vapply(states, function(A) {
    tree <- vb_A_to_tree(A, fit$coal_times, labels, check = FALSE)
    as.numeric(phangorn::pml(tree, fit$sequences, bf = rep(0.25, 4L),
                            Q = c(1, 2, 1, 1, 2, 1), rate = fit$config$rate)$logLik)
  }, numeric(1L))
  if (any(!is.finite(loglik))) stop("Non-finite exact likelihood.")
  logposterior <- loglik - exact_A_logsumexp(loglik)
  posterior <- exp(logposterior)
  initial_score <- as.numeric(A_values %*% as.vector(fit$theta_initial))
  learned_score <- as.numeric(A_values %*% as.vector(fit$theta))
  models <- list(prior = rep(1 / length(states), length(states)),
                 initial_A = exp(initial_score - exact_A_logsumexp(initial_score)),
                 learned = exp(learned_score - exact_A_logsumexp(learned_score)))
  clades <- lapply(states, function(A) vapply(seq.int(2L, nrow(A)),
    function(i) paste(sort(labels[A[i, ] == 1L]), collapse = "|"), ""))
  keys <- sort(unique(unlist(clades, use.names = FALSE)))
  incidence <- t(vapply(clades, function(x) as.integer(keys %in% x),
                        integer(length(keys))))
  exact_clade <- as.numeric(crossprod(posterior, incidence))
  learned_clade <- as.numeric(crossprod(models$learned, incidence))
  summary <- data.frame(n = n, states = length(states),
                        beta_init = fit$config$beta_init,
                        total_variation = sum(abs(posterior - models$learned)) / 2,
                        clade_mae = mean(abs(exact_clade - learned_clade)),
                        kl_prior = exact_A_kl(models$prior, logposterior),
                        kl_initial = exact_A_kl(models$initial_A, logposterior),
                        kl_to_posterior = exact_A_kl(models$learned, logposterior),
                        exact_logZ_A = exact_A_logsumexp(learned_score))
  result <- list(summary = summary, fit_file = result_file, plot_file = plot_file,
                 posterior = posterior, logposterior = logposterior,
                 models = models, A_values = A_values, clade_keys = keys,
                 exact_clade = exact_clade, learned_clade = learned_clade)
  if (!is.null(plot_file)) plot_exact_compare_VB_master_A(result, plot_file)
  result
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  fit_file <- if (length(args)) args[1L] else exact_A_find_fit()
  plot_file <- if (length(args) >= 2L) args[2L] else
    file.path(dirname(fit_file), "exact_posterior_comparison.png")
  result <- run_exact_compare_VB_master_A(fit_file, plot_file)
  print(result$summary, row.names = FALSE)
  cat("Fit:", normalizePath(result$fit_file, mustWork = TRUE), "\n")
  cat("Plot:", normalizePath(result$plot_file, mustWork = TRUE), "\n")
}
