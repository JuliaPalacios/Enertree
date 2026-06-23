library(progress)
library(ggplot2)


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
  seqgen(opts="-mHKY -t2.0 -f0.25,0.25,0.25,0.25 -l1000 -s0.01",
         newick.tree=M_true_newick_string,
         temp.file="temp.fasta")
  data2<-read.phylip("temp.fasta")
  
  data<-dat2fasta(data2,outfile="sequences.fasta")
  mydna<-as.phyDat(read.FASTA("sequences.fasta"))
  write.tree(M_true_tree$newick, file="true_tree.newick")
  return(list(M_true_tree = M_true_tree$newick,
              sequences = mydna))
}





ZigZag <- function(n) {
  fact <- list()
  zig <- list()
  fact[[1]] <- 1
  for (i in 1:n) {
    fact[[i + 1]] <- fact[[i]] * i
  }

  zig[[1]] <- 1
  zig[[2]] <- 1

  for (i in 2:(n - 1)) {
    sum <- 0
    for (k in 0:(i - 1)) {
      sum <- sum + (fact[[i]] / (fact[[i - k]] * fact[[k + 1]])) * zig[[k + 1]] * zig[[i - k]]
    }
    zig[[i + 1]] <- sum / 2
  }

  zig[[length(zig)]]
}



Logdistance_prob <- function(fmat, M, beta=1, diam=1, N=1){
  n <- dim(fmat)[1] + 1
  -(beta/diam) * distance_Fmat(fmat, M, dist="l2")^2/N^4
}


log_likelihood_given_tree <- function(tree_fmat, coal_times, sequences, mode = "max", R = 1) {
  best_ll <- -Inf
  best_tree <- NULL
  ll_sum <- 0
  for (r in seq_len(R)) {
    rooted_tree <- mytree_from_F(tree_fmat, coal_times)
    ll <- pml(
      rooted_tree,
      sequences,
      bf = c(0.25, 0.25, 0.25, 0.25),
      Q  = c(1, 2, 1, 1, 2, 1)
    )$log
    ll_sum <- ll_sum + ll
    if (ll > best_ll) {
      best_ll <- ll
      best_tree <- rooted_tree
    }
  }
  log_likelihood <- if (mode == "max") best_ll else ll_sum / R
  list(log_likelihood = log_likelihood, rooted_tree = best_tree)
}


sampleF <- function(M, beta, iter, diam, startF = NULL, take_every = 1){
  acceptance_rate <- 0
  n <- nrow(M) + 1

  chainF <- vector("list", iter/take_every)

  if (is.null(startF)) {
    current_tree <- rcoal(n)
    chainF[[1]]  <- gen_Fmat(current_tree, tol=13)
  } else {
    chainF[[1]] <- startF
  }

  current_Encod <- my_encod(chainF[[1]])
  f_current     <- Fmat_from_myencod(current_Encod)
  Logdist_curr  <- Logdistance_prob(f_current, M, beta=beta, diam=diam)

  MeanM <- matrix(0, nrow=n-1, ncol=n-1)
  
  
  for (j in 1:iter){
    
    proposed    <- proposal_myencod(current_Encod)
    f_proposed  <- Fmat_from_myencod(proposed)
    Logdist_prop <- Logdistance_prob(f_proposed, M, beta=beta, diam=diam)

    prob <- exp(Logdist_prop - Logdist_curr)

    if (runif(1) < prob){
      if(j%%take_every == 0){
        chainF[[j/take_every]] <- f_proposed
      }
      
      current_Encod   <- proposed
      f_current       <- f_proposed
      Logdist_curr    <- Logdist_prop
      acceptance_rate <- acceptance_rate + 1
    } else {
      if(j%%take_every == 0){
        chainF[[j/take_every]] <- f_current
      }
    }
    if(j%%take_every == 0){
       MeanM <- MeanM + chainF[[j/take_every]]
      }
    
  }

  MeanM <- MeanM/iter
  list(
    chainF = chainF,
    acceptance_rate = acceptance_rate/iter,
    MeanM = MeanM
  )
}




tree_MDS_comparison_plot <- function(trees_gen, trees_data, true_tree, M_tree = NULL, title = "Tree Comparison", save_path = NULL) {
  # 1) Combine trees
  if (!is.null(M_tree)) {
    trees_all <- c(trees_data, trees_gen, list(true_tree), list(M_tree))
  } else {
    trees_all <- c(trees_data, trees_gen, list(true_tree))
  }
  n <- length(trees_all)  # should be 2001 if 1000+1000+1

  # 2) Compute full symmetric distance matrix
  D <- matrix(0, n, n)
  
  pb <- progress_bar$new(
    format = "  computing distances [:bar] :percent eta: :eta",
    total = n - 1,
    clear = FALSE,
    width = 60,
    show_after = 0
  )
  
  for (i in 2:n) {
    ti <- trees_all[[i]]
    for (j in 1:(i - 1)) {
      d <- distance_Fmat(ti, trees_all[[j]], dist="l2")^2
      D[i, j] <- d
      D[j, i] <- d
    }
    pb$tick()
  }
  pb$terminate()
  
  # 3) Convert to 'dist' and run MDS
  D_dist <- as.dist(D)
  mds <- cmdscale(D_dist, k = 2, eig = TRUE)
  
  coords <- as.data.frame(mds$points)
  colnames(coords) <- c("Dim1", "Dim2")
  
  # 4) Add group labels
  sub_samp <- length(trees_data)
  if (!is.null(M_tree)) {
    labels <- c(rep("Data", sub_samp), rep("Generated", sub_samp), "True", "M Tree")
    coords$group <- factor(labels, levels = c("Data", "Generated", "True", "M Tree"))
  } else {
  labels <- c(rep("Data", sub_samp), rep("Generated", sub_samp), "True")
  coords$group <- factor(labels, levels = c("Data", "Generated", "True"))
  }
  
  
  
  # 5) Compute proportion of variance captured
  pos_eig <- mds$eig[mds$eig > 0]
  prop_2d <- sum(pos_eig[1:2]) / sum(pos_eig)
  print(paste("Proportion of variance captured in 2D:", round(prop_2d, 4)))
  
  # 6) Plot with ggplot2
  if (is.null(M_tree)) {
    p <- ggplot(coords, aes(x = Dim1, y = Dim2, color = group)) +
          geom_point(size = 2) +
          scale_color_manual(values = c("steelblue", "tomato", "green", "purple")) +
          labs(title = title, x = "MDS 1", y = "MDS 2", color = "Group") +
          theme_minimal(base_size = 14)
  }else {
    p <- ggplot(coords, aes(x = Dim1, y = Dim2, color = group)) +
          geom_point(size = 2) +
          scale_color_manual(values = c("steelblue", "tomato", "green", "purple")) +
          labs(title = title, x = "MDS 1", y = "MDS 2", color = "Group") +
          theme_minimal(base_size = 14)
  }

  
  print(p)
  
  # 7) Save if path provided
  if (!is.null(save_path)) {
    ggsave(filename = paste0("plots/", save_path), plot = p, width = 8, height = 6, dpi = 150)
  }
  }


matrix_MDS_comparison_plot <- function(matrices, labels, chain_matrices = NULL, sub_samp = NULL, seed = 42,
                                       title = "MDS Comparison", save_path = NULL,
                                       colors = c("green", "steelblue", "tomato", "grey70"),
                                       color_by_log_likelihood = FALSE, coal_times = NULL,
                                       sequences = NULL, log_likelihood_mode = "average",
                                       log_likelihood_R = 10) {
  if (length(matrices) != length(labels)) {
    stop("matrices and labels must have the same length")
  }

  if (color_by_log_likelihood) {
    if (is.null(coal_times) || is.null(sequences)) {
      stop("coal_times and sequences must be provided when color_by_log_likelihood = TRUE")
    }
  }

  chain_samples <- list()
  if (!is.null(chain_matrices)) {
    chain_samples <- Filter(Negate(is.null), chain_matrices)
    if (!is.null(sub_samp)) {
      sub_samp <- min(sub_samp, length(chain_samples))
      set.seed(seed)
      idx <- sample(length(chain_samples), sub_samp)
      chain_samples <- chain_samples[idx]
    }
  }

  if (color_by_log_likelihood) {
    chain_samples <- lapply(chain_samples, nearby_Fmat)
    matrices <- lapply(matrices, nearby_Fmat)
  }

  all_matrices <- c(chain_samples, matrices)
  n <- length(all_matrices)
  D <- matrix(0, n, n)

  pb <- progress_bar$new(
    format = "  computing matrix distances [:bar] :percent eta: :eta",
    total = max(1, n - 1),
    clear = FALSE,
    width = 60,
    show_after = 0
  )

  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      D[i, j] <- distance_Fmat(all_matrices[[i]], all_matrices[[j]], dist = "l2")^2
      D[j, i] <- D[i, j]
    }
    pb$tick()
  }
  pb$terminate()

  mds <- cmdscale(as.dist(D), k = 2, eig = TRUE)
  coords <- as.data.frame(mds$points)
  colnames(coords) <- c("Dim1", "Dim2")

  if (length(chain_samples) > 0) {
    group_labels <- c(rep("M Chain", length(chain_samples)), labels)
    group_levels <- c("M Chain", labels)
  } else {
    group_labels <- labels
    group_levels <- labels
  }
  coords$group <- factor(group_labels, levels = group_levels)
  coords$label <- c(rep("", max(0, length(chain_samples))), labels)

  if (color_by_log_likelihood) {
    ll_pb <- progress_bar$new(
      format = "  computing log likelihoods [:bar] :percent eta: :eta",
      total = n,
      clear = FALSE,
      width = 60,
      show_after = 0
    )

    coords$log_likelihood <- vapply(
      all_matrices,
      function(tree_fmat) {
        ll_pb$tick()
        log_likelihood_given_tree(
          tree_fmat = tree_fmat,
          coal_times = coal_times,
          sequences = sequences,
          mode = log_likelihood_mode,
          R = log_likelihood_R
        )$log_likelihood
      },
      numeric(1)
    )
    ll_pb$terminate()
  }

  chain_coords <- coords[coords$group == "M Chain", , drop = FALSE]
  matrix_coords <- coords[coords$label != "", , drop = FALSE]

  if (color_by_log_likelihood) {
    matrix_coords$matrix_group <- factor(matrix_coords$label, levels = labels)
    matrix_labels_with_ll <- paste0(
      labels,
      " (LL=",
      formatC(coords$log_likelihood[coords$label != ""], format = "f", digits = 2),
      ")"
    )

    p <- ggplot()
    if (nrow(chain_coords) > 0) {
      p <- p +
        geom_point(
          data = chain_coords,
          aes(x = Dim1, y = Dim2, color = log_likelihood),
          size = 1.5,
          alpha = 0.5
        )
    }
    p <- p +
      geom_point(
        data = matrix_coords,
        aes(x = Dim1, y = Dim2, color = log_likelihood, shape = matrix_group),
        size = 4.5
      ) +
      geom_text(
        data = matrix_coords,
        aes(x = Dim1, y = Dim2, label = label),
        color = "black",
        vjust = -0.8,
        show.legend = FALSE
      ) +
      scale_color_gradientn(colors = c("navy", "skyblue", "gold", "firebrick")) +
      scale_shape_manual(values = rep(8, length(labels)), labels = matrix_labels_with_ll, name = "Matrix") +
      guides(shape = guide_legend(override.aes = list(color = "black", size = 4.5))) +
      labs(title = title, x = "MDS 1", y = "MDS 2", color = "Log Likelihood") +
      theme_minimal(base_size = 14)
  } else {
    p <- ggplot(coords, aes(x = Dim1, y = Dim2, color = group))
    if (nrow(chain_coords) > 0) {
      p <- p + geom_point(data = chain_coords, size = 1.5, alpha = 0.5)
    }
    p <- p +
      geom_point(data = matrix_coords, size = 3) +
      geom_text(data = matrix_coords, aes(label = label), vjust = -0.8, show.legend = FALSE) +
      scale_color_manual(values = colors[seq_len(length(group_levels))]) +
      labs(title = title, x = "MDS 1", y = "MDS 2", color = "Matrix") +
      theme_minimal(base_size = 14)
  }

  print(p)

  if (!is.null(save_path)) {
    ggsave(filename = paste0("plots/", save_path), plot = p, width = 8, height = 6, dpi = 150)
  }
}


tree_histogram_comparison_plot <- function(chain_gen, chain_data, M_true_fmat, burn_in = 0,
                                           title = "Distance Comparison", save_path = NULL) {
  # Compute distances for data chain
  distances_data <- sapply(chain_data[(burn_in+1):length(chain_data)], function(tree) {
    distance_Fmat(tree, M_true_fmat, dist = "l2")^2
  })

  # Compute distances for generated chain
  distances_gen <- sapply(chain_gen[(burn_in+1):length(chain_gen)], function(tree) {
    distance_Fmat(tree, M_true_fmat, dist = "l2")^2
  })

  # Combine into one data frame for ggplot
  df <- data.frame(
    Distance = c(distances_data, distances_gen),
    Group = factor(c(rep("Data", length(distances_data)), rep("Generated", length(distances_gen))),
                   levels = c("Data", "Generated"))
  )

  # Plot
  p <- ggplot(df, aes(x = Distance, fill = Group, y = after_stat(density))) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 30, color = "black") +
    scale_fill_manual(values = c("steelblue", "tomato")) +
    labs(title = title, x = "Squared L2 Distance", y = "Density", fill = "Group") +
    theme_minimal(base_size = 14)

  print(p)

  # Save if path provided
  if (!is.null(save_path)) {
    ggsave(filename = paste0("plots/", save_path), plot = p, width = 8, height = 6, dpi = 150)
  }
}


matrix_chain_ess <- function(matrix_chain, normalize = FALSE, per_n_samples = 1000) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for ESS calculations. Please install it first.")
  }

  chain_samples <- Filter(Negate(is.null), matrix_chain)

  if (length(chain_samples) == 0) {
    stop("matrix_chain must contain at least one non-NULL matrix.")
  }

  first_dim <- dim(chain_samples[[1]])
  if (length(first_dim) != 2) {
    stop("Each element of matrix_chain must be a matrix.")
  }

  same_dim <- vapply(
    chain_samples,
    function(mat) is.matrix(mat) && identical(dim(mat), first_dim),
    logical(1)
  )
  if (!all(same_dim)) {
    stop("All matrices in matrix_chain must have the same dimensions.")
  }

  if (!is.numeric(per_n_samples) || length(per_n_samples) != 1 || per_n_samples <= 0) {
    stop("per_n_samples must be a single positive number.")
  }

  n_iter <- length(chain_samples)
  chain_matrix <- t(vapply(chain_samples, c, numeric(prod(first_dim))))
  ess_values <- vapply(
    seq_len(ncol(chain_matrix)),
    function(idx) as.numeric(coda::effectiveSize(chain_matrix[, idx])),
    numeric(1)
  )
  ess_matrix <- matrix(ess_values, nrow = first_dim[1], ncol = first_dim[2])
  dimnames(ess_matrix) <- dimnames(chain_samples[[1]])

  if (normalize) {
    ess_matrix <- ess_matrix / n_iter * per_n_samples
  }

  ess_matrix
}


matrix_chain_ess_heatmap <- function(matrix_chain, normalize = FALSE, per_n_samples = 1000,
                                     title = NULL, save_path = NULL) {
  ess_matrix <- matrix_chain_ess(
    matrix_chain = matrix_chain,
    normalize = normalize,
    per_n_samples = per_n_samples
  )

  ess_df <- as.data.frame(as.table(ess_matrix), responseName = "ESS")
  colnames(ess_df) <- c("Row", "Column", "ESS")

  if (is.null(title)) {
    if (normalize) {
      title <- paste0("ESS Heatmap (per ", per_n_samples, " samples)")
    } else {
      title <- "ESS Heatmap"
    }
  }

  fill_label <- if (normalize) {
    paste0("ESS per ", per_n_samples)
  } else {
    "ESS"
  }

  fill_scale <- if (normalize) {
    scale_fill_gradient(
      low = "grey95",
      high = "darkgreen",
      trans = "log10",
      limits = c(1, 1000),
      oob = scales::squish
    )
  } else {
    scale_fill_gradient(
      low = "grey95",
      high = "darkgreen",
    )
  }

  p <- ggplot(ess_df, aes(x = Column, y = Row, fill = ESS)) +
    geom_tile(color = "white", linewidth = 0.3) +
    fill_scale +
    labs(title = title, x = "Column", y = "Row", fill = fill_label) +
    theme_minimal(base_size = 14) +
    scale_y_discrete(limits = rev(levels(ess_df$Row))) +
    coord_fixed()

  print(p)

  if (!is.null(save_path)) {
    ggsave(filename = paste0("plots/", save_path), plot = p, width = 8, height = 6, dpi = 150)
  }

  invisible(list(plot = p, ess_matrix = ess_matrix))
}



#library("TreeTools")
update_time <- function(tree, coal_times) {
  # coal_times <- cumsum(coalescent.intervals(tree)$interval.length)
  # tiplabels(); nodelabels()
  ## sort before update to avoid problems
  # xx1 <- sort(n.t, index.return = T)
  
  old.edge <- tree$edge
  n.sample <- tree$Nnode + 1
  t.tot <- max(ape::node.depth.edgelength(tree))
  n.t <- t.tot - ape::node.depth.edgelength(tree)  ## gives the node length
  n.t[1:n.sample] <- 0
  new.n.t <- n.t
  
  # order nodes according to length, then coalescent times are in reverse order
  xx1 <- sort(new.n.t, index.return = TRUE)
  index <- (2 * n.sample - 1):(n.sample + 1)
  
  for (j in (n.sample + 1):(2 * n.sample - 1)) {
    old.edge[which(tree$edge[,1] == xx1$ix[j]), 1] <- index[j - n.sample]
    old.edge[which(tree$edge[,2] == xx1$ix[j]), 2] <- index[j - n.sample]
  }
  
  new.n.t[(n.sample + 1):(2 * n.sample - 1)] <- rev(coal_times)
  
  # If we sort them, we can get the correspondence between leaves and coal. times
  xx <- sort(new.n.t, index.return = TRUE)
  new.edge <- old.edge
  
  for (j in (n.sample + 1):(2 * n.sample - 1)) {
    new.edge[which(old.edge[,1] == xx$ix[j]), 1] <- index[j - n.sample]
    new.edge[which(old.edge[,2] == xx$ix[j]), 2] <- index[j - n.sample]
  }
  
  # check <- new.edge[new.edge[,2] > n,]
  # check2 <- check[,1] - check[,2]
  # while (sum(check2[check2 > 0]) > 0) {
  #   maxcon <- which.max(check[,1] - check[,2])
  #   changeto <- check[maxcon,2]
  #   changefrom <- check[maxcon,1]
  #   new.edge2 <- new.edge
  #   new.edge2[new.edge[,1] == changefrom, 1] <- changeto
  #   new.edge2[new.edge[,1] == changeto, 1] <- changefrom
  #   new.edge2[new.edge[,2] == changefrom, 2] <- changeto
  #   new.edge2[new.edge[,2] == changeto, 2] <- changefrom
  #   new.edge <- new.edge2 
  #   check <- new.edge[new.edge[,2] > n,]
  #   check2 <- check[,1] - check[,2]
  # }
  
  new.edge.length <- new.n.t[new.edge[,1]] - new.n.t[new.edge[,2]]
  new.tree <- tree
  new.tree$edge <- new.edge
  new.tree$edge.length <- new.edge.length
  # new.tree$tip.label <- rev(tree$tip.label)
  
  tree2 <- write.tree(new.tree)
  new.tree2 <- read.tree(text = tree2)
  
  trees <- c(tree, new.tree2) 
  trees <- .compressTipLabel(trees)
  t2 <- trees[[2]]
  
  return(t2)
}

  
