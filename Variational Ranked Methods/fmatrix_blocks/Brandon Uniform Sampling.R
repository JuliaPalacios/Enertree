uniform_f_matrix <- function(n_tips) {
  m <- 2*n_tips - 2
  fmat <- matrix(0, nrow = m, ncol = m)
  
  fixed_diag <- c(2:n_tips, (n_tips - 1):1)
  
  for (i in 1:m) {
    fmat[i, i] <- fixed_diag[i]
  }
  
  for (i in 2:m) {
    D <- numeric(i - 1)
    
    D[1] <- fmat[i - 1, 1]
    
    if (i > 2) {
      for (j in 2:(i - 1)) {
        D[j] <- fmat[i - 1, j] - fmat[i - 1, j - 1]
      }
    }
    
    L <- base::sample(seq_along(D), size = 1, prob = D)
    
    for (j in 1:(i - 1)) {
      if (j < L) {
        fmat[i, j] <- fmat[i - 1, j]
      } else {
        fmat[i, j] <- fmat[i - 1, j] - 1
      }
    }
  }
  
  fmat
}

enumerate_labeled_isochronous_fmats <- function(diagonal) {
  m <- length(diagonal)
  F_start <- matrix(0L, nrow = m, ncol = m)
  diag(F_start) <- diagonal
  
  states <- list(list(F = F_start, probability = 1))
  
  for (i in 2:m) {
  next_states <- list()
    for (state in states){
      fmat <- state$F
      prob = state$probability
      cols <- seq_len(i - 1L)
      previous_row <- fmat[i - 1L, cols]
      D <- previous_d_counts(previous_row)
      for (L in which(D > 0L)){
        new_fmat <- fmat
        if(D[L] != 0){
          new_fmat[i, cols] <- previous_row - as.integer(cols >= L)
          new_prob<- prob *  D[L]/sum(D)
          next_states[[length(next_states) + 1L]] <- list(
            F = new_fmat,
            probability = new_prob
          )
        }
      }
    }
  states <- next_states
  
  }
  
  keys <- vapply(states, function(state) matrix_key(state$F), character(1))
  if (any(duplicated(keys))) {
    stop("Enumeration produced duplicate F-matrices.")
  }
  
  states
  
  
  
  
}