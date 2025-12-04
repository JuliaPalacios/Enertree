library(fmatrix)

Logdistance_prob <- function(fmat, centre=NULL, beta=1,diam=1){
  ## upto proportionality constant
  n <- dim(fmat)[1] + 1
  
  if(is.null(centre)){
    centre <- kingman_m(n)
  }
  
  return(-(beta/diam) * distance_Fmat(fmat, centre,dist="l2")^2)
}

distance_prob <- function(fmat, centre=NULL, beta=1,diam=1){
  ## upto proportionality constant
  n <- dim(fmat)[1] + 1
  
  if(is.null(centre)){
    centre <- kingman_m(n)
  }
  
  return(exp(-(beta/diam) * distance_Fmat(fmat, centre,dist="l2")^2))
}
##Generates samples of F-matrices via Metropolist-Hastings
##from the Gibbs distribution
##Samyak's chain
sampleF<-function(M,beta,iter,diam){
  print(M)
  acceptance_rate<-0
  n<-nrow(M)+1
  current_tree<-rcoal(n)
  chainF<-list()
  chainF[[1]]<-gen_Fmat(current_tree,tol=13)
  MeanM<-matrix(0,nrow=n-1,ncol=n-1)
  current_Encod<-my_encod(chainF[[1]])
  f_current <- Fmat_from_myencod(current_Encod)
  Logdist_curr<-Logdistance_prob(f_current,M,beta=beta,diam=diam)
  MeanDist<-0
  MeanEDist<-0
  for (j in 1:iter){
    proposed <- proposal_myencod(current_Encod)
    f_proposed <- Fmat_from_myencod(proposed)
    Logdist_prop<-Logdistance_prob(f_proposed,M,beta=beta,diam=diam)
    prob <- exp(Logdist_prop-Logdist_curr)
    if (runif(1) < prob){
      chainF[[j]] <- f_proposed
      current_Encod<-proposed
      f_current<-f_proposed
      Logdist_curr<-Logdist_prop
      acceptance_rate<-acceptance_rate+1
    } else{
      chainF[[j]]<-f_current
    }
    MeanDist<-MeanDist-Logdist_curr/beta
    MeanEDist<-MeanEDist+exp(Logdist_curr)
    
    MeanM<-MeanM+chainF[[j]]
  }
  MeanM<-MeanM/iter
  MeanDist<-MeanDist/iter
  MeanEDist<-MeanEDist/iter
  return(list(chainF=chainF,acceptance_rate=acceptance_rate/iter,MeanM=MeanM,MeanDist=MeanDist,MeanEDist=MeanEDist))
}

normalizing<-function(beta,maxv,grid.size){
  #min and max of distance squared
  #grid.size number of grid points
  delta<-maxv/(grid.size-1)
  grid<-delta/2+seq(0,maxv,by=delta)
  return(sum(delta*exp(-beta*grid)))
}



#Returns the cardinality of the space of ranked tree shapes
ZigZag<-function(n){
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


labels_prob<-function(F){
  nr<-nrow(F)
  c<-sum(diff(F[nr,])==2)
  return(log(2^(c)/factorial(nr+1)))
}


gen_caterpillar <- function(n){
  
  "Return caterpillar tree"
  F_mat <- matrix(rep(seq(1,n-1),n-1),nrow=n-1,byrow = T)
  F_mat[upper.tri(F_mat)] <- 0
  diag(F_mat) <- seq(2,n)
  return(F_mat)
}

nearby_Fmat_sampling <- function(fmat) {
  samples <- generate_sample_markov_chain(M=fmat, b=5, num_samps=10,burnin=1000)
  dist_from_fmat <- sapply(samples, function(sample) norm(sample-fmat, type="F"))
  closest_samples <- samples[dist_from_fmat == min(dist_from_fmat)]
  norms <- sapply(closest_samples, function(sample) norm(sample, type="F"))
  smallest_norm <- closest_samples[norms == min(norms)]
  return(closest_samples[[1]])
}
