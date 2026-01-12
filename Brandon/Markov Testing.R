library(fmatrix)
library(progress)

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


sampleF<-function(M,beta,iter,diam, startF = NULL){
  acceptance_rate<-0
  n<-nrow(M)+1
  
  chainF<-list()
  if (is.null(startF)){
    current_tree<-rcoal(n)
    chainF[[1]]<-gen_Fmat(current_tree,tol=13)
    
  }else{
    chainF[[1]]<- startF
  }
  
  

  
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
    #print(prob)
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


n = 25
sample_tree <- rcoal(n)
beta = 1
diam = 1
iter = 100000

M <- gen_Fmat(sample_tree)

samples <-sampleF(M,beta,iter,diam)

chain <-samples$chainF
chain <- chain[500:length(chain)]

#mean_tree <- sa_mean(chain)
mean_tree <- M

iter = iter -500
l2_dist <- rep(0, iter)

pb = progress_bar$new("[:bar] :percent", total=iter)
for (i in 1:iter){
  pb$tick()
  l2_dist[i] <- distance_Fmat(chain[[i]], mean_tree, dist = "l2")
  
}

l2_dist^2 #<- l2_dist


beta_test <- 1 / (2 * mean(l2_dist^2))
l0 = unique(l2_dist^2)
y<-exp(-beta_test*sort(l0))
y2 <-exp(-beta*sort(l0))

# histogram as density
x1 <- sort(l0)
p1 <- y  / sum(y)
p2 <- y2 / sum(y2)

# draw histogram
hist(l2_dist^2,
     probability = TRUE,
     breaks = 50,
     
     xlab = "L2 distance squared")

# add curves
lines(x1, p1, lwd = 2, col = "red")
lines(x1, p2, lwd = 2, col = "purple")







centre <- M
diam <- 1

k <- 10
trials <- 10000
batch_size <- 50
print_every <- 1000

beta_star <- 6
step0 <- 0.05

n_data <- length(chain)

for (t in 1:trials) {
  
  idx <- sample.int(n_data, batch_size, replace = TRUE)
  
  s_pos_sum <- 0
  s_neg_sum <- 0
  acc_sum   <- 0
  
  for (b in 1:batch_size) {
    
    x0 <- chain[[ idx[b] ]]  
    
    d0 <- distance_Fmat(x0, centre, dist = "l2")
    s_pos <- d0^2
    
    out <- sampleF(M = centre,
                   beta = beta_star,
                   iter = k + 1,
                   diam = diam,
                   startF = x0)
    
    xk <- out$chainF[[ k + 1 ]]
    dk <- distance_Fmat(xk, centre, dist = "l2")
    s_neg <- dk^2
    
    s_pos_sum <- s_pos_sum + s_pos
    s_neg_sum <- s_neg_sum + s_neg
    acc_sum   <- acc_sum + out$acceptance_rate
  }
  
  E_pos <- s_pos_sum / batch_size
  E_neg <- s_neg_sum / batch_size
  acc   <- acc_sum / batch_size
  
  g_beta <- (E_neg - E_pos) / diam
  
  step <- step0 #/ sqrt(t)              
  beta_star <- beta_star + step * g_beta
  if (beta_star < 1e-8) beta_star <- 1e-8
  
  if (t %% print_every == 0) {
    cat("iter", t,
        "beta", beta_star,
        "E_pos", E_pos,
        "E_neg", E_neg,
        "acc", acc, "\n")
  }
}

beta_hat <- beta_star
beta_hat

test_tree <- rcoal(6)
plot(test_tree)

my_encod_test <- my_encod(gen_Fmat(test_tree,tol=13))




beta_test <- 1 / (2 * mean(l2_dist^2))
l0 = unique(l2_dist^2)
y<-exp(-beta_hat*sort(l0))
y2 <-exp(-beta*sort(l0))

# histogram as density
x1 <- sort(l0)
p1 <- y  / sum(y)
p2 <- y2 / sum(y2)

# draw histogram
hist(l2_dist^2,
     probability = TRUE,
     breaks = 50,
     xlab = "L2 distance squared")

# add curves
lines(x1, p1, lwd = 2, col = "red")
lines(x1, p2, lwd = 2, col = "purple")
