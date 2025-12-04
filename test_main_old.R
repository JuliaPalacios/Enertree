#Will focus only on topology having everything else the truth

source("~/Documents/Phylo-VB/Phylo-VB/all_functions.R")
#####################R Ne, genealogy simulation
library(phylodyn)
library(fmatrix)
bottleneck_traj<-function(t){
  result=rep(0,length(t))
  result[t<=0.5]<-1
  result[t>0.5 & t<1]<-.1
  result[t>=1]<-1
  return(result)
}

exp_traj<-function(t){
  result=rep(0,length(t))
  result[t<=0.1]<-10
  result[t>0.1 & t<0.25]<-10*exp(2-20*t[t>0.1 & t<0.25]) #50 times smaller
  result[t>=0.25]<-0.5
  return(result)
}
set.seed(123)
samp_times<-0
n_sampled<-n<-15


 
#ZigZag(n)
##We are thinking of isochronous trees first
library("ape")
simul1<-coalsim(samp_times = samp_times, n_sampled = n_sampled, traj = exp_traj,method="tt",val_upper=11)

#simul1<-coalsim(samp_times = samp_times, n_sampled = n_sampled, traj = bottleneck_traj,method="tt",val_upper=11)
tree<-generate_newick((simul1))
plot(ladderize(tree$newick),show.tip.label = T)
axisPhylo()
trueM<-gen_Fmat(tree$newick,tol=8)
est<-BNPR(tree$newick)
plot_BNPR(est)
points(seq(0,2.5,by=.1),exp_traj(seq(0,2.5,by=.1)),type="l",col="blue")
###################Given the tree, simulate data from JC
library("phyclust",quiet=TRUE)
tree.newick<-write.tree(tree$newick)
seqgen(opts="-mHKY -t0.5 -f0.25,0.25,0.25,0.25 -l100",newick.tree=tree.newick,temp.file="temp.fasta")
#install.packages("phylotools")
library("phylotools")
data2<-read.phylip("temp.fasta")
data<-dat2fasta(data2,outfile="out.fasta")
library("phangorn")
mydna<-as.phyDat(read.FASTA("out.fasta"))
#tree length
#treedata<-coalescent.intervals(tree$newick)
#tree_length<-sum(treedata$interval.length*treedata$lineages)
#tree_height<-treedata$total.depth

#############Given the data, start with the UPGMA, I am ignoring mutation rate atm
dm <- as.matrix(dist.ml(mydna,"JC69"))
upgma<-upgma(dm)

#compare truth with rough estimate of tree
par(mfrow=c(1,2))
plot(ladderize(tree$newick),show.tip.label = T, main="Truth")
axisPhylo()

plot(ladderize(upgma),main="UPGMA")
axisPhylo()
inter_coal_times<-coalescent.intervals(upgma)$interval.length
inter_coal_times[which(inter_coal_times<=0.001)]<-.01
coal_times<-cumsum(inter_coal_times)
upgma<-update_time(upgma,coal_times)


##We will need to solve multifurcations with UPGMA

##If we endup coding the heterochronous likelihood, look at Jonathan's code. He considers missing and other types of missingness


#let's compute the log-likelihood of the tree
#fit=pml(upgma,data=mydna,model="JC69")
#fit=pml(tree,data=mydna,model="JC69")
#let's find the branches that optimize the likelihood
# fitJC=optim.pml(fit,optRooted=TRUE,optNni=FALSE,optBf=FALSE,rearrangement=FALSE)
# logLik(fitJC)
# plot(ladderize(fitJC$tree))
# axisPhylo()

###Estimate Ne with INLA
Ne<-BNPR(upgma)
plot_BNPR(Ne)
points(seq(0,2.5,by=.1),exp_traj(seq(0,2.5,by=.1)),type="l",col="blue")

#Coalescent log prior

#install.packages("TreeTools")
library("TreeTools")
##This function needs to be validated, updates the branch lengths of the tree
update_time<-function(tree,coal_times){
  #coal_times<-cumsum(coalescent.intervals(tree)$interval.length)
  #tiplabels();nodelabels()
  ##sort before update to avoid problems
  #xx1<-sort(n.t,index.return=T)
  old.edge<-tree$edge
  n.sample <- tree$Nnode + 1
  t.tot <- max(ape::node.depth.edgelength(tree))
  n.t <- t.tot - ape::node.depth.edgelength(tree) ##gives the node length
  new.n.t<-n.t
  #order nodes according to length, then coalescent times are in reverse order
  xx1<-sort(new.n.t,index.return=T)
  index<-(2*n.sample-1):(n.sample+1)
  for (j in (n.sample+1):(2*n.sample-1)){
    old.edge[which(tree$edge[,1]==xx1$ix[j]),1]<-index[j-n.sample]
    old.edge[which(tree$edge[,2]==xx1$ix[j]),2]<-index[j-n.sample]
  }
  new.n.t[(n.sample+1):(2*n.sample-1)]<-rev(coal_times)
  #If we sort them, we can get the correspondence between leaves and coal. times
  xx<-sort(new.n.t,index.return=T)
  new.edge<-old.edge
  for (j in (n.sample+1):(2*n.sample-1)){
    new.edge[which(old.edge[,1]==xx$ix[j]),1]<-index[j-n.sample]
    new.edge[which(old.edge[,2]==xx$ix[j]),2]<-index[j-n.sample]
  }
  new.edge.length<-new.n.t[new.edge[,1]]-new.n.t[new.edge[,2]]
  new.tree<-tree
  new.tree$edge<-new.edge
  new.tree$edge.length<-new.edge.length
  #new.tree$tip.label<-rev(tree$tip.label)
  tree2<-write.tree(new.tree)
  new.tree2<-read.tree(text=tree2)
  trees <- c(tree, new.tree2) 
  trees <- .compressTipLabel(trees)
  t2 <- trees[[2]]
  return(t2)
}


  
##Check I am inputing the right grid
#coal_init<-phylodyn:::coal_lik_init(samp_times=samp_times, n_sampled=n_sampled, coal_times=summary$coal_times, grid=Ne$grid)
#coal_lik<-phylodyn:::coal_loglik(init=coal_init, f=log(Ne$effpop))



    
##Start with M matrix corresponding to Kingman coalescent
##https://github.com/RSamyak/fmatrix
#install_github("RSamyak/fmatrix")

##Sample uniform ranked tree shapes and compute the expectation of e^beta*d^2
#this function should be replaced with whatever Leon thinks is best, atm uses a 
#Markov chain
sampleU<-function(M,beta,iter,diam){
  n<-nrow(M)+1
  current_tree<-rcoal(n)
  chainF<-list()
  chainF[[1]]<-gen_Fmat(current_tree,tol=13)
  MeanEDist<-0
  current_Encod<-my_encod(chainF[[1]])
  for (j in 1:iter){
    proposed <- proposal_myencod(current_Encod)
    chainF[[j]] <- Fmat_from_myencod(proposed)
    current_Encod<-proposed
    dist_prop<-distance_prob(chainF[[j]],M,beta,diam=diam)
    MeanEDist<-MeanEDist+dist_prop
    }
  MeanEDist<-MeanEDist/iter
  return(list(MeanEDist=MeanEDist))
}


# acceptance_prob <- function(proposed, current){
#   f_proposed <- Fmat_from_myencod(proposed)
#   f_current <- Fmat_from_myencod(current)
#   return(list(accep=distance_prob(f_proposed)/distance_prob(f_current),old=f_current,new=f_proposed))
# }

##using the diameter with the hope that there should be some scaling
#proportional to entropy

diameter<-function(n){
  return(n*(n-2)*(n^2-4*n+6)/48 + ((-1)^n - 1)*(2*n-3)/32)
}

diam<-diameter(n)
diam<-1
##Start VB on topologies, Algorithm 1 with two modifications
##Modification 1: the gradient of the variational distributions is approximated via MCMC
##Modification 2: I am ignoring the normalizing constant in the last expression of ELBO
#mytree_from_F
coal<-simul1$coal_times  #For now, we assume we know the true times
#M<-kingman_m(n_sampled)
S<-30  #number of iterations for gradient estimation per iteration
M<-gen_Fmat(upgma,tol=8) #start with upgma ranked tree shape
Mupgma<-gen_Fmat(upgma,tol=8)
beta<-runif(1)
#coal<-expected_time
#coal<-coal_times
iter<-1000
elbo_tau<-rep(0,iter)
#Traditional
rho <-1/seq(2,iter+1)
#Adagrad (we initialize the gradient somehow)
ada <- TRUE
eta <- 1/50
chain<- sampleF(M,beta,100,diam) #after thinning 100 samples, start the chain according to Kingman 
chain2<-sampleU(M,beta,100,diam)
z<-chain$chainF[[100]]
c<-sum(diff(z[n-1,])==2) #number of cherries of current tree ranked tree shape
z.tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),z)),coal)
##need to check this function in detail and see how the tip assignment is done
##ignoring the normalizing constant
distSq<-Logdistance_prob(z,M,beta,diam)
z.tree$tip.label<-sample(new.tree$tip.label) #here we permute the labels
conditional<-pml(z.tree,data=mydna,model="JC69")$logLik+log(2^(n-1-c)/factorial(n-1))-distSq+log(chain2$MeanEDist)+log(ZigZag(n))
grad_beta<-(chain$MeanDist+distSq/beta)*conditional
grad_M<-(2*beta*(z-chain$MeanM))*conditional
adagrad_beta <- grad_beta^2
adagrad_M<-grad_M^2

##Draw s samples from q
for (i in 1:iter){
  chain<-sampleF(M,beta,100*S,diam) #inter-coal times
  grad_beta<-0
  grad_M<-0
  for (j in 1:S){
    z<-chain$chainF[[100*j]]
    chain2<-sampleU(M,beta,100,diam)
    c<-sum(diff(z[n-1,])==2)
    z.tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),z)),coal)
    distSq<-Logdistance_prob(z,M,beta,diam)
    z.tree$tip.label<-sample(new.tree$tip.label) #here we permute the labels
    conditional<-pml(z.tree,data=mydna,model="JC69")$logLik+log(2^(n-1-c)/factorial(n-1))-distSq+log(chain2$MeanEDist)+log(ZigZag(n))
    grad_beta<-grad_beta+(chain$MeanDist+distSq/beta)*conditional
    grad_M<-grad_M+(2*(beta/diam)*(z-chain$MeanM))*conditional
    elbo_tau[i]<-elbo_tau[i]+conditional
  }
  grad_beta<-grad_beta/S
  grad_M<-grad_M/S
  elbo_tau[i]<-elbo_tau[i]/S
  print("elbo")
  print(elbo_tau[i])
  #print(paste("Average grad",grad,sep=""))
  if (ada==FALSE){
    beta<-beta+rho[i]*grad_beta
    M<-M+rho[i]*grad_M
  } else {
    adagrad_beta <- adagrad_beta+grad_beta^2
    beta <- beta+eta*adagrad_beta^(-1/2)*grad_beta
    adagrad_M<-adagrad_M+grad_M^2
    temp<-adagrad_M^(-1/2)
    temp[upper.tri(temp,T)]<-0
    temp[cbind(2:nrow(M),1:(nrow(M)-1))]<-0
    temp2<-eta*temp*grad_M
    temp2[grad_M==0]<-0
    M<-M+temp2
    #M[M<0]<-0
  }
  if (max(M)>n){break}
  print(paste("Iter",i,sep=""))
  print(beta)
  print(sqrt(sum((M-trueM)^2)))
}
M[M<0]<-0

#assume you found M and beta, then what? How do you find the best topology?
#sample from variational distribution and only accept those that increase the unnormalized posterior

beta2<-beta*2
#Sample from VB
iter<-1
old_like<-pml(z.tree,data=mydna,model="JC69")$logLik+log(2^(n-1-c)/factorial(n-1))
accp<-0
for (i in 1:iter){
  chain<-sampleF(M,beta2,5000,diam) #inter-coal times
  for (j in 1:5000){
    z<-chain$chainF[[j]]
    c<-sum(diff(z[n-1,])==2)
    z.tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),z)),coal)
    for (k in 1:1000){
      z.tree$tip.label<-sample(new.tree$tip.label) #here we permute the labels
      trees <- c(z.tree, tree$newick) 
      trees <- .compressTipLabel(trees)
      t2 <- trees[[1]]
      new_like<-pml(z.tree,data=mydna,model="JC69")$logLik+log(2^(n-1-c)/factorial(n-1))
      if (new_like>old_like){
        current<-z.tree
        accp<-accp+1
        old_like<-new_like
        print(i)
        print(new_like)
      }
    }
  }
}



plot(ladderize(current))
tiplabels()
axisPhylo()
plot(ladderize(tree$newick),show.tip.label = T)
tiplabels()
axisPhylo()

pml(current,data=mydna,model="JC69")$logLik
pml(tree$newick,data=mydna,model="JC69")$logLik
pml(upgma,data=mydna,model="JC69")$logLik


##Now sample from the variational distribution and compare distributions on
#ranked tree shapes with BEAST and the tree distance distribution

trees_variational<-sampleF(M,beta,10001,diam=1)

trees_beast<-read.nexus("~/Documents/Beast/out.trees")
#Distance to true M
dist_beast<-rep(0,10001)
dist_var<-rep(0,10001)
for (j in 1:length(trees_beast)){
  Fmat<-gen_Fmat(trees_beast[[j]],tol=8)
  dist_beast[j]<-distance_Fmat(Fmat, trueM,dist="l2")
  dist_var[j]<-distance_Fmat(trees_variational$chainF[[j]],trueM,dist="l2")
}

p1<-hist(dist_beast,main="MCMC")
p2<-hist(dist_var)

plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,3),main="MCMC vs VI",xlab="Distance to true M")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,3), add=T)  # second

##Now for each Fmatrix in population, sample 100 labels and average the posterior
#then compute distance to true and sum

####Deprecated
# beta2<-beta*10
# #Sample from VB
# iter<-1
# old_like<-pml(z.tree,data=mydna,model="JC69")$logLik+log(2^(n-1-c)/factorial(n-1))
# accp<-0
# for (i in 1:iter){
#     chain<-sampleF(M,beta2,1000) #inter-coal times
#     for (j in 1:1000){
#     z<-chain$chainF[[j]]
#     c<-sum(diff(z[n-1,])==2)
#     z.tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),z)),coal)
#     for (k in 1:10){
#     z.tree$tip.label<-sample(new.tree$tip.label) #here we permute the labels
#     trees <- c(z.tree, tree$newick) 
#     trees <- .compressTipLabel(trees)
#     t2 <- trees[[1]]
#     new_like<-pml(z.tree,data=mydna,model="JC69")$logLik+log(2^(n-1-c)/factorial(n-1))
#     if (new_like>old_like){
#       current<-z.tree
#       accp<-accp+1
#       old_like<-new_like
#       print(i)
#       print(new_like)
#     }
#     }
#     }
# }

#Found M, then fix and only sample labels
iter<-100
old_like<-pml(current,data=mydna,model="JC69")$logLik+log(2^(n-1-c)/factorial(n-1))
accp<-0
for (i in 1:iter){
  chain<-sampleF(M,beta,1000) #inter-coal times
  for (j in 1:1000){
    z<-chain$chainF[[j]]
    c<-sum(diff(z[n-1,])==2)
    z.tree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),z)),coal)
    z.tree$tip.label<-sample(new.tree$tip.label) #here we permute the labels
    trees <- c(z.tree, tree$newick) 
    trees <- .compressTipLabel(trees)
    t2 <- trees[[1]]
    new_like<-pml(z.tree,data=mydna,model="JC69")$logLik+log(2^(n-1-c)/factorial(n-1))
    if (new_like>old_like){
      current<-z.tree
      accp<-accp+1
      old_like<-new_like
      print(i)
      print(new_like)
    }
  }
}

##Comparing normalizing constants
data("F-lists")
beta_seq<-seq(0.1,2,by=0.01)
beta_true<-rep(0,length(beta_seq))
beta_approx<-rep(0,length(beta_seq))
beta_approx2<-rep(0,length(beta_seq))

for (j in 1:length(beta_seq)){
  Flist_probs <- sapply(F.list7, distance_prob,centre=F.list7[[length(F.list7)]],beta=beta_seq[j])
  variational<-Flist_probs/sum(Flist_probs)
  beta_true[j]<-sum(Flist_probs)
  beta_approx[j]<-ZigZag(7)*sum(Flist_probs)
  #beta_approx[j]<-ZigZag(7)*sum(Flist_probs*variational)
  beta_approx2[j]<-ZigZag(7)*(exp(sum(log(Flist_probs)*variational)))
  
}

plot(beta_seq,beta_true,type="l",main="true Z vs approx Z, n=7")
points(beta_seq,beta_approx,type="l",col="red")
points(beta_seq,beta_approx2,type="l",col="blue")
text(3,20,"approx Z",col="red")
text(3,30,"Geo approx Z",col="blue")
abline(v=1)
abline(v=0.5)


data("F-lists")
beta_seq<-seq(0.1,5,by=0.05)
beta_true<-rep(0,length(beta_seq))
beta_approx<-rep(0,length(beta_seq))

beta_approx2<-rep(0,length(beta_seq))
for (j in 1:length(beta_seq)){
Flist_probs <- sapply(F.list10, distance_prob,centre=F.list10[[length(F.list10)]],beta=beta_seq[j])
variational<-Flist_probs/sum(Flist_probs)
beta_true[j]<-sum(Flist_probs)
beta_approx[j]<-ZigZag(10)*sum(Flist_probs*variational)
beta_approx2[j]<-ZigZag(10)*(exp(sum(log(Flist_probs)*variational)))
}

plot(beta_seq,beta_true,type="l",main="true Z vs approx Z, n=10")
points(beta_seq,beta_approx,type="l",col="red")
points(beta_seq,beta_approx2,type="l",col="blue")
text(3,20,"approx Z",col="red")

###############################
Flist_probs <- sapply(F.list7, distance_prob,centre=M,beta=beta)
variational<-Flist_probs/sum(Flist_probs)
plot(sort(variational),type="l")

Flist_probs <- sapply(F.list10, distance_prob,centre=F.list10[[length(F.list10)]],beta=0.5)
variational<-Flist_probs/sum(Flist_probs)
plot(sort(variational),type="l")

Flist_probs <- sapply(F.list10, distance_prob,centre=F.list10[[length(F.list10)]],beta=0.05)
variational<-Flist_probs/sum(Flist_probs)
plot(variational,main="beta .05")
plot(sort(variational),type="l")


Flist_probs <- sapply(F.list10, distance_prob,centre=F.list10[[length(F.list10)]],beta=0.5)
variational<-Flist_probs/sum(Flist_probs)
plot(variational,type="l")
plot(sort(variational),type="l")

posterior<-rep(0,length(F.list5))
for (j in 1:length(F.list5)){
  for (i in 1:1000){
  c<-sum(diff(F.list5[[j]][nrow(F.list5[[1]]),])==2)
  mytree<-tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),F.list5[[j]])),coal)
  mytree$tip.label<-sample(tree$newick$tip.label)
  trees <- c(mytree, tree$newick) 
  trees <- .compressTipLabel(trees)
  t2 <- trees[[1]]
  posterior[j]<-posterior[j]+pml(t2,data=mydna,model="JC69")$logLik+log(2^(n-c-1)/factorial(n-1))
  }
  posterior[j]<-posterior[j]/1000
}
barplot(variational,main="Variational")

barplot(exp(posterior)/sum(exp(posterior)),main="posterior")

Flist_probs <- sapply(F.list7, distance_prob,centre=kingman_m(7),beta=.4)
variational<-Flist_probs/sum(Flist_probs)


# 
# n <- 7
# chain <- matrix(0,nrow=MAX,ncol=n-2) #each row is a realization, Samyak's encoding
# 
# set.seed(1)
# current <- init <- chain[[1]] <- rEncod(1, n=n)[[1]]
# counter <- 1
# ##i would rather start from a random coalescent
# MAX <- 1000
#   
# 
#  proposed <- proposal_myencod(current)
#   prob <- acceptance_prob(proposed, current)
#   
#   if(runif(1) < prob){
#     current <- chain[[counter + 1]] <- proposed
#   } else{
#     current <- chain[[counter + 1]] <- current
#   }
#   counter <- counter + 1
# }
# 
# 
# 
# Flist_probs <- sapply(F.list7, distance_prob)
# 
# encod_list <- lapply(F.list7, my_encod)
# 
# 
# empirical<-rep(0,length(encod_list))
# 
# for (j in 1:length(empirical)){
#   empirical[j]<-sum(chain)
# }
# 
# match_fn <- function(encod, list){
#   which(sapply(list, function(u){identical(u, encod)}))
# }
# 
# subchain <- 1*(MAX)
# 
# chain_matches <- sapply(chain[subchain], match_fn, list = encod_list)
# 
# empirical_frac <- table(chain_matches)/length(chain_matches)
# 
# plot(as.numeric(empirical_frac), Flist_probs/sum(Flist_probs))
# abline(0,1)
# 
# 
# ##Understanding distance-based distribution
# 
# 
# data("F-lists")
# Flist_probs <- sapply(F.list7, distance_prob)
# M<-kingman_m(7)
# dist<-rep(0,length(F.list7))
# for (j in 1:length(F.list7)){
#   dist[j]<-distance_Fmat(F.list7[[j]],M)
# }
# 
# Expected_distance<-sum(dist*Flist_probs/sum(Flist_probs))
# Z<-sum(Flist_probs)

find the Frechet
dist_to_all<-rep(0,length(F.list7))
for (j in 1:length(F.list7)){
  for (i in 1:length(F.list7)){
    dist_to_all[j]<-dist_to_all[j]+sum((F.list7[[j]]-F.list7[[i]])^2)^.5
  }
}

fmean<-which.min(dist_to_all)

Flist_probs<-sapply(F.list7, distance_prob,centre=F.list7[[40]],beta=.1)
variational<-Flist_probs/sum(Flist_probs)
plot(sort(variational),type="S")

hist(variational,breaks=8)
mean<-matrix(0,nrow=nrow(F.list7[[1]]),ncol=ncol(F.list7[[1]]))
probs<-0
for (j in 1:length(F.list7)){
  mean<-mean+F.list7[[j]]*distance_prob(F.list7[[j]],centre=F.list7[[40]],beta=.1)
  probs<-probs+distance_prob(F.list7[[j]],centre=F.list7[[40]],beta=.1)
}
mean<-mean/probs

Flist_probs2<-sapply(F.list7, distance_prob,centre=mean,beta=.1)
variational<-Flist_probs2/sum(Flist_probs2)
points(sort(variational),type="l",col="blue")

The proof can be similar to softmax Gumbel Ryan Adams 2013.

1. Notes: transform both beta and M entries to ensure positiveness -after conference
2. code the elbo and plot trajectory
3. Validate new code with n=7.
4. Simulated annealing using Mackenzie's chain. But how can we use our Variational distribution that we just learned?
5. This is actually ok, since remember we want to estiamte the posterior,
so we are now sampling from the variational distribution and hopefully
we will get the right Frechet mean with higher probability.

and here, we can compare with MCMC done with BEAST.

6. Add in discussion a comment on incorporating constraints. This "could" reduce the complexity
7. The proposed variational distribution can also be used in an importance sampling approach for model comparison and marginal prob. approx. (Bayes factors)
