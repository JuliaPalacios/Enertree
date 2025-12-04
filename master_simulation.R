setwd("~/My Drive/Statistics/Phylo-VB/")
setwd("~/Documents/Phylo-VB/Phylo-VB")


library("phylodyn")
library("ape")
library("phyclust",quiet=TRUE)
library("phylotools")
library("phytools")

source("estimate_coal_utils.R")
source("./estimate_topology_utils.R")
source("./expectations.R")
source("./estimate_topology_coaltimes_utils.R")

exp_traj<-function(t){
  result=rep(0,length(t))
  result[t<=0.1]<-10
  result[t>0.1 & t<0.25]<-10*exp(2-20*t[t>0.1 & t<0.25]) #50 times smaller
  result[t>=0.25]<-0.5
  return(result)
}


unif_traj<-function(t){
  result=rep(0.5,length(t))
  return(result)
}



set.seed(123)
samp_times<-0
num_tips<- 10


#Data generation
init_results <- generate_true_M_and_data(num_tips,exp_traj)
inter_coal_times <- coalescent.intervals(init_results$M_true_tree)$interval.length
inter_coal_times[which(inter_coal_times<=0.000001)] <- 0.0000001
coal_times <- cumsum(inter_coal_times)


#Common input
input_prefix <- "trial"
output_prefix <- "output"
data_filename="out.fasta"


#######################################################
########## INFERENCE TIMES ONLY (TRUE TREE)  ##########
#######################################################

####Inputs for times
write.tree(init_results$M_true_tree, 
           file=paste0(input_prefix, '.newick'))
tree_filename <- paste0(input_prefix, '.newick')
file.remove(paste0(output_prefix, '.stan'))

#Inference times
infer_coal_times(data_filename, tree_filename,output_prefix)

##For Julia:
#build_command <- paste0("/Users/juliapr/opt/anaconda3/bin/phylostan build -s ",
                        output_prefix,
                        ".stan -m JC69 -C 4 --clock strict --coalescent skyride --compile")
#system(build_command)
# run_command <- paste0("/Users/juliapr/opt/anaconda3/bin/phylostan run -s ",
#                       output_prefix,
#                       ".stan -m JC69 -C 4 --clock strict --coalescent skyride -i ",
#                       data_filename,
#                       " -t ",
#                       tree_filename,
#                       " -o ",
#                       output_prefix,
#                       " -q meanfield")

#system(run_command)


####Estimated times
coal_times_est <- read_coal_times(filename=output_prefix,
                                  num_tips=num_tips)

###Estimated eff_times
Ne <- read_eff_pop_size(filename=output_prefix, coal_times_est)

#Plot effpop size
res<-BNPR(init_results$M_true_tre,prec_alpha=1)
plot_BNPR(res)
lines(Ne$mid,Ne$medianeff,lwd=2,type="l",ylim=c(0.9*min(Ne$q2.5),1.2*max(Ne$q97.5)),col="red")
lines(Ne$mid,Ne$q97.5,lwd=2,col="grey")                
lines(Ne$mid,Ne$q2.5,lwd=2,col="grey")                
points(seq(0,2.5,by=.1),exp_traj(seq(0,2.5,by=.1)),type="l",col="blue")


#######################################################
########## INFERENCE TIMES AND TOPOLOGY    ############
#######################################################


#Inputs
num_samps<-100
data <- init_results$sequences
g_est <- log(1.5)
burnin <- round(num_samps*0.1)
num_grad_desc_steps <- 50
num_tip_label_iters <- 30
num_gibbs_iters <- 50
loops <- 2


#Run the algorithm

out <- estimate_topology_coaltimes(g_est, #Topology parameters 
                                        num_tips, 
                                        num_samps,
                                        burnin,
                                        data,
                                        step_size,
                                        num_grad_desc_steps,
                                        num_tip_label_iters,
                                        joint=TRUE,
                                        upgma_init=TRUE,
                                        data_filename, #Times parameters
                                        tree_filename,
                                        output_prefix,
                                        loops)#Joint parameters


####Estimated times
M_est <- out$M_est

####Estimated times
coal_times_est <- out$coal_times_est
M_estimated_tree <- mytree_from_F(nearby_Fmat(M_est), coal_times_est)

###Estimated eff_times
Ne <- read_eff_pop_size(filename=output_prefix, coal_times_est)

#Plot effpop size
plot(Ne$mid,Ne$medianeff,lwd=2,type="l",ylim=c(0.9*min(Ne$q2.5),1.2*max(Ne$q97.5)))
lines(Ne$mid,Ne$q97.5,lwd=2,col="grey")                
lines(Ne$mid,Ne$q2.5,lwd=2,col="grey")                
points(seq(0,2.5,by=.1),exp_traj(seq(0,2.5,by=.1)),type="l",col="blue")





