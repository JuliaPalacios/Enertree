# CONTAINS HELPERS TO INTERACT WITH PHYLOSTAN AND
# MANIPULATE ITS OUTPUT.


library("ape")
library(phylodyn)
library("phyclust",quiet=TRUE)
library("phylotools")
library("phangorn")
library(phylotools)
library(devtools)
#updated April 2023
update_time<-function(tree, coal_times){
  #coal_times<-cumsum(coalescent.intervals(tree)$interval.length)
  #tiplabels();nodelabels()
  ##sort before update to avoid problems
  #xx1<-sort(n.t,index.return=T)
  old.edge<-tree$edge
  n.sample <- tree$Nnode + 1
  t.tot <- max(ape::node.depth.edgelength(tree))
  n.t <- t.tot - ape::node.depth.edgelength(tree) ##gives the node length
  #n.t[1:n.sample]<-0
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
  check<-new.edge[new.edge[,2]>n,]
  check2<-check[,1]-check[,2]
  while (sum(check2[check2>0])>0){
    maxcon<-which.max(check[,1]-check[,2])
    changeto<-check[maxcon,2]
    changefrom<-check[maxcon,1]
    new.edge2<-new.edge
    new.edge2[new.edge[,1]==changefrom,1]<-changeto
    new.edge2[new.edge[,1]==changeto,1]<-changefrom
    new.edge2[new.edge[,2]==changefrom,2]<-changeto
    new.edge2[new.edge[,2]==changeto,2]<-changefrom
    new.edge<-new.edge2 
    check<-new.edge[new.edge[,2]>n,]
    check2<-check[,1]-check[,2]
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
# update_time<-function(tree,coal_times){
#   old.edge<-tree$edge
#   n.sample <- tree$Nnode + 1
#   t.tot <- max(ape::node.depth.edgelength(tree))
#   n.t <- t.tot - ape::node.depth.edgelength(tree) ##gives the node length
#   new.n.t<-n.t
#   #order nodes according to length, then coalescent times are in reverse order
#   xx1<-sort(new.n.t,index.return=T)
#   index<-(2*n.sample-1):(n.sample+1)
#   for (j in (n.sample+1):(2*n.sample-1)){
#     old.edge[which(tree$edge[,1]==xx1$ix[j]),1]<-index[j-n.sample]
#     old.edge[which(tree$edge[,2]==xx1$ix[j]),2]<-index[j-n.sample]
#   }
#   new.n.t[(n.sample+1):(2*n.sample-1)]<-rev(coal_times)
#   #If we sort them, we can get the correspondence between leaves and coal. times
#   xx<-sort(new.n.t,index.return=T)
#   new.edge<-old.edge
#   for (j in (n.sample+1):(2*n.sample-1)){
#     new.edge[which(old.edge[,1]==xx$ix[j]),1]<-index[j-n.sample]
#     new.edge[which(old.edge[,2]==xx$ix[j]),2]<-index[j-n.sample]
#   }
#   new.edge.length<-new.n.t[new.edge[,1]]-new.n.t[new.edge[,2]]
#   new.tree<-tree
#   new.tree$edge<-new.edge
#   new.tree$edge.length<-new.edge.length
#   tree2<-write.tree(new.tree)
#   new.tree2<-read.tree(text=tree2)
#   trees <- c(tree, new.tree2) 
#   trees <- .compressTipLabel(trees)
#   t2 <- trees[[2]]
#   return(t2)
# }

exp_traj<-function(t){
  result=rep(0,length(t))
  result[t<=0.1]<-10
  result[t>0.1 & t<0.25]<-10*exp(2-20*t[t>0.1 & t<0.25]) #50 times smaller
  result[t>=0.25]<-0.5
  return(result)
}
#
# Makes call to phylostan command-line interface. 
# Implicitly assumes isochronous sampling, Jukes-Cantor substitution
# model, strict clock, constant coalescent.
#
# @param data_filename: filename of the sequence data.
# @param tree_filename: filename of the tree topology in newick format. With branch lengths also acceptable.
# @param output_prefix: arbitrary name of the generated files.
#
# @returns nothing.

infer_coal_times <- function(data_filename,
                             tree_filename,
                             output_prefix) {
  
  "To do
    -  determine how to set mutation rate
    -  make the phylostan command more flexible (e.g. prior on Ne and mutation rate"
  
  build_command <- paste0("phylostan build -s ",
                          output_prefix,
                          ".stan -m JC69 -C 4 --clock strict --coalescent skyride --compile")
  system(build_command)
  run_command <- paste0("phylostan run -s ",
                        output_prefix,
                        ".stan -m JC69 -C 4 --clock strict --coalescent skyride -i ",
                        data_filename,
                        " -t ",
                        tree_filename,
                        " -o ",
                        output_prefix,
                        " -q meanfield")
                        
  system(run_command)
}

#
#
# Reads inferred coalescent times from the output of phylostan
# command-line interface.
#
# @param filename: name of file containing VB samples.
# @param num_tips: number of tips.
#
# @returns: a sorted list of coalescent times of length num_tips-1.
#
read_coal_times <- function(filename, num_tips) {
  samples_txt <- scan(file=filename,
                      what=character(),
                      sep="\n")
  filtered <- Filter(f=function(line) !grepl("#", line),
                     x=samples_txt)
  # take the final sample and legend.
  legend <- filtered[1]
  sample <- filtered[length(filtered)]
  
  legend <- unlist(strsplit(legend, ","))
  sample <- unlist(strsplit(sample, ","))
  
  # this takes everything from the output.
  # from it we extract heights.(num_tips + 1) 
  # through height.(2*num_tips - 1)
  all_fields <- setNames(sample, legend)
  heights <- c()
  for (i in (num_tips+1):(2*num_tips-1)) {
    heights <- c(heights, paste0("heights.", i))
  }
  coal_times <- as.double(unname(all_fields[heights]))
  sorted_coal_times <- sort(coal_times)
  return(sorted_coal_times)
}

read_eff_pop_size <- function(filename, coal_times){
  n<- length(coal_times)+1
  "Same as read_coal_times but for the effective pop size"
  "To do
    -  verify where is exactly the grid
    -  at the moment is located at the coal_times half point"

  data<-read.delim(file=filename,header=TRUE,skip=20,sep=",")
  data <- data[-c(1,2),] #remove 2 lines of NA
  logeff <- data[ , grepl( "thetas." , names(data) )] 
  eff <- exp(logeff)
  #Grid (quit arbitrary, I don't believe it is the correct one)
  n.grid <- dim(eff)[2]
  grid.length <- coal_times[n-1]/n.grid
  mid <- (c(0,coal_times[1:(n.grid-1)])+coal_times)/2
  #mid <- coal_times
  
  medianeff <- apply(eff,2, median,na.rm=TRUE)
  q97.5 <- apply(eff,2, quantile,0.975,na.rm=TRUE)
  q2.5 <- apply(eff,2, quantile,0.025,na.rm=TRUE)
  return(list(medianeff=medianeff,q97.5=q97.5,q2.5=q2.5,mid=mid))
}





