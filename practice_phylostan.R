"Practice to see if phylostan works "
setwd("~/My Drive/Statistics/Phylo-VB/")



"Delete everything  that was generated in previous runs (in particular the .stan
 and .pkl file) in the folder, otherwise you may need to recompile the software"


library(phylodyn)
library("ape")
library("phyclust",quiet=TRUE)
library("phylotools")
library("phytools")

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
n_sampled<- 30

file.remove("trial-GTR-W4.stan")
file.remove("out.fasta")
file.remove("inputtree.tree")
file.remove("output")
simul1<-coalsim(samp_times = samp_times, n_sampled = n_sampled, traj = unif_traj,method="tt",val_upper=11)
maxco <- simul1$coal_times[n_sampled-1]
simul1$coal_times[n_sampled-1]
#simul1$coal_times <- simul1$coal_times/simul1$coal_times[n_sampled-1]
#simul1$coal_times <- simul1$coal_times*maxco
#simul1$intercoal_times <- diff(c(0,simul1$coal_times))
tree<-generate_newick(simul1)
tree.newick<-write.tree(tree$newick)

seqgen(opts="-mHKY -t0.5 -f0.25,0.25,0.25,0.25 -l100",newick.tree=tree.newick,temp.file="temp.fasta")
#install.packages("phylotools")
data2<-read.phylip("temp.fasta")
data<-dat2fasta(data2,outfile="out.fasta")
write.tree(tree$newick,file="inputtree.tree")

#Phylostan part 

"Two options: 
1. Move to the terminal  and run  the two lines below

phylostan build -s trial-GTR-W4.stan  -m HKY -C 4  --clock strict --coalescent skyride --compile
phylostan run -s trial-GTR-W4.stan  -m HKY -C 4  --clock strict --coalescent skyride -i out.fasta -t inputtree.tree -o output -q meanfield



known rate 
phylostan build -s trial-GTR-W4.stan -m HKY -C 4 --heterochronous --estimate_rate --clock strict --coalescent skyride --compile
phylostan run -s trial-GTR-W4.stan -m HKY -C 4 --heterochronous --estimate_rate --clock strict --coalescent skyride -i out.fasta -t inputtree.tree -o output -q meanfield


2. Run from R. You need to substitute /opt/homebrew/bin/ with the location of phylostan in your computer. You can 
get it typing in the terminal

command -v phylostan stan

For me, it only works the first line and not the second. It seems that it is related to some messed up python installation in my mac"

system("/opt/homebrew/bin/phylostan build -s trial-GTR-W4.stan -m HKY -C 4 --clock strict --coalescent skyride")

system("/opt/homebrew/bin/phylostan run -s trial-GTR-W4.stan -m HKY -C 4 --clock strict --coalescent skyride -i out.fasta -t inputtree.tree -o output -q meanfield")

"Regardless if you did 1 or 2, you should have the output file"

#Analysis of the output
data<-read.delim(file="output",header=TRUE,skip=20,sep=",")
data <- data[-c(1,2),] #remove 2 lines of NA

#Coalescent times (from VB samples)
n <- sum(n_sampled)
heights <- data[ , grepl( "heights." , names(data) )] 
medianheigh <- apply(heights,2, median,na.rm=TRUE)
medianheigh <- medianheigh[(n+1):(n+n-1)]
coal_times <- sort(as.numeric(medianheigh))
#save <- coal_times
intercoal_times <- c(coal_times[1],diff(coal_times))

plot(simul1$coal_times,coal_times)
abline(0,1) 
#Recall that we are estimating mu but we cannot because of isochronous data





#Effective pop size
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
matplot(mid,t(eff[1:200,]),type="l",
        col = rgb(red = 0, green = 0, blue = 1, alpha = 0.1),
        ylab="eff Pop",xlab="Time",log="y") 
lines(mid,medianeff,lwd=2)
lines(mid,q97.5,lwd=2,col="grey")                
lines(mid,q2.5,lwd=2,col="grey")                
points(seq(0,2.5,by=.1),exp_traj(seq(0,2.5,by=.1)),type="l",col="blue")

      

###INLA comparison
est<-BNPR(tree$newick)
## Phylostan
phylo <- list(coal_times=coal_times,samp_times=samp_times,n_sampled=n_sampled)
est2<-BNPR(phylo)
#Plot
plot_BNPR(est,xlim=c(0,0.4))
points(seq(0,2.5,by=.1),unif_traj(seq(0,2.5,by=.1)),type="l",col="blue",xlim=c(0,1))
lines(est2$x,est2$effpop,col="red",lwd=2)
lines(est2$x,est2$effpop025,col="red",lty="dashed")
lines(est2$x,est2$effpop975,col="red",lty="dashed")
lines(mid,medianeff,lwd=2,col="green")
lines(mid,q97.5,lwd=2,col="green")                
lines(mid,q2.5,lwd=2,col="green")   

# 
# ####Heights
# heights <- data[ , grepl( "heights." , names(data) )] 
# medianheigh <- apply(heights,2, median,na.rm=TRUE)
# medianheigh
# 
# 
# 
# treedecoding <- read.tree(file="inputtree.tree")
# plot(treedecoding)
# treedecoding$edge
# 
# 
# 
# 
# 
# #####
# n <- 40
# medianheigh[(n+1):(n+n-1)]
# 
# for (i in 2:100){
# sample <- data[i,]
# sam.prop<- sample[ , grepl( "props." , names(sample) )] 
# sam.hei <- sample[ , grepl( "heights." , names(sample) )] 
# sam.hei <- sam.hei[(n+1):(n+n-1)]
# print(sum(diff(as.numeric(sam.hei))<0))
# }
# 
# sam.hei
# 0.72943*sam.prop
# 
# ###Heights from fluA
# data<-read.delim(file="../fluA/fluA",header=TRUE,skip=20,sep=",")
# heights <- data[ , grepl( "heights." , names(data) )] 
# medianheigh <- apply(heights,2, median,na.rm=TRUE)
# medianheigh
# n <- 69
# diff(medianheigh[(n+1):(n+n-1)])
# 
# 
# ###Heights from sim6
# data<-read.delim(file="../sim6/output",header=TRUE,skip=20,sep=",")
# props<- data[ , grepl( "props." , names(data) )] 
# med.props <- apply(props,2, median,na.rm=TRUE)
# med.props
# 
# 
# blens<- data[ , grepl( "blensU" , names(data) )] 
# med.blens <- apply(blens,2, median,na.rm=TRUE)
# med.blens
# 
# 
# med.heigh <- medianheigh[7:11]
# 
# 
# 
# sample <- data[4,]
# sam.prop<- sample[ , grepl( "props." , names(sample) )] 
# sam.hei <- sample[ , grepl( "heights." , names(sample) )] 
# sam.hei <- sam.hei[7:11]
# 
# as.numeric(sam.hei)*as.numeric(sam.prop[1])
# 
# 
# 
# treedecoding <- read.tree(file="../sim6/sim6.tree")
# plot(treedecoding)
# nodeHeights(treedecoding)
# 
# 
# 
