





estimate_topology_coaltimes <- function(g_est, #Topology parameters 
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
                                         loops) #Joint parameter
{
 
  
    #Initialization
  if (upgma_init==TRUE){
    dm <- dist.ml(data,"JC69") ##check, this is giving the problems
    dm[is.nan(dm)] <- 0.00001
    dm[dm==Inf] <- 0.00001
    upgma<-upgma(dm)
    #plot(upgma)
    #Initial coalescent times
    inter_coal_times_est <-coalescent.intervals(upgma)$interval.length
    inter_coal_times_est[which(inter_coal_times_est<=0.001)]<-.01
    coal_times_est <- cumsum(inter_coal_times_est)
    #Initial topology
    #Julia: we shouldn't need rounding in this part
    M_est<-gen_Fmat(upgma,tol=8)
   # M_est <- round(gen_Fmat(upgma, tol=8),0)
  } #Possibly add input coal_times and M_est
    
    for (i in 1:loops){#This should be substituted with a stopping rule
    "here start the topology estimation"
    M_est_tree <- mytree_from_F(nearby_Fmat(M_est), coal_times_est)
    out_topology <- estimate_topology2(M_est,
                                       g_est, 
                                      num_tips, 
                                      num_samps,
                                      burnin,
                                      data,
                                      coal_times,
                                      step_size,
                                      num_grad_desc_steps,
                                      num_tip_label_iters,
                                      joint=TRUE)
    #New estimated tree
    g_est <- out_topology$g_est
    M_est <- out_topology$M_est
    M_estimated_tree <- mytree_from_F(nearby_Fmat(M_est), coal_times_est)
    
    "Here starts the times estimate"
    ### Coalescent times estimations
    write.tree(M_estimated_tree, 
               file=paste0(input_prefix, '.newick'))
    file.remove(paste0(output_prefix, '.stan'))
    #Inference times
    infer_coal_times(data_filename, tree_filename,output_prefix)
    #NewEstimated times
    coal_times_est <- read_coal_times(filename=output_prefix,
                                      num_tips=num_tips)
    
    }
    
    return(list(M_est=M_est,g_est=g_est,coal_times_est=coal_times_est))

}



