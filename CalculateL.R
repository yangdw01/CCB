###################################################################################################
# calculate L for ARL0 = 100 (Wilson score confidence interval)
# -------------------------------------------------------------------------------------------------
library(abind)
source("DataGen.R")
load("Part0_result//PART0result_500.RData")
load("Part0_result//PART0result_500_freq.RData")
# -------------------------------------------------------------------------------------------------

# functions ---------------------------------------------------------------------------------------
calculate_RL <- function(LCL, UCL, PS){
  
  tmp_rl <- ((LCL > PS) | (UCL < PS)) %>% which
  if(length(tmp_rl) == 0)
  { 
    tmp_rl <- length(PS)
    
  }else{
    
    tmp_rl <- tmp_rl[1]
  }
  return(tmp_rl)
}

calculate_CL <- function(probs, tmpL, tmpN, StdVal = 10^(-3), case){
  
  if(case == 1){ # Wald interval
    
    ps <- probs 
    qs <- 1-probs
    
    tmpUCL <- ps + tmpL * sqrt(ps * qs / tmpN)
    tmpLCL <- ps - tmpL * sqrt(ps * qs / tmpN)
    
  }else if(case == 2){ # Agresti-Coulli interval
    
    adj_Ns <- tmpN + tmpL^2
    adj_Xs <- tmpN * probs + tmpL^2/2
    adj_ps <- adj_Xs / adj_Ns
    adj_qs <- 1 - adj_ps
    
    tmpUCL <- adj_ps + tmpL * sqrt(adj_ps * adj_qs / adj_Ns)
    tmpLCL <- adj_ps - tmpL * sqrt(adj_ps * adj_qs / adj_Ns)
    
  }else if(case == 3){ # Jeffrey interval
    
    tmpa <- tmpN * probs + 1/2
    tmpb <- tmpN - tmpN * probs + 1/2
    
    tmpUCL <- qbeta(  pnorm(tmpL), tmpa, tmpb)
    tmpLCL <- qbeta(1-pnorm(tmpL), tmpa, tmpb)
    
    close0_idx <- probs < StdVal
    close1_idx <- probs > 1- StdVal
    
    tmpLCL[close0_idx] <- 0
    tmpUCL[close1_idx] <- 1
    
  }else if(case == 4){ # Wilson score interval
    
    ps <- probs 
    qs <- 1-probs
    
    adj_Ns <- tmpN + tmpL^2
    adj_Xs <- tmpN * probs + tmpL^2/2
    
    tmpUCL <- adj_Xs/adj_Ns + tmpL * (sqrt(tmpN)/adj_Ns) * sqrt(ps*qs + tmpL^2/(4*tmpN))
    tmpLCL <- adj_Xs/adj_Ns - tmpL * (sqrt(tmpN)/adj_Ns) * sqrt(ps*qs + tmpL^2/(4*tmpN))
    
    close0_idx <- probs < StdVal
    close1_idx <- probs > 1- StdVal
    
    tmpLCL[close0_idx] <- 0
    tmpUCL[close1_idx] <- 1
    
  }else if(case == 5){ # Clopper–Pearson interval
    
    xs <- tmpN * probs
    
    tmpUCL <- 1/(1 + (tmpN - xs)/(xs+1)/qf(pnorm(tmpL), 2*(xs+1), 2*(tmpN - xs)))
    tmpLCL <- 1/(1 + (tmpN - xs + 1)/(xs)/qf(1-pnorm(tmpL), 2*xs, 2*(tmpN - xs + 1)))
    
    close0_idx <- probs < StdVal
    close1_idx <- probs > 1- StdVal
    
    tmpLCL[close0_idx] <- 0
    tmpUCL[close1_idx] <- 1
  }
  
  return(rbind(tmpLCL, tmpUCL))
}

calculate_CLall <- function(probs, tmpL, tmpN, StdVal = 10^(-3)){
  
  Result <- list()
  for(i in 1:5){ Result[[i]] <- calculate_CL(probs=probs, tmpL=tmpL, tmpN=tmpN, StdVal=StdVal, case=i) }
  return(Result)
}
# -------------------------------------------------------------------------------------------------

# calculate L -------------------------------------------------------------------------------------
Lvalues <- seq(1, 10, 0.01)

Bayesian_ARLs <- Frequentist_ARLs <- NULL
Phase1_clustering_Bayesian <- Phase1_clustering_freq <- NULL
for(tmp_L in Lvalues)
{
  cat(tmp_L, "    ")
  tmp_Bayesian_RLs <- tmp_Frequentist_RLs <- NULL
  for(i in seq(output))
  {
    tmp_output  <- output[[i]]
    tmp_output2 <- output2[[i]]
    
    # baeysian 
    tmp_CLs <-  calculate_CLall(probs=tmp_output$Bayesian[[1]],    tmpL=tmp_L, tmpN=Subgroup_SS)
    tmp_rls <-  sapply(1:length(tmp_CLs), function(j){ sapply(seq(ncol(tmp_CLs[[j]])), 
                                                             function(x){ calculate_RL(tmp_CLs[[j]][1,x], tmp_CLs[[j]][2,x], tmp_output$Bayesian[[3]][,x]) }) })
    # frequentist
    tmp_prob <- (tmp_output2$centerline %>% sort(decreasing = T))[tmp_output2[[3]] %>% apply(2,mean) %>% order(decreasing = T)]
    tmp_CLs2 <- calculate_CLall(probs=tmp_prob,      tmpL=tmp_L, tmpN=Subgroup_SS)
    tmp_rls2 <- sapply(1:length(tmp_CLs2), function(j){ sapply(seq(ncol(tmp_CLs2[[j]])), 
                                                             function(x){ calculate_RL(tmp_CLs2[[j]][1,x], tmp_CLs2[[j]][2,x], tmp_output2[[3]][,x]) }) })
    # save
    tmp_Bayesian_RLs    <- abind(tmp_Bayesian_RLs,    tmp_rls, along=3)
    tmp_Frequentist_RLs <- abind(tmp_Frequentist_RLs, tmp_rls2, along=3)
  }
  Bayesian_ARLs       <- abind(Bayesian_ARLs,       tmp_Bayesian_RLs    %>% apply(c(2,3), min),  along=3)
  Frequentist_ARLs    <- abind(Frequentist_ARLs,    tmp_Frequentist_RLs %>% apply(c(2,3), min),  along=3)
}
Bayesian_result    <- Bayesian_ARLs    %>% apply(c(1,3), mean)
Frequentist_result <- Frequentist_ARLs %>% apply(c(1,3), mean)

target_ARL0 <- 100

Bayesian_L    <- sapply(1:nrow(Bayesian_result),    function(x){ Lvalues[(Bayesian_result[x,]>target_ARL0)][1] })
Frequentist_L <- sapply(1:nrow(Frequentist_result), function(x){ Lvalues[(Frequentist_result[x,]>target_ARL0)][1] })
cat("Bayesian : ", Bayesian_L)
cat("Frequentist : ", Frequentist_L)
# -------------------------------------------------------------------------------------------------

# clustering --------------------------------------------------------------------------------------
Bayesian_clustering    <- rep(1, length(output))
Frequentist_clustering <- sapply(seq(output2), function(i){ output2[[i]]$clustering[1] })
# -------------------------------------------------------------------------------------------------
#
###################################################################################################

rbind(Bayesian_clustering, Frequentist_clustering) %>% t %>% boxplot
