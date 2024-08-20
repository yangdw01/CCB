###################################################################################################
###################################################################################################
####
####
####      JQT         by DY
####
####      (Simulation 2) - PART0. Selection of Control limit L
####
####
###################################################################################################
###################################################################################################
#
# Code explanations :
# 
#   Choose an appropriate control limit L with ARL0 = 100
#
###################################################################################################



###################################################################################################
# library
# -------------------------------------------------------------------------------------------------
library(dplyr)
library(reshape2)
library(tidyr)
library(fda)
library(LaplacesDemon)
library(fossil)
library(foreach)
library(doParallel)
library(combinat)
# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# source
# -------------------------------------------------------------------------------------------------
source("DataGen.R")
source("Functions.R")
# -------------------------------------------------------------------------------------------------
#
###################################################################################################

###################################################################################################
# parallel computing
# -------------------------------------------------------------------------------------------------
# basic setting
niter <- 500

# check # of cores
detectCores()

# create the cluster
myCluster <- makeCluster(15)

# register the cluster with the foreach package
registerDoParallel(myCluster)

# parallel computing
output2 <- foreach(i = 1:niter, .packages = c("dplyr","tidyr","fda","LaplacesDemon","fossil","combinat")) %dopar% {
  
  cat(i,"th iteration !!! \n")
  
  ###################################################################################################
  # data generation
  # -------------------------------------------------------------------------------------------------
  isPART0 <- T
  seed_num <- i
  DATA <- DataGen(isPART0=isPART0, seednum=seed_num, delta=NULL, scenario=NULL)
  # -------------------------------------------------------------------------------------------------
  #
  ###################################################################################################
  
  
  ###################################################################################################
  # (Step1) FDA - with regularization
  # -------------------------------------------------------------------------------------------------
  # lambda selection
  ## reshaping data
  tmp_DATA <- DATA %>% filter(time <= Phase1_SS)
  tmp_data <- tmp_DATA %>% pivot_wider(names_from = domain1:domain2, values_from = y) %>% dplyr::select(-c("time","CurveNUM","clusterNUM")) %>% as.matrix() %>% t()
  
  ## obtain lambda
  loglam = seq(-6, 0, 0.25)
  gcv_result <- NULL
  for(i in 1:length(loglam)){
    
    tmp_lambda <- 10^loglam[i]
    tmp_basis_coefficient <- solve(t(basis_matrix) %*% basis_matrix + tmp_lambda * Om) %*% t(basis_matrix) %*% tmp_data
    tmp_fhat <- basis_matrix %*% tmp_basis_coefficient 
    tmp_Smat <- basis_matrix %*%solve(t(basis_matrix) %*% basis_matrix + tmp_lambda * Om) %*% t(basis_matrix)
    tmp_gcv <- sapply(1:ncol(tmp_fhat), function(x){ ((tmp_data[,x] - tmp_fhat[,x])/(1-diag(tmp_Smat)))^2 %>% mean() }) %>% sum()
    
    gcv_result <- c(gcv_result, tmp_gcv)
  }
  
  lambda <- 10^(loglam[gcv_result %>% which.min])
  basis_coefficient <- solve(t(basis_matrix) %*% basis_matrix + lambda * Om) %*% t(basis_matrix) %*% tmp_data
  fhat <- basis_matrix %*% basis_coefficient 
  Smat <- basis_matrix %*%solve(t(basis_matrix) %*% basis_matrix + lambda * Om) %*% t(basis_matrix)
  gcv <- sapply(1:ncol(fhat), function(x){ ((tmp_data[,x] - fhat[,x])/(1-diag(Smat)))^2 %>% mean() }) %>% sum()
  # -------------------------------------------------------------------------------------------------
  
  # calculate coefficients --------------------------------------------------------------------------
  tmp_data <- DATA %>% pivot_wider(names_from = domain1:domain2, values_from = y) %>% dplyr::select(-c("time","CurveNUM","clusterNUM")) %>% as.matrix() %>% t()
  basis_coefficient <- solve(t(basis_matrix) %*% basis_matrix + lambda * Om) %*% t(basis_matrix) %*% tmp_data
  fhat <- basis_matrix %*% basis_coefficient 
  Smat <- basis_matrix %*%solve(t(basis_matrix) %*% basis_matrix + lambda * Om) %*% t(basis_matrix)
  
  ## obtain coefficients of fd object
  Step1_result <- basis_coefficient %>% t()
  # -------------------------------------------------------------------------------------------------
  #
  ###################################################################################################
  
  
  # ###################################################################################################
  # # (Step 2) K-means algorithm
  # # -------------------------------------------------------------------------------------------------
  # # K-means algorithm
  # set.seed(1)
  # 
  # Step1_result_Kmeans <- kmeans(Step1_result, centers=3, iter.max=10000)
  # true_clusters <- DATA %>% distinct(time, CurveNUM, clusterNUM) %>% pull(clusterNUM)
  # 
  # rd_result0 <- rand.index(Step1_result_Kmeans$cluster, true_clusters)
  # ard_result0 <- adj.rand.index(Step1_result_Kmeans$cluster, true_clusters)
  # rd_result0
  # ard_result0
  # # -------------------------------------------------------------------------------------------------
  # #
  # ###################################################################################################
  
  
  # ###################################################################################################
  # # (Step 2) SUGS algorithm
  # # -------------------------------------------------------------------------------------------------
  # # basic setting
  # Y = Step1_result
  # M <- 15
  # SUGS_ordering_num <- 10
  # 
  # # priors
  # a = 1
  # m = 0
  # psi = 1
  # c = 1
  # d = 10
  # alpha_stars <- c(0.01, 0.05, c(0.1 + 0.2 * (0:19)))
  # eta_probs <- rep(1/length(alpha_stars), length(alpha_stars))
  # 
  # # online learning
  # Final_ordering <- NULL
  # multinomial_prob <- NULL
  # tmp_SSlist <- rep(Subgroup_SS, Phase1_SS + MC_ARL0_SS)
  # for(w in 1:(Phase1_SS + MC_ARL0_SS)){
  #   
  #   tmp_idx <- ((w-1) * Subgroup_SS) + 1:Subgroup_SS
  #   
  #   tmp_result_list <- list()
  #   tmp_ordering_list <- NULL 
  #   
  #   for(ww in 1:SUGS_ordering_num){
  #     
  #     cat(w, "th /", Phase1_SS + MC_ARL0_SS, " : ", ww, "th ordering \n")
  #     
  #     # ordering
  #     if(ww == 1){
  #       
  #       tmp_resample <- seq(tmp_idx)
  #       
  #     }else{
  #       
  #       set.seed(ww + w*(Phase1_SS + MC_ARL0_SS + SUGS_ordering_num))  
  #       tmp_resample <- sample(seq(tmp_idx), length(tmp_idx))
  #     }
  #     tmp_Y <- Y[tmp_idx[tmp_resample],]
  #     
  #     # SUGS algorithm
  #     if(w==1){
  #       
  #       tmp_result <- SUGS_algorithm_sequential_fast(Y=tmp_Y, a=a, m=m, psi=psi, c=c, d=d, alpha_stars=alpha_stars, eta_probs=eta_probs, first_iter=TRUE,
  #                                                    prev_bhat_j=NULL, prev_Gamma_i=NULL, prev_phi_t=NULL, prev_psi_hj=NULL, prev_m_hj=NULL, prev_a_hj=NULL, prev_b_hj=NULL)
  #     }else{
  #       
  #       tmp_result <- SUGS_algorithm_sequential_fast(Y=tmp_Y, a=a, m=m, psi=psi, c=c, d=d, alpha_stars=alpha_stars, eta_probs=eta_probs, first_iter=FALSE,
  #                                                    prev_bhat_j=tmp_result__bhat_j, prev_Gamma_i=tmp_result__Gamma_i, prev_phi_t=tmp_result__phi_t, 
  #                                                    prev_psi_hj=tmp_result__psi_hj, prev_m_hj=tmp_result__m_hj, prev_a_hj=tmp_result__a_hj, prev_b_hj=tmp_result__b_hj)
  #     }
  #     
  #     tmp_result_list[[ww]] <- tmp_result
  #     tmp_ordering_list <- rbind(tmp_ordering_list, tmp_idx[tmp_resample])
  #   }
  #   max_idx <- sapply(1:SUGS_ordering_num, function(x){ tmp_result_list[[x]][[1]] }) %>% which.max()
  #   
  #   # SUGS algorithm result
  #   tmp_result <- tmp_result_list[[max_idx]]
  #   tmp_result__bhat_j   <- tmp_result[[2]]
  #   tmp_result__Gamma_i  <- tmp_result[[3]]
  #   tmp_result__phi_t    <- tmp_result[[4]]
  #   tmp_result__psi_hj   <- tmp_result[[5]]
  #   tmp_result__m_hj     <- tmp_result[[6]]
  #   tmp_result__a_hj     <- tmp_result[[7]]
  #   tmp_result__b_hj     <- tmp_result[[8]]
  #   
  #   # for fast calculation
  #   if(w > Phase1_SS + 1)
  #   {
  #     delete_idx <- (1:Subgroup_SS) + Phase1_SS*Subgroup_SS
  #     
  #     tmp_result__bhat_j  <- tmp_result__bhat_j[-delete_idx,]
  #     tmp_result__phi_t   <- tmp_result__phi_t[-delete_idx,]
  #     tmp_result__psi_hj  <- sapply(seq(tmp_result__psi_hj), function(x){ tmp_result__psi_hj[[x]][-delete_idx,] }, simplify = F)
  #     tmp_result__m_hj    <- sapply(seq(tmp_result__m_hj),   function(x){   tmp_result__m_hj[[x]][-delete_idx,] }, simplify = F)
  #     tmp_result__a_hj    <- sapply(seq(tmp_result__a_hj),   function(x){   tmp_result__a_hj[[x]][-delete_idx,] }, simplify = F)
  #     tmp_result__b_hj    <- sapply(seq(tmp_result__b_hj),   function(x){   tmp_result__b_hj[[x]][-delete_idx,] }, simplify = F)
  #   }
  #   
  #   # ordering
  #   Final_ordering <- c(Final_ordering, tmp_ordering_list[max_idx,])
  # }
  # 
  # # SUGS algorithm result
  # Final_result___bhat_j  <- tmp_result__bhat_j
  # Final_result___Gamma_i <- tmp_result__Gamma_i
  # Final_result___phi_t   <- tmp_result__phi_t
  # Final_result___psi_hj  <- tmp_result__psi_hj
  # Final_result___m_hj    <- tmp_result__m_hj
  # Final_result___a_hj    <- tmp_result__a_hj
  # Final_result___b_hj    <- tmp_result__b_hj
  # Gamma_i <- Final_result___Gamma_i
  # 
  # # clustering result
  # true_clusters <- DATA %>% distinct(time, CurveNUM, clusterNUM) %>% pull(clusterNUM)
  # rd_result  <- rand.index(Gamma_i, true_clusters[Final_ordering])
  # ard_result <- adj.rand.index(Gamma_i, true_clusters[Final_ordering])
  # # -------------------------------------------------------------------------------------------------
  # #
  # ###################################################################################################
  
  
  ###################################################################################################
  # (Step 3) p-chart
  # -------------------------------------------------------------------------------------------------
  # # p chart based on Bayesian approach
  # ## parameter estimation
  # Phase1_idx <- tmp_SSlist[1:Phase1_SS] %>% sum %>% seq
  # Phase1_alpha_prob <- (Final_result___phi_t[Phase1_idx,] %>% tail(1) %>% c())
  # Phase1_alpha <- ((Final_result___phi_t[Phase1_idx,] %>% tail(1) %>% c()) * alpha_stars) %>% sum()
  # Phase1_Gamma_i <- Final_result___Gamma_i[Phase1_idx]
  # Phase1_Gamma_params <- sapply(1:M, function(k){ c(sum(Phase1_Gamma_i == k),sum(Phase1_Gamma_i > k)) })
  # Phase1_pks <- sapply(1:M, function(k){ rbeta(1, 1+Phase1_Gamma_params[1,k], Phase1_Gamma_params[2,k]+Phase1_alpha) })
  # Phase1_piks <- sapply(1:M, function(k){ Phase1_pks[k] * prod(c(1,1-Phase1_pks)[1:k]) }) %>% unlist()
  # Phase1_piks <- c(Phase1_piks, 1-sum(Phase1_piks))
  # 
  # ## plotting statistics
  # multinomial_result <- NULL
  # for(w in 1:MC_ARL0_SS){
  #   
  #   tmp_idx <- ((Phase1_SS + w-1) * Subgroup_SS) + 1:Subgroup_SS
  #   tmp_Gamma <- Final_result___Gamma_i[tmp_idx]
  #   tmp_multinomial_result <- sapply(1:(M+1), function(k){ sum(tmp_Gamma == k)/Subgroup_SS })
  #   
  #   multinomial_result <- rbind(multinomial_result, tmp_multinomial_result)
  # }
  # 
  # ## sigma
  # Phase1_confints_sigma <- sapply(1:(M+1), function(x){sqrt(Phase1_piks[x] * (1-Phase1_piks[x]) / Subgroup_SS)})
  
  # p chart based on Frequentist approach
  ## Phase 1 - parameter estimation
  set.seed(1)
  Phase1_idx <- seq(Subgroup_SS * Phase1_SS)
  Step1_result_Kmeans <- kmeans(Step1_result, centers=3, iter.max=10000, nstart=10) # Kmeans_custom(Step1_result)
  Phase1_cluster_props <- (Step1_result_Kmeans$cluster[Phase1_idx] %>% table())/(Subgroup_SS * Phase1_SS)
  Phase1_cluster_props_sigma <- sqrt(Phase1_cluster_props * (1-Phase1_cluster_props) / Subgroup_SS)
  true_clusters <- DATA %>% distinct(time, CurveNUM, clusterNUM) %>% pull(clusterNUM)
  
  ## plotting statistics
  freq_PS <- NULL
  est_clusters <- Step1_result_Kmeans$cluster[Phase1_idx]
  perm_list <- permn(3)
  tmp_clusters_chk <- est_clusters
  tmp_ordering_result <- NULL
  for(i in seq(perm_list))
  {
    tmp_perm <- perm_list[[i]]
    tmp_tmp_clusters_chk <- tmp_clusters_chk
    for(j in 1:3){ tmp_tmp_clusters_chk[tmp_clusters_chk == j] <- tmp_perm[j] }
    tmp_ordering_result <- c(tmp_ordering_result, mean(true_clusters[Phase1_idx] == tmp_tmp_clusters_chk))
  }
  cluster_ordering <- perm_list[[tmp_ordering_result %>% which.max]]
  tmp_final_clusters <- est_clusters
  for(j in 1:3){ tmp_final_clusters[est_clusters == j] <- cluster_ordering[j] }
  est_clusters <- tmp_final_clusters
  
  for(w in 1:(MC_ARL0_SS)){
    
    # clustering
    tmp_idx <- ((Phase1_SS + w-1) * Subgroup_SS) + 1:Subgroup_SS
    tmp_Step1_result <- Step1_result[tmp_idx,]
    tmp_Step1_result_Kmeans <- kmeans(tmp_Step1_result, centers=3, iter.max=10000, nstart=10) # Kmeans_custom(tmp_Step1_result)
    tmp_clusters <- tmp_Step1_result_Kmeans$cluster
    
    # check cluster ordering
    perm_list <- permn(3)
    tmp_clusters_chk <- tmp_clusters
    tmp_ordering_result <- NULL
    for(i in seq(perm_list))
    {
      tmp_perm <- perm_list[[i]]
      tmp_tmp_clusters_chk <- tmp_clusters_chk
      for(j in 1:3){ tmp_tmp_clusters_chk[tmp_clusters_chk == j] <- tmp_perm[j] }
      tmp_ordering_result <- c(tmp_ordering_result, mean(true_clusters[tmp_idx] == tmp_tmp_clusters_chk))
    }
    cluster_ordering <- perm_list[[tmp_ordering_result %>% which.max]]
    tmp_final_clusters <- tmp_clusters
    for(j in 1:3){ tmp_final_clusters[tmp_clusters == j] <- cluster_ordering[j] }
    freq_PS <- rbind(freq_PS,
                     c(mean(tmp_final_clusters == 1),
                       mean(tmp_final_clusters == 2),
                       mean(tmp_final_clusters == 3)))
    #(tmp_final_clusters[-seq(Phase1_idx)] %>% table())/(Subgroup_SS))
    est_clusters <- c(est_clusters,tmp_final_clusters)
  }
  
  rd_result0 <- rand.index(est_clusters, true_clusters)
  ard_result0 <- adj.rand.index(est_clusters, true_clusters)
  
  # save result
  Part0_for_L_freq <- list(centerline = Phase1_cluster_props, sigma=Phase1_cluster_props_sigma, PS=freq_PS, clustering=c(rd_result0, ard_result0))
  # -------------------------------------------------------------------------------------------------
  #
  ###################################################################################################
  Part0_for_L_freq
}

# stop the cluster 
stopCluster(myCluster)
# -------------------------------------------------------------------------------------------------
#
###################################################################################################
save(output2, file="PART0result_500_freq.RData")



# ###################################################################################################
# # calculate L for ARL0 = 100
# # -------------------------------------------------------------------------------------------------
# # load("DY_Simulation2//PART0result.RData")
# 
# calculate_RL <- function(LCL, UCL, PS){
#   
#   tmp_rl <- ((LCL > PS) | (UCL < PS)) %>% which
#   if(length(tmp_rl) == 0)
#   { 
#     tmp_rl <- length(PS)
#     
#   }else{
#     
#     tmp_rl <- tmp_rl[1]
#   }
#   
#   return(tmp_rl)
# }
# 
# Lvalues <- seq(1,4, 0.01)
# Bayesian_ARLs <- Frequentist_ARLs <- NULL
# for(tmp_L in Lvalues)
# {
#   cat(tmp_L, "    ")
#   tmp_Bayesian_RLs <- tmp_Frequentist_RLs <- NULL
#   for(i in seq(output))
#   {
#     tmp_output <- output[[i]]
#     
#     # baeysian 
#     tmp_UCL <- tmp_output$Bayesian[[1]] + tmp_L * tmp_output$Bayesian[[2]]
#     tmp_LCL <- tmp_output$Bayesian[[1]] - tmp_L * tmp_output$Bayesian[[2]]
#     tmp_rl  <- sapply(seq(tmp_UCL), function(x){ calculate_RL(tmp_LCL[x], tmp_UCL[x], tmp_output$Bayesian[[3]][,x]) })
#     
#     # frequentist
#     tmp_UCL2 <- tmp_output$Frequentist[[1]] + tmp_L * tmp_output$Frequentist[[2]]
#     tmp_LCL2 <- tmp_output$Frequentist[[1]] - tmp_L * tmp_output$Frequentist[[2]]
#     tmp_rl2  <- sapply(seq(tmp_UCL2), function(x){ calculate_RL(tmp_LCL2[x], tmp_UCL2[x], tmp_output$Frequentist[[3]][,x]) })
#     
#     # save
#     tmp_Bayesian_RLs    <- rbind(tmp_Bayesian_RLs,    tmp_rl)
#     tmp_Frequentist_RLs <- rbind(tmp_Frequentist_RLs, tmp_rl2)
#   }
#   Bayesian_ARLs    <- rbind(Bayesian_ARLs,    tmp_Bayesian_RLs       %>% apply(2, mean))
#   Frequentist_ARLs <- rbind(Frequentist_ARLs, tmp_Frequentist_RLs    %>% apply(2, mean))
# }
# Bayesian_result <- Bayesian_ARLs %>% apply(1, min)
# Frequentist_result <- Frequentist_ARLs %>% apply(1, min)
# 
# target_ARL0 <- 100
# Bayesian_L <- Lvalues[((Bayesian_result > target_ARL0) %>% which)[1]]
# Frequentist_L <- Lvalues[((Frequentist_result > target_ARL0) %>% which)[1]]
# Bayesian_L
# Frequentist_L
# # -------------------------------------------------------------------------------------------------
# #
# ###################################################################################################