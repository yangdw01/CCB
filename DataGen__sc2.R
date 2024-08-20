###################################################################################################
###################################################################################################
####
####
####      JQT         by DY
####
####      (Simulation 2) - PART0. Data generation
####
####
###################################################################################################
###################################################################################################
#
# Code explanations :
# 
#   Generate the data for Simulation 1
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
library(plot3D)
# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# Data generation setting
# -------------------------------------------------------------------------------------------------
# true curve
## domain
M <- mesh(seq(-1, 1, length.out=26), seq(-1, 1, length.out=26))
t1 <- M$x
t2 <- M$y

## true coef
TRUE_COEF <- rbind(c(1,1,1), c(-1,2,0), c(-1,-1,-0.5), c(1,-1,-1))
TRUE_MEAN_2d <- sapply(1:nrow(TRUE_COEF), function(i){ TRUE_COEF[i,1] * sin(t1) + TRUE_COEF[i,2] * cos(t2) + TRUE_COEF[i,3] * (t1*t2) }, simplify=F)

## plot
cols <- hcl.colors(10, "YlOrRd")
contour(x=seq(-1, 1, length.out=26), y=seq(-1, 1,length.out=26), z=TRUE_MEAN_2d[[1]], col=cols)
contour(x=seq(-1, 1, length.out=26), y=seq(-1, 1,length.out=26), z=TRUE_MEAN_2d[[2]], col=cols)
contour(x=seq(-1, 1, length.out=26), y=seq(-1, 1,length.out=26), z=TRUE_MEAN_2d[[3]], col=cols)
contour(x=seq(-1, 1, length.out=26), y=seq(-1, 1,length.out=26), z=TRUE_MEAN_2d[[4]], col=cols)
# -------------------------------------------------------------------------------------------------

# basis -------------------------------------------------------------------------------------------
x <- y <- seq(-1, 1, length.out=26)
x_bsobj=create.bspline.basis(range(x),norder=4,breaks=quantile(x,prob=seq(0,1,length = 4)))
y_bsobj=create.bspline.basis(range(y),norder=4,breaks=quantile(y,prob=seq(0,1,length = 4)))
xbs=eval.basis(x,x_bsobj)
ybs=eval.basis(y,y_bsobj)

# row-wise kronecker product (Z matrix)
basis_matrix <- do.call('rbind', sapply(seq(x), function(ii){ do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[ii,i]*ybs)) }, simplify = F))

# penalty
Fmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om <- Fmat+Gmat+2*Hmat
# -------------------------------------------------------------------------------------------------

# SPC setting  ------------------------------------------------------------------------------------
## basic setting
Subgroup_SS <- 100
Phase1_SS <- 20

## setting for Part 0
MC_ARL0_SS <- 200

## setting for Part 1 & 2
Phase2_IC_SS <- 10
Phase2_OOC_SS <- 30
SSlist <- rep(Subgroup_SS, Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS)
# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# Data generation 
# -------------------------------------------------------------------------------------------------
DataGen_sc2 <- function(isPART0, seednum, delta){
  
  cat(seednum, "th data !!! \n")
  set.seed(seednum)
  
  if(isPART0 == T)
  {
    # data generation for part 0
    tmp_param_dirichlet <- c(5.5/10, 4.5/10-delta, delta, 0)*400 # 200 #
    
    # probability
    Prob_matrix <- matrix(NA, nrow=Phase1_SS + MC_ARL0_SS, ncol=4)
    for(i in seq(Phase1_SS + MC_ARL0_SS)){ Prob_matrix[i,] <- rdirichlet(1, tmp_param_dirichlet) }
    
    # cluster structure
    cluster_matrix <- matrix(NA, nrow=Phase1_SS + MC_ARL0_SS, ncol=4)
    for(i in seq(Phase1_SS + MC_ARL0_SS)){ cluster_matrix[i,] <- c(rmultinom(1, Subgroup_SS, Prob_matrix[i,])) }
    
    # data generation
    sig <- 1
    DATAlist <- list()
    for(i in seq(Phase1_SS + MC_ARL0_SS))
    {
      tmp_cluster <- sapply(1:4, function(x){ rep(x,cluster_matrix[i,x]) }) %>% unlist()
      tmp_DATA <- NULL
      for(j in 1:Subgroup_SS)
      {
        tmp_COEF <- TRUE_COEF + rnorm(length(TRUE_COEF), 0, 0.1)
        
        tmp_time <- i
        tmp_CurveNUM <- j
        tmp_clusterNUM <- tmp_cluster[j]
        tmp_y <- c(cbind(sin(c(t1)), cos(c(t2)), c(t1)*c(t2)) %*% tmp_COEF[tmp_clusterNUM,]) + rnorm(length(t1), 0, sig)
        
        tmp_data <- tibble(time=tmp_time, CurveNUM=tmp_CurveNUM, clusterNUM=tmp_clusterNUM, domain1=c(t1), domain2=c(t2), y=tmp_y)
      
        # save data
        tmp_DATA <- tmp_DATA %>% bind_rows(tmp_data)
      }
      DATAlist[[i]] <- tmp_DATA
    }
    DATA <- do.call(bind_rows, DATAlist)
  }
  return(DATA)
}
# -------------------------------------------------------------------------------------------------
#
###################################################################################################

# ###################################################################################################
# # Data generation for Part0 - choose L with ARL0 = 100
# # -------------------------------------------------------------------------------------------------
# # following variables are required to be defined : 
# # isPART0, seed_num
# # -------------------------------------------------------------------------------------------------
# if(isPART0 == T)
# {
#   set.seed(seed_num)
#   tmp_param_dirichlet <- c(5/10, 4/10, 1/10, 0)*400
# 
#   # probability
#   Prob_matrix <- matrix(NA, nrow=Phase1_SS + MC_ARL0_SS, ncol=4)
#   for(i in seq(Phase1_SS + MC_ARL0_SS)){ Prob_matrix[i,] <- rdirichlet(1, tmp_param_dirichlet) }
# 
#   # cluster structure
#   cluster_matrix <- matrix(NA, nrow=Phase1_SS + MC_ARL0_SS, ncol=4)
#   for(i in seq(Phase1_SS + MC_ARL0_SS)){ cluster_matrix[i,] <- c(rmultinom(1, Subgroup_SS, Prob_matrix[i,])) }
# 
#   # data generation
#   sig <- 1
#   DATAlist <- list()
#   for(i in seq(Phase1_SS + MC_ARL0_SS))
#   {
#     tmp_cluster <- sapply(1:4, function(x){ rep(x,cluster_matrix[i,x]) }) %>% unlist()
#     tmp_DATA <- NULL
#     for(j in 1:Subgroup_SS)
#     {
#       tmp_coef1 <- true_coef1 + rnorm(trueK, 0, 0.1)
#       tmp_coef2 <- true_coef2 + rnorm(trueK, 0, 0.1)
#       tmp_coef3 <- true_coef3 + rnorm(trueK, 0, 0.1)
#       tmp_coef4 <- true_coef4 + rnorm(trueK, 0, 0.1)
#       tmp_COEF <- rbind(tmp_coef1, tmp_coef2, tmp_coef3, tmp_coef4)
# 
#       tmp_time <- i
#       tmp_CurveNUM <- j
#       tmp_clusterNUM <- tmp_cluster[j]
#       tmp_y <- c(basis_matrix %*% tmp_COEF[tmp_clusterNUM,]) + rnorm(length(t), 0, sig)
#       tmp_data <- tibble(time=tmp_time, CurveNUM=tmp_CurveNUM, clusterNUM=tmp_clusterNUM, domain=t, y=tmp_y)
# 
#       # save data
#       tmp_DATA <- tmp_DATA %>% bind_rows(tmp_data)
#     }
#     DATAlist[[i]] <- tmp_DATA
#   }
#   DATA <- do.call(bind_rows, DATAlist)
# }
# # -------------------------------------------------------------------------------------------------
# #
# ###################################################################################################

# ###################################################################################################
# # Data generation for Part1 and 2 - ARL1 & rand index
# # -------------------------------------------------------------------------------------------------
# # following variables are required to be defined :
# # isPART0, delta, scenario, seed_num
# # -------------------------------------------------------------------------------------------------
# # probability
# ## phase 1 & IC & scenario 1,3  :   rdirichlet(1, c(5/10, 4/10, 1/10, 0)*400)
# ## phase 1 & IC & scenario 2    :   rdirichlet(1, c(5.5/10-delta/2, 4.5/10-delta/2, delta, 0)*400)
# ## OOC & scenario 1             :   rdirichlet(1, c(5/10-delta/2, 4/10-delta/2, 1/10, delta)*400)
# ## OOC & scenario 2             :   rdirichlet(1, c(5.5/10, 4.5/10, 0, 0)*400)
# ## OOC & scenario 3             :   rdirichlet(1, c(5/10-delta/2, 4/10-delta/2, 1/10+delta, 0)*400)
# if(isPART0 == F)
# {
#   set.seed(seed_num)
# 
#   # parameters for dirichlet distribution
#   if(scenario == 1)
#   {
#     tmp_param_dirichlet1 <- c(5/10, 4/10, 1/10, 0)*400
#     tmp_param_dirichlet2 <- c(5/10-delta/2, 4/10-delta/2, 1/10, delta)*400
#   }else if(scenario == 2)
#   {
#     tmp_param_dirichlet1 <- c(5.5/10-delta/2, 4.5/10-delta/2, delta, 0)*400
#     tmp_param_dirichlet2 <- c(5.5/10, 4.5/10, 0, 0)*400
#   }else if(scenario == 3)
#   {
#     tmp_param_dirichlet1 <- c(5/10, 4/10, 1/10, 0)*400
#     tmp_param_dirichlet2 <- c(5/10-delta/2, 4/10-delta/2, 1/10+delta, 0)*400
#   }
# 
#   # data generation
#   ## probability
#   Prob_matrix <- matrix(NA, nrow=Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS, ncol=4)
#   for(i in seq(Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS))
#   {
#     if(i <= Phase1_SS + Phase2_IC_SS)
#     {
#       Prob_matrix[i,] <- rdirichlet(1, tmp_param_dirichlet1)
#     }else{
#       Prob_matrix[i,] <- rdirichlet(1, tmp_param_dirichlet2)
#     }
#   }
# 
#   ## cluster structure
#   cluster_matrix <- matrix(NA, nrow=Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS, ncol=4)
#   for(i in seq(Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS))
#   {
#     cluster_matrix[i,] <- c(rmultinom(1, Subgroup_SS, Prob_matrix[i,]))
#   }
# 
#   ## data generation
#   sig <- 1
#   DATAlist <- list()
#   for(i in seq(Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS))
#   {
#     tmp_cluster <- sapply(1:4, function(x){ rep(x,cluster_matrix[i,x]) }) %>% unlist()
#     tmp_DATA <- NULL
#     for(j in 1:Subgroup_SS)
#     {
#       tmp_coef1 <- true_coef1 + rnorm(trueK, 0, 0.1)
#       tmp_coef2 <- true_coef2 + rnorm(trueK, 0, 0.1)
#       tmp_coef3 <- true_coef3 + rnorm(trueK, 0, 0.1)
#       tmp_coef4 <- true_coef4 + rnorm(trueK, 0, 0.1)
#       tmp_COEF <- rbind(tmp_coef1, tmp_coef2, tmp_coef3, tmp_coef4)
# 
#       tmp_time <- i
#       tmp_CurveNUM <- j
#       tmp_clusterNUM <- tmp_cluster[j]
#       tmp_y <- c(basis_matrix %*% tmp_COEF[tmp_clusterNUM,]) + rnorm(length(t), 0, sig)
#       tmp_data <- tibble(time=tmp_time, CurveNUM=tmp_CurveNUM, clusterNUM=tmp_clusterNUM, domain=t, y=tmp_y)
# 
#       ### save data
#       tmp_DATA <- tmp_DATA %>% bind_rows(tmp_data)
#     }
#     DATAlist[[i]] <- tmp_DATA
#   }
#   DATA <- do.call(bind_rows, DATAlist)
# }
# # -------------------------------------------------------------------------------------------------
# #
# ###################################################################################################


# ###################################################################################################
# # plot for checking data
# # -------------------------------------------------------------------------------------------------
# plot(100, 100, xlim=c(0,1), ylim=c(-6,12), xlab='t', ylab="y", main="Generated data for Simulation 1",
#      cex.axis=1.5, cex.lab=2, cex.main=2)
# sig <- 1
# for(i in 1:4)
# {
#   for(j in 1:50)
#   {
#     tmp_coef <- true_COEF[i,] + rnorm(trueK, 0, 0.1)
#     tmp_y <- c(basis_matrix %*% tmp_coef) + rnorm(length(t), 0, sig)
#     lines(t, tmp_y, lty=1, lwd=1, col=i)
#   }
# }
# legend("topright", c("cluster 1", "cluster 2", "cluster 3", "cluster 4"), lty=1, col=1:4, lwd=2)
# # -------------------------------------------------------------------------------------------------
# #
# ###################################################################################################