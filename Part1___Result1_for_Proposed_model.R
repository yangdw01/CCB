###################################################################################################
###################################################################################################
####
####
####      JQT         by DY
####
####      (Simulation 2) - PART1. Result for the proposed model
####
####
###################################################################################################
###################################################################################################
#
# Code explanations :
# 
#   Obtain the result for the proposed model 
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
library(ggplot2)
library(gridExtra)
library(plotly)
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
# data generation
# -------------------------------------------------------------------------------------------------
isPART0 <- F
seed_num <- 2
scenario <- 1 # 2 # 3 #
delta <- 0.15 # 0.05, 0.1, 0.15, 0.2
L <- 3
DATA <- DataGen(isPART0=isPART0, seednum=seed_num, delta=delta, scenario=scenario)
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



###################################################################################################
# (Step 2) SUGS algorithm
# -------------------------------------------------------------------------------------------------
# basic setting
Y = Step1_result
M <- 15
SUGS_ordering_num <- 10

# priors
a = 1
m = 0
psi = 1
c = 1
d = 10
alpha_stars <- c(0.01, 0.05, c(0.1 + 0.2 * (0:19)))
eta_probs <- rep(1/length(alpha_stars), length(alpha_stars))

# target time
target_time <- c(Phase1_SS, Phase1_SS + Phase2_IC_SS/2, # in-control
                 Phase1_SS + Phase2_IC_SS+1, Phase1_SS + Phase2_IC_SS+Phase2_OOC_SS) # out-of-control
target_parameters <- list()

# online learning
Final_ordering <- NULL
multinomial_prob <- NULL
tmp_SSlist <- rep(Subgroup_SS, Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS)
for(w in 1:(Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS)){
  
  tmp_idx <- ((w-1) * Subgroup_SS) + 1:Subgroup_SS
  
  tmp_result_list <- list()
  tmp_ordering_list <- NULL 
  
  for(ww in 1:SUGS_ordering_num){
    
    cat(w, "th /", Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS, " : ", ww, "th ordering \n")
    
    # ordering
    if(ww == 1){
      
      tmp_resample <- seq(tmp_idx)
      
    }else{
      
      set.seed(ww + w*(Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS + SUGS_ordering_num))  
      tmp_resample <- sample(seq(tmp_idx), length(tmp_idx))
    }
    tmp_Y <- Y[tmp_idx[tmp_resample],]
    
    # SUGS algorithm
    if(w==1){
      
      tmp_result <- SUGS_algorithm_sequential_fast(Y=tmp_Y, a=a, m=m, psi=psi, c=c, d=d, alpha_stars=alpha_stars, eta_probs=eta_probs, first_iter=TRUE,
                                                   prev_bhat_j=NULL, prev_Gamma_i=NULL, prev_phi_t=NULL, prev_psi_hj=NULL, prev_m_hj=NULL, prev_a_hj=NULL, prev_b_hj=NULL)
    }else{
      
      tmp_result <- SUGS_algorithm_sequential_fast(Y=tmp_Y, a=a, m=m, psi=psi, c=c, d=d, alpha_stars=alpha_stars, eta_probs=eta_probs, first_iter=FALSE,
                                                   prev_bhat_j=tmp_result__bhat_j, prev_Gamma_i=tmp_result__Gamma_i, prev_phi_t=tmp_result__phi_t, 
                                                   prev_psi_hj=tmp_result__psi_hj, prev_m_hj=tmp_result__m_hj, prev_a_hj=tmp_result__a_hj, prev_b_hj=tmp_result__b_hj)
    }
    
    tmp_result_list[[ww]] <- tmp_result
    tmp_ordering_list <- rbind(tmp_ordering_list, tmp_idx[tmp_resample])
  }
  max_idx <- sapply(1:SUGS_ordering_num, function(x){ tmp_result_list[[x]][[1]] }) %>% which.max()
  
  # SUGS algorithm result
  tmp_result <- tmp_result_list[[max_idx]]
  tmp_result__bhat_j   <- tmp_result[[2]]
  tmp_result__Gamma_i  <- tmp_result[[3]]
  tmp_result__phi_t    <- tmp_result[[4]]
  tmp_result__psi_hj   <- tmp_result[[5]]
  tmp_result__m_hj     <- tmp_result[[6]]
  tmp_result__a_hj     <- tmp_result[[7]]
  tmp_result__b_hj     <- tmp_result[[8]]
  
  # for fast calculation
  if(w > Phase1_SS + 1)
  {
    delete_idx <- (1:Subgroup_SS) + Phase1_SS*Subgroup_SS
    
    tmp_result__bhat_j  <- tmp_result__bhat_j[-delete_idx,]
    tmp_result__phi_t   <- tmp_result__phi_t[-delete_idx,]
    tmp_result__psi_hj  <- sapply(seq(tmp_result__psi_hj), function(x){ tmp_result__psi_hj[[x]][-delete_idx,] }, simplify = F)
    tmp_result__m_hj    <- sapply(seq(tmp_result__m_hj),   function(x){   tmp_result__m_hj[[x]][-delete_idx,] }, simplify = F)
    tmp_result__a_hj    <- sapply(seq(tmp_result__a_hj),   function(x){   tmp_result__a_hj[[x]][-delete_idx,] }, simplify = F)
    tmp_result__b_hj    <- sapply(seq(tmp_result__b_hj),   function(x){   tmp_result__b_hj[[x]][-delete_idx,] }, simplify = F)
  }
  
  # ordering
  Final_ordering <- c(Final_ordering, tmp_ordering_list[max_idx,])
  
  # save for cluster-specific mean curve or surface
  if(w %in% target_time){
    
    tmp_idx <- which(target_time %in% w)
    
    target_parameters[[tmp_idx]] <- list( sapply(tmp_result__m_hj,   function(x){ x %>% tail(1) %>% c }),
                                          sapply(tmp_result__psi_hj, function(x){ x %>% tail(1) %>% c }),
                                          sapply(tmp_result__a_hj,   function(x){ x %>% tail(1) %>% c }),
                                          sapply(tmp_result__b_hj,   function(x){ x %>% tail(1) %>% c })   )
  }
}

# SUGS algorithm result
Final_result___bhat_j  <- tmp_result__bhat_j
Final_result___Gamma_i <- tmp_result__Gamma_i
Final_result___phi_t   <- tmp_result__phi_t
Final_result___psi_hj  <- tmp_result__psi_hj
Final_result___m_hj    <- tmp_result__m_hj
Final_result___a_hj    <- tmp_result__a_hj
Final_result___b_hj    <- tmp_result__b_hj
Gamma_i <- Final_result___Gamma_i

# clustering result
true_clusters <- DATA %>% distinct(time, CurveNUM, clusterNUM) %>% pull(clusterNUM)
rd_result  <- rand.index(Gamma_i, true_clusters[Final_ordering])
ard_result <- adj.rand.index(Gamma_i, true_clusters[Final_ordering])
# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# (Step 3) p-chart
# -------------------------------------------------------------------------------------------------
# p chart based on Bayesian approach
## parameter estimation
Phase1_idx <- tmp_SSlist[1:Phase1_SS] %>% sum %>% seq
Phase1_alpha_prob <- (Final_result___phi_t[Phase1_idx,] %>% tail(1) %>% c())
Phase1_alpha <- ((Final_result___phi_t[Phase1_idx,] %>% tail(1) %>% c()) * alpha_stars) %>% sum()
Phase1_Gamma_i <- Final_result___Gamma_i[Phase1_idx]
Phase1_Gamma_params <- sapply(1:M, function(k){ c(sum(Phase1_Gamma_i == k),sum(Phase1_Gamma_i > k)) })
Phase1_pks <- sapply(1:M, function(k){ rbeta(1, 1+Phase1_Gamma_params[1,k], Phase1_Gamma_params[2,k]+Phase1_alpha) })
Phase1_piks <- sapply(1:M, function(k){ Phase1_pks[k] * prod(c(1,1-Phase1_pks)[1:k]) }) %>% unlist()
Phase1_piks <- c(Phase1_piks, 1-sum(Phase1_piks))

## plotting statistics
multinomial_result <- NULL
for(w in 1:(Phase1_SS + Phase2_IC_SS + Phase2_OOC_SS)){
  
  tmp_idx <- ((w-1) * Subgroup_SS) + 1:Subgroup_SS
  tmp_Gamma <- Final_result___Gamma_i[tmp_idx]
  tmp_multinomial_result <- sapply(1:(M+1), function(k){ sum(tmp_Gamma == k)/Subgroup_SS })
  
  multinomial_result <- rbind(multinomial_result, tmp_multinomial_result)
}

## Control limits
Ls <-  c(3.38, 3.85, 3.79, 4.60, 3.63)
Centerline <- Phase1_piks
CLs <- sapply(1:3, function(j){ calculate_CL(probs=Centerline, tmpL=Ls[j], tmpN=Subgroup_SS, case=j) }, simplify = F) 

## plot
j <- 3

Titleneams <- c("(a)", "(b)", "(c)", "(d)")
Figurenames <- paste0("Fig7",c("a","b","c","d"),".png")

plotting_clusters <- 1:4
PSs <- multinomial_result[,plotting_clusters]
CENTERs <- Centerline[plotting_clusters]
CONTROLLs <- sapply(CLs, function(x){ x[,plotting_clusters] }, simplify = F)
yvals <- c(0.70, 0.585, 0.24, 0.24)

cluster <- 1
tmp_df <- tibble(t = seq(PSs[,cluster]), 
                 statistic = PSs[,cluster], 
                 CenterLine = rep(CENTERs[cluster],   length(PSs[,cluster])), 
                 LCL = rep(CONTROLLs[[j]][1,cluster], length(PSs[,cluster])), 
                 UCL = rep(CONTROLLs[[j]][2,cluster], length(PSs[,cluster])))
tmp_df %>% 
  ggplot() +
  geom_line(aes(x=t, y=statistic),                  linewidth=1.5) + 
  geom_point(aes(x=t, y=statistic),                 size=3) + 
  geom_line(aes(x=t, y=LCL, colour="Control Limit"), linewidth=1.5, linetype=2) +
  geom_line(aes(x=t, y=UCL, colour="Control Limit"), linewidth=1.5, linetype=2) +
  geom_point(data = tmp_df %>% slice(target_time), aes(x=t, y=statistic), size=4, color=5) +
  geom_vline(aes(xintercept = Phase1_SS + 1/2,              ), linewidth=1.5, linetype=3, colour=3) +
  geom_vline(aes(xintercept = Phase1_SS + Phase2_IC_SS + 1/2), linewidth=1.5, linetype=3, colour=4) + 
  theme_bw() +
  labs(title=Titleneams[cluster], y = "Cluster proportion", color = "") +
  theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
        legend.title= element_text(size=27),
        legend.text = element_text(size=12),
        legend.box="horizontal",
        axis.text   = element_text(size=20),
        axis.title  = element_text(size=20,face="bold"),
        legend.key=element_rect(fill=NA)) +
  geom_text(aes(x=Phase1_SS + 1/2-4, y= yvals[cluster], label = "Phase I"),    colour = 3, vjust = 1, size = 6) +
  geom_text(aes(x=Phase1_SS + Phase2_IC_SS + 1/2-4, y= yvals[cluster], label = "In-control"), colour = 4, vjust = 1, size = 6) +
  geom_text(aes(x=50, y= yvals[cluster], label = "Out-of-control"), colour = 1, vjust = 1, size = 6)



figs <- list()
cluster <- 1
tmp_df <- tibble(t = seq(PSs[,cluster]), 
                 statistic = PSs[,cluster], 
                 CenterLine = rep(CENTERs[cluster],   length(PSs[,cluster])), 
                 LCL = rep(CONTROLLs[[j]][1,cluster], length(PSs[,cluster])), 
                 UCL = rep(CONTROLLs[[j]][2,cluster], length(PSs[,cluster])))
figs[[1]] <- tmp_df %>% 
  ggplot() +
  geom_line(aes(x=t, y=statistic),                  linewidth=1.5) + 
  geom_point(aes(x=t, y=statistic),                 size=3) + 
  geom_line(aes(x=t, y=LCL, colour="Control Limit"), linewidth=1.5, linetype=2) +
  geom_line(aes(x=t, y=UCL, colour="Control Limit"), linewidth=1.5, linetype=2) +
  geom_point(data = tmp_df %>% slice(target_time), aes(x=t, y=statistic), size=4, color=5) +
  geom_vline(aes(xintercept = Phase1_SS + 1/2,              ), linewidth=1.5, linetype=3, colour=3) +
  geom_vline(aes(xintercept = Phase1_SS + Phase2_IC_SS + 1/2), linewidth=1.5, linetype=3, colour=4) + 
  theme_bw() +
  labs(title=Titleneams[cluster], y = "Cluster proportion", color = "") +
  theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
        legend.title= element_text(size=27),
        legend.text = element_text(size=16),
        legend.box="horizontal",
        axis.text   = element_text(size=20),
        axis.title  = element_text(size=20,face="bold"),
        legend.key=element_rect(fill=NA)) +
  geom_text(aes(x=Phase1_SS + 1/2-4,                y= yvals[1], label = "Phase I"),    colour = 3, vjust = 1, size = 6) +
  geom_text(aes(x=Phase1_SS + Phase2_IC_SS + 1/2-4, y= yvals[1], label = "In-control"), colour = 4, vjust = 1, size = 6) +
  geom_text(aes(x=50,                               y= yvals[1], label = "Out-of-control"), colour = 1, vjust = 1, size = 6)

cluster <- 2
tmp_df <- tibble(t = seq(PSs[,cluster]), 
                 statistic = PSs[,cluster], 
                 CenterLine = rep(CENTERs[cluster],   length(PSs[,cluster])), 
                 LCL = rep(CONTROLLs[[j]][1,cluster], length(PSs[,cluster])), 
                 UCL = rep(CONTROLLs[[j]][2,cluster], length(PSs[,cluster])))
figs[[2]] <- tmp_df %>% 
  ggplot() +
  geom_line(aes(x=t, y=statistic),                  linewidth=1.5) + 
  geom_point(aes(x=t, y=statistic),                 size=3) + 
  geom_line(aes(x=t, y=LCL, colour="Control Limit"), linewidth=1.5, linetype=2) +
  geom_line(aes(x=t, y=UCL, colour="Control Limit"), linewidth=1.5, linetype=2) +
  geom_point(data = tmp_df %>% slice(target_time), aes(x=t, y=statistic), size=4, color=5) +
  geom_vline(aes(xintercept = Phase1_SS + 1/2,              ), linewidth=1.5, linetype=3, colour=3) +
  geom_vline(aes(xintercept = Phase1_SS + Phase2_IC_SS + 1/2), linewidth=1.5, linetype=3, colour=4) + 
  theme_bw() +
  labs(title=Titleneams[cluster], y = "Cluster proportion", color = "") +
  theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
        legend.title= element_text(size=27),
        legend.text = element_text(size=16),
        legend.box="horizontal",
        axis.text   = element_text(size=20),
        axis.title  = element_text(size=20,face="bold"),
        legend.key=element_rect(fill=NA)) +
  geom_text(aes(x=Phase1_SS + 1/2-4,                y= yvals[2], label = "Phase I"),    colour = 3, vjust = 1, size = 6) +
  geom_text(aes(x=Phase1_SS + Phase2_IC_SS + 1/2-4, y= yvals[2], label = "In-control"), colour = 4, vjust = 1, size = 6) +
  geom_text(aes(x=50,                               y= yvals[2], label = "Out-of-control"), colour = 1, vjust = 1, size = 6)

cluster <- 3
tmp_df <- tibble(t = seq(PSs[,cluster]), 
                 statistic = PSs[,cluster], 
                 CenterLine = rep(CENTERs[cluster],   length(PSs[,cluster])), 
                 LCL = rep(CONTROLLs[[j]][1,cluster], length(PSs[,cluster])), 
                 UCL = rep(CONTROLLs[[j]][2,cluster], length(PSs[,cluster])))
figs[[3]] <- tmp_df %>% 
  ggplot() +
  geom_line(aes(x=t, y=statistic),                  linewidth=1.5) + 
  geom_point(aes(x=t, y=statistic),                 size=3) + 
  geom_line(aes(x=t, y=LCL, colour="Control Limit"), linewidth=1.5, linetype=2) +
  geom_line(aes(x=t, y=UCL, colour="Control Limit"), linewidth=1.5, linetype=2) +
  geom_point(data = tmp_df %>% slice(target_time), aes(x=t, y=statistic), size=4, color=5) +
  geom_vline(aes(xintercept = Phase1_SS + 1/2,              ), linewidth=1.5, linetype=3, colour=3) +
  geom_vline(aes(xintercept = Phase1_SS + Phase2_IC_SS + 1/2), linewidth=1.5, linetype=3, colour=4) + 
  theme_bw() +
  labs(title=Titleneams[cluster], y = "Cluster proportion", color = "") +
  theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
        legend.title= element_text(size=27),
        legend.text = element_text(size=16),
        legend.box="horizontal",
        axis.text   = element_text(size=20),
        axis.title  = element_text(size=20,face="bold"),
        legend.key=element_rect(fill=NA)) +
  geom_text(aes(x=Phase1_SS + 1/2-4,                y= yvals[3], label = "Phase I"),    colour = 3, vjust = 1, size = 6) +
  geom_text(aes(x=Phase1_SS + Phase2_IC_SS + 1/2-4, y= yvals[3], label = "In-control"), colour = 4, vjust = 1, size = 6) +
  geom_text(aes(x=50,                               y= yvals[3], label = "Out-of-control"), colour = 1, vjust = 1, size = 6)

cluster <- 4
tmp_df <- tibble(t = seq(PSs[,cluster]), 
                 statistic = PSs[,cluster], 
                 CenterLine = rep(CENTERs[cluster],   length(PSs[,cluster])), 
                 LCL = rep(CONTROLLs[[j]][1,cluster], length(PSs[,cluster])), 
                 UCL = rep(CONTROLLs[[j]][2,cluster], length(PSs[,cluster])))
figs[[4]] <- tmp_df %>% 
  ggplot() +
  geom_line(aes(x=t, y=statistic),                  linewidth=1.5) + 
  geom_point(aes(x=t, y=statistic),                 size=3) + 
  geom_line(aes(x=t, y=LCL, colour="Control Limit"), linewidth=1.5, linetype=2) +
  geom_line(aes(x=t, y=UCL, colour="Control Limit"), linewidth=1.5, linetype=2) +
  geom_point(data = tmp_df %>% slice(target_time), aes(x=t, y=statistic), size=4, color=5) +
  geom_vline(aes(xintercept = Phase1_SS + 1/2,              ), linewidth=1.5, linetype=3, colour=3) +
  geom_vline(aes(xintercept = Phase1_SS + Phase2_IC_SS + 1/2), linewidth=1.5, linetype=3, colour=4) + 
  theme_bw() +
  labs(title=Titleneams[cluster], y = "Cluster proportion", color = "") +
  theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
        legend.title= element_text(size=27),
        legend.text = element_text(size=16),
        legend.box="horizontal",
        axis.text   = element_text(size=20),
        axis.title  = element_text(size=20,face="bold"),
        legend.key=element_rect(fill=NA)) +
  geom_text(aes(x=Phase1_SS + 1/2-4,                y= yvals[4], label = "Phase I"),    colour = 3, vjust = 1, size = 6) +
  geom_text(aes(x=Phase1_SS + Phase2_IC_SS + 1/2-4, y= yvals[4], label = "In-control"), colour = 4, vjust = 1, size = 6) +
  geom_text(aes(x=50,                               y= yvals[4], label = "Out-of-control"), colour = 1, vjust = 1, size = 6)

png("Fig9.png", width = 600, height = 300, units='mm', res = 900)
myplot <- grid.arrange(
  grobs = figs,
  widths = c(1,1),
  layout_matrix = rbind(c(1,2), c(3,4))
)
print(myplot)
dev.off()





###### Monte Carlo Simulation 
I <- 10000

basis <- basis_matrix

df_plot <- NULL
for(w in 1:length(target_time))
{
  tmp_m   <- target_parameters[[w]][[1]]
  tmp_psi <- target_parameters[[w]][[2]]
  tmp_a   <- target_parameters[[w]][[3]]
  tmp_b   <- target_parameters[[w]][[4]]
  
  tmp_clusters <- 1:ncol(tmp_m)
  tmp_mean_function <- NULL
  for(cluster in tmp_clusters){
    
    tmp_MCsample <- NULL
    for(i in 1:I){
      
      tmp_tau <- rgamma(nrow(tmp_a), shape=tmp_a[,cluster], rate=tmp_b[,cluster])
      tmp_mu  <- rnorm(nrow(tmp_m), tmp_m[,cluster], sqrt(tmp_psi[,cluster]/tmp_tau))
      tmp_MCsample <- rbind(tmp_MCsample, tmp_mu)
    }
    
    tmp_pararmeter <- tmp_MCsample %>% apply(2,mean)
    tmp_mean_function <- rbind(tmp_mean_function, c(basis %*% tmp_pararmeter))
  }
  
  for(cluster in tmp_clusters){
    
    tmp_tibble <- tibble(Time=w, Cluster=cluster, v1=c(t1), v2=c(t2), f=tmp_mean_function[cluster,])
    df_plot <- df_plot %>% bind_rows(tmp_tibble)
  }
}
df_plot <- df_plot %>% mutate(Cluster = as.factor(Cluster))

Titleneams <- c("<b>(a) t=20<b>", "<b>(b) t=20<b>", "<b>(c) t=20<b>", 
                "<b>(d) t=25<b>", "<b>(e) t=25<b>", "<b>(f) t=25<b>",
                "<b>(g) t=31<b>", "<b>(h) t=31<b>", "<b>(i) t=31<b>", "<b>(j) t=31<b>",
                "<b>(k) t=70<b>", "<b>(l) t=70<b>", "<b>(m) t=70<b>", "<b>(n) t=70<b>")
w <- 4
cc <- 4
tmp_tmp_df <- df_plot %>% filter(Time == w) %>% filter(Cluster == cc) %>% pivot_wider(names_from = v2, values_from = f) %>% select(-v1, -Time, -Cluster) %>% as.matrix

plot_ly() %>% add_surface(x=t1, y=t2, z=~tmp_tmp_df, zmin=-5,zmax=5,
                          colorbar=list(title='', tickfont=list(size=20), x=1, y=0.8)) %>%
  layout(scene = list(aspectmode = "manual",
                      aspectratio = list(x=1.1, y=1.1, z=1.1),
                      zaxis=list(range=c(-5,5), title="f(v1,v2)", titlefont = list(size = 30), tickfont = list(size = 14)),
                      xaxis=list(title="v1", titlefont = list(size = 22), tickfont = list(size = 14)),
                      yaxis=list(title="v2", titlefont = list(size = 22), tickfont = list(size = 14))),
         title = Titleneams[14], titlefont = list(size = 50, face="bold"),
         margin = list(t=100))




for(w in 1:length(target_time)){
  
  tmp_df <- df_plot %>% filter(Time == w)
  tmp_cluster_num <- tmp_df$Cluster %>% unique %>% length

  
  myplot <- tmp_df %>%
    ggplot( aes(x=t, y=f, group = Cluster, color=Cluster)) +
    geom_line(linewidth=2) + 
    theme_bw() +
    labs(title=Titleneams[w]) +
    scale_x_continuous(name = "v") +
    scale_y_continuous(name = "f(v)") +
    theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
          legend.title= element_text(size=20),
          legend.text = element_text(size=20),
          axis.text   = element_text(size=20),
          axis.title  = element_text(size=25,face="bold"))
  figs2[[w]] <- myplot
}

png("Fig8.png", width = 600, height = 300, units='mm', res = 900)
myplot <- grid.arrange(
  grobs = figs2,
  widths = c(1,1),
  layout_matrix = rbind(c(1,2), c(3,4))
)
print(myplot)
dev.off()