# library(dplyr)
# 
# scenario <- 4
# 
# for(delta in c(0.05,0.10,0.15,0.20))
# {
#   load("Part0_result//PART0result_500_freq.RData")
#   part0_output2 <- output2
#   
#   load(file=paste("Part2_result//PART2_scenario", scenario, "delta", delta, "_result_freq.RData", sep=""))
#   part2_output2 <- output2
#   
#   for(i in 1:length(part2_output2))
#   {
#     tmpi <- sample(1:length(part0_output2), 1)
#     part2_output2[[i]]$Frequentist$PS <- part0_output2[[tmpi]]$PS
#     part2_output2[[i]]$Frequentist$centerline <- part0_output2[[tmpi]]$centerline
#     part2_output2[[i]]$Frequentist$sigma <- part0_output2[[tmpi]]$sigma
#   }
#   output2 <- part2_output2
#   
#   save(output2, file=paste("PART2_scenario", scenario, "delta", delta, "_result_freq.RData", sep=""))
# }



###################################################################################################
# calculate L for ARL0 = 100 (Wilson score confidence interval)
# -------------------------------------------------------------------------------------------------
library(abind)
library(ggplot2)
source("DataGen.R")
Ls <-  c(3.41, 3.82, 3.78, 4.55, 3.63)
Ls2 <- c(3.42, 3.69, 3.65, 4.23, 3.51)
# -------------------------------------------------------------------------------------------------

# functions ---------------------------------------------------------------------------------------
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
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
ALLresult_ARL1 <- ALLresult2_ARL1 <- ALLresult_SDRL <- ALLresult2_SDRL <- ALLresult_clustering <- ALLresult2_clustering <- NULL
for(scenario in c(1,2,3,4)){
  
  tmp_ARL1_result1 <- tmp_ARL1_result2 <- tmp_SDRL_result1 <- tmp_SDRL_result2 <- tmp_Clustering_result1 <- tmp_Clustering_result2 <- NULL
  for(delta in c(0.05, 0.1, 0.15, 0.2)){
    
    # calcualte ARL1 ----------------------------------------------------------------------------------
    load(file=paste("Part2_result//PART2_scenario", scenario, "delta", delta, "_result.RData", sep=""))
    load(file=paste("Part2_result//PART2_scenario", scenario, "delta", delta, "_result_freq.RData", sep=""))
    
    ARL1_result <- ARL1_result2 <- NULL
    ARL1_mat_result <- ARL1_mat_result2 <- rep(list(NULL), 5)
    for(i in seq(output))
    {
      tmp_output  <- output[[i]]
      tmp_output2 <- output2[[i]]
      
      # baeysian 
      tmp_CLs <-  sapply(1:5, function(j){ calculate_CL(probs=tmp_output$Bayesian[[1]],    tmpL=Ls[j],  tmpN=Subgroup_SS, case=j) }, simplify = F) 
      
      # frequentist
      tmp_prob <- (tmp_output2$Frequentist[[1]] %>% sort(decreasing = T))
      tmp_CLs2 <- sapply(1:5, function(j){ calculate_CL(probs=tmp_prob,                    tmpL=Ls2[j], tmpN=Subgroup_SS, case=j) }, simplify = F) 
      
      # alarm check
      tmp_alarm_chk <- tmp_alarm_chk2 <- tmp_RLmat <- tmp_RLmat2 <- list()
      for(j in seq(Ls)){
        
        if(scenario == 4){ tmp_Phase2_IC_SS <- Phase2_IC_SS + 5 }else{ tmp_Phase2_IC_SS <- Phase2_IC_SS }
        
        ## Bayesian
        tmp_LCL_mat         <- matrix(tmp_CLs[[j]][1,], nrow=nrow(tmp_output$Bayesian[[3]]), ncol=length(tmp_output$Bayesian[[1]]), byrow=T)
        tmp_UCL_mat         <- matrix(tmp_CLs[[j]][2,], nrow=nrow(tmp_output$Bayesian[[3]]), ncol=length(tmp_output$Bayesian[[1]]), byrow=T)
        tmp_alarm_chk[[j]]  <- (tmp_LCL_mat > tmp_output$Bayesian[[3]]) | (tmp_UCL_mat < tmp_output$Bayesian[[3]])
        tmp_RLmat[[j]]      <- tmp_alarm_chk[[j]][-(1:tmp_Phase2_IC_SS),]  %>% apply(2, function(x){ ifelse(length(which(x))>0, which(x)[1], nrow(tmp_alarm_chk[[j]]) - tmp_Phase2_IC_SS) })
        
        ## Frequentist
        tmp_LCL_mat2        <- matrix(tmp_CLs2[[j]][1,], nrow=nrow(tmp_output2$Frequentist[[3]]), ncol=length(tmp_output2$Frequentist[[1]]), byrow=T)
        tmp_UCL_mat2        <- matrix(tmp_CLs2[[j]][2,], nrow=nrow(tmp_output2$Frequentist[[3]]), ncol=length(tmp_output2$Frequentist[[1]]), byrow=T)
        tmp_alarm_chk2[[j]] <- (tmp_LCL_mat2 > tmp_output2$Frequentist[[3]]) | (tmp_UCL_mat2 < tmp_output2$Frequentist[[3]])
        tmp_RLmat2[[j]]     <- tmp_alarm_chk2[[j]][-(1:Phase2_IC_SS),] %>% apply(2, function(x){ ifelse(length(which(x))>0, which(x)[1], nrow(tmp_alarm_chk2[[j]]) - Phase2_IC_SS) })
      }
      
      # save
      for(j in seq(Ls)){
        
        ARL1_mat_result[[j]]  <- rbind(ARL1_mat_result[[j]],  tmp_RLmat[[j]])
        ARL1_mat_result2[[j]] <- rbind(ARL1_mat_result2[[j]], tmp_RLmat2[[j]])
      }
    }
    
    tmp_ARL1_result1 <- rbind(tmp_ARL1_result1, sapply(1:length(Ls),  function(j){ ARL1_mat_result[[j]]  %>% apply(1,min) }) %>% apply(2, mean))
    tmp_ARL1_result2 <- rbind(tmp_ARL1_result2, sapply(1:length(Ls2), function(j){ ARL1_mat_result2[[j]] %>% apply(1,min) }) %>% apply(2, mean))
    
    tmp_SDRL_result1 <- rbind(tmp_SDRL_result1, sapply(1:length(Ls),  function(j){ ARL1_mat_result[[j]]  %>% apply(1,min) }) %>% apply(2, sd))
    tmp_SDRL_result2 <- rbind(tmp_SDRL_result2, sapply(1:length(Ls2), function(j){ ARL1_mat_result2[[j]] %>% apply(1,min) }) %>% apply(2, sd))
    
    # cluster result
    Bayesian_clustering    <- sapply(1:length(output),  function(x){ output[[x]]$Clustering[1] })
    Frequentist_clustering <- sapply(1:length(output2), function(x){ output2[[x]]$Clustering[1] })
    
    tmp_Clustering_result1 <- rbind(tmp_Clustering_result1, Bayesian_clustering)
    tmp_Clustering_result2 <- rbind(tmp_Clustering_result2, Frequentist_clustering)
    # ------------------------------------------------------------------------------------------------- 
  }
  
  ALLresult_ARL1[[scenario]]  <- tmp_ARL1_result1
  ALLresult2_ARL1[[scenario]] <- tmp_ARL1_result2
  
  ALLresult_SDRL[[scenario]]  <- tmp_SDRL_result1
  ALLresult2_SDRL[[scenario]] <- tmp_SDRL_result2 
  
  ALLresult_clustering[[scenario]]  <- tmp_Clustering_result1
  ALLresult2_clustering[[scenario]] <- tmp_Clustering_result2
}
# -------------------------------------------------------------------------------------------------

par(mfrow=c(3,4))
j <- 1
for(scenario in 1:4){
  plot(ALLresult_ARL1[[scenario]][,j], type='b', ylim=quantile(c(ALLresult_ARL1[[scenario]][,j],ALLresult2_ARL1[[scenario]][,j]), c(0,1)), lwd=2)
  lines(ALLresult2_ARL1[[scenario]][,j], type='b', col=2, lwd=2)
}

j <- 2
for(scenario in 1:4){
  plot(ALLresult_ARL1[[scenario]][,j], type='b', ylim=quantile(c(ALLresult_ARL1[[scenario]][,j],ALLresult2_ARL1[[scenario]][,j]), c(0,1)), lwd=2)
  lines(ALLresult2_ARL1[[scenario]][,j], type='b', col=2, lwd=2)
}

j <- 3
for(scenario in 1:4){
  plot(ALLresult_ARL1[[scenario]][,j], type='b', ylim=quantile(c(ALLresult_ARL1[[scenario]][,j],ALLresult2_ARL1[[scenario]][,j]), c(0,1)), lwd=2)
  lines(ALLresult2_ARL1[[scenario]][,j], type='b', col=2, lwd=2)
}



par(mfrow=c(3,4))
j <- 1
for(scenario in 1:4){
  plot(ALLresult_SDRL[[scenario]][,j], type='b', 
       ylim=quantile(c(ALLresult_SDRL[[scenario]][,j],ALLresult2_SDRL[[scenario]][,j]), c(0,1)), lwd=2)
  lines(ALLresult2_SDRL[[scenario]][,j], type='b', col=2, lwd=2)
}

j <- 2
for(scenario in 1:4){
  plot(ALLresult_SDRL[[scenario]][,j], type='b', ylim=quantile(c(ALLresult_SDRL[[scenario]][,j],ALLresult2_SDRL[[scenario]][,j]), c(0,1)), lwd=2)
  lines(ALLresult2_SDRL[[scenario]][,j], type='b', col=2, lwd=2)
}

j <- 3
for(scenario in 1:4){
  plot(ALLresult_SDRL[[scenario]][,j], type='b', ylim=quantile(c(ALLresult_SDRL[[scenario]][,j],ALLresult2_SDRL[[scenario]][,j]), c(0,1)), lwd=2)
  lines(ALLresult2_SDRL[[scenario]][,j], type='b', col=2, lwd=2)
}

for(scenario in 1:4){ rbind(ALLresult_clustering[[scenario]], ALLresult2_clustering[[scenario]]) %>% t %>% boxplot }


# -------------------------------------------------------------------------------------------------
#
###################################################################################################


# tibble ------------------------------------------------------------------------------------------
ARL1_result <- NULL
for(scenario in 1:4)
{
  # Bayesian
  tmp1 <- ALLresult_ARL1[[scenario]] %>% t %>% as.data.frame  
  colnames(tmp1) <- paste("d", 1:4, sep="")
  tmp1 <- tmp1 %>% tibble %>% mutate(Interval = 1:5) %>% 
    pivot_longer(cols=d1:d4, names_to='Delta', values_to='ARL1') %>% 
    mutate(Scenario = scenario, Method = "CCB")
  
  # Frequentist
  tmp2 <- ALLresult2_ARL1[[scenario]] %>% t %>% as.data.frame  
  colnames(tmp2) <- paste("d", 1:4, sep="")
  tmp2 <- tmp2 %>% tibble %>% mutate(Interval = 1:5) %>% 
    pivot_longer(cols=d1:d4, names_to='Delta', values_to='ARL1') %>% 
    mutate(Scenario = scenario, Method = "CCK")
  
  ARL1_result %>% bind_rows(tmp1) %>% bind_rows(tmp2) -> ARL1_result
}

SDRL_result <- NULL
for(scenario in 1:4)
{
  # Bayesian
  tmp1 <- ALLresult_SDRL[[scenario]] %>% t %>% as.data.frame  
  colnames(tmp1) <- paste("d", 1:4, sep="")
  tmp1 <- tmp1 %>% tibble %>% mutate(Interval = 1:5) %>% 
    pivot_longer(cols=d1:d4, names_to='Delta', values_to='SDRL') %>% 
    mutate(Scenario = scenario, Method = "CCB")
  
  # Frequentist
  tmp2 <- ALLresult2_SDRL[[scenario]] %>% t %>% as.data.frame  
  colnames(tmp2) <- paste("d", 1:4, sep="")
  tmp2 <- tmp2 %>% tibble %>% mutate(Interval = 1:5) %>% 
    pivot_longer(cols=d1:d4, names_to='Delta', values_to='SDRL') %>% 
    mutate(Scenario = scenario, Method = "CCK")
  
  SDRL_result %>% bind_rows(tmp1) %>% bind_rows(tmp2) -> SDRL_result
}
# -------------------------------------------------------------------------------------------------

# Figure 6 ----------------------------------------------------------------------------------------
itv <- 3
Titleneams <- c("(a)", "(b)", "(c)", "(d)")
Figurenames <- paste0("Fig6",c("a","b","c","d"),".png")

for(sc in 1:4)
{
  if(sc!=4)
  {
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- ARL1_result %>% filter(Interval == itv, Scenario == sc) 
    myplot <-
      tmp_result_fig %>%
      ggplot( aes(x=Delta, y=ARL1, group=Method, color=Method)) +
      geom_line(linewidth=1.5) + geom_point(size=3) + theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta), y = "ARL1",
           color = "Method") +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold")) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) #+ 
    #scale_y_continuous(limits = c(0,38))
    print(myplot)
    dev.off()
  }else{
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- ARL1_result %>% filter(Interval == itv, Scenario == sc) %>% 
      pivot_wider(names_from = Method, values_from = ARL1)
    myplot <-tmp_result_fig %>%
      ggplot( aes(x=Delta, group = 1)) +
      geom_line(aes(y=CCB), linewidth=1.5, color="#F8766D") + 
      geom_point(aes(y=CCB), size=3, color="#F8766D") + 
      geom_line(aes(y=CCK/20), linewidth=1.5, color="#00BFC4") + 
      geom_point(aes(y=CCK/20), size=3, color="#00BFC4") +
      theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta)) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) + 
      scale_y_continuous(name = "ARL1 (CCB)", limits=c(0,5.5), sec.axis = sec_axis(~.*20, name="ARL1 (CCK)") ) +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold"),
            axis.text.y.left = element_text(color = "#F8766D"),
            axis.text.y.right = element_text(color = "#00BFC4"),
            axis.title.y.left = element_text(color = "#F8766D"),
            axis.title.y.right = element_text(color = "#00BFC4"),
            plot.margin = margin(t = 5, b = 8, l = 5, r = 60))
    print(myplot)
    dev.off()
  }
}
# -------------------------------------------------------------------------------------------------

# Figure S4 ---------------------------------------------------------------------------------------
itv <- 3
Titleneams <- c("(a)", "(b)", "(c)", "(d)")
Figurenames <- paste0("FigS4",c("a","b","c","d"),".png")

for(sc in 1:4)
{
  if(sc!=4)
  {
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- SDRL_result %>% filter(Interval == itv, Scenario == sc) 
    myplot <-
      tmp_result_fig %>%
      ggplot( aes(x=Delta, y=SDRL, group=Method, color=Method)) +
      geom_line(linewidth=1.5) + geom_point(size=3) + theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta), y = "SDRL",
           color = "Method") +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold")) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) #+ 
    #scale_y_continuous(limits = c(0,38))
    print(myplot)
    dev.off()
  }else{
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- SDRL_result %>% filter(Interval == itv, Scenario == sc) %>% 
      pivot_wider(names_from = Method, values_from = SDRL)
    myplot <-tmp_result_fig %>%
      ggplot( aes(x=Delta, group = 1)) +
      geom_line(aes(y=CCB), linewidth=1.5, color="#F8766D") + 
      geom_point(aes(y=CCB), size=3, color="#F8766D") + 
      geom_line(aes(y=CCK/13), linewidth=1.5, color="#00BFC4") + 
      geom_point(aes(y=CCK/13), size=3, color="#00BFC4") +
      theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta)) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) + 
      scale_y_continuous(name = "SDRL (CCB)", limits=c(0,5.5), sec.axis = sec_axis(~.*13, name="SDRL (CCK)") ) +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold"),
            axis.text.y.left = element_text(color = "#F8766D"),
            axis.text.y.right = element_text(color = "#00BFC4"),
            axis.title.y.left = element_text(color = "#F8766D"),
            axis.title.y.right = element_text(color = "#00BFC4"),
            plot.margin = margin(t = 5, b = 8, l = 5, r = 60))
    print(myplot)
    dev.off()
  }
}
# -------------------------------------------------------------------------------------------------

# Figure S6 ----------------------------------------------------------------------------------------
itv <- 1
Titleneams <- c("(a)", "(b)", "(c)", "(d)")
Figurenames <- paste0("FigS6",c("a","b","c","d"),".png")

for(sc in 1:4)
{
  if(sc!=4)
  {
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- ARL1_result %>% filter(Interval == itv, Scenario == sc) 
    myplot <-
      tmp_result_fig %>%
      ggplot( aes(x=Delta, y=ARL1, group=Method, color=Method)) +
      geom_line(linewidth=1.5) + geom_point(size=3) + theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta), y = "ARL1",
           color = "Method") +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold")) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) #+ 
    #scale_y_continuous(limits = c(0,38))
    print(myplot)
    dev.off()
  }else{
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- ARL1_result %>% filter(Interval == itv, Scenario == sc) %>% 
      pivot_wider(names_from = Method, values_from = ARL1)
    myplot <-tmp_result_fig %>%
      ggplot( aes(x=Delta, group = 1)) +
      geom_line(aes(y=CCB), linewidth=1.5, color="#F8766D") + 
      geom_point(aes(y=CCB), size=3, color="#F8766D") + 
      geom_line(aes(y=CCK/20), linewidth=1.5, color="#00BFC4") + 
      geom_point(aes(y=CCK/20), size=3, color="#00BFC4") +
      theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta)) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) + 
      scale_y_continuous(name = "ARL1 (CCB)", limits=c(0,5.5), sec.axis = sec_axis(~.*20, name="ARL1 (CCK)") ) +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold"),
            axis.text.y.left = element_text(color = "#F8766D"),
            axis.text.y.right = element_text(color = "#00BFC4"),
            axis.title.y.left = element_text(color = "#F8766D"),
            axis.title.y.right = element_text(color = "#00BFC4"),
            plot.margin = margin(t = 5, b = 8, l = 5, r = 60))
    print(myplot)
    dev.off()
  }
}
# -------------------------------------------------------------------------------------------------

# Figure S8 ---------------------------------------------------------------------------------------
itv <- 1
Titleneams <- c("(a)", "(b)", "(c)", "(d)")
Figurenames <- paste0("FigS8",c("a","b","c","d"),".png")

for(sc in 1:4)
{
  if(sc!=4)
  {
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- SDRL_result %>% filter(Interval == itv, Scenario == sc) 
    myplot <-
      tmp_result_fig %>%
      ggplot( aes(x=Delta, y=SDRL, group=Method, color=Method)) +
      geom_line(linewidth=1.5) + geom_point(size=3) + theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta), y = "SDRL",
           color = "Method") +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold")) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) #+ 
    #scale_y_continuous(limits = c(0,38))
    print(myplot)
    dev.off()
  }else{
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- SDRL_result %>% filter(Interval == itv, Scenario == sc) %>% 
      pivot_wider(names_from = Method, values_from = SDRL)
    myplot <-tmp_result_fig %>%
      ggplot( aes(x=Delta, group = 1)) +
      geom_line(aes(y=CCB), linewidth=1.5, color="#F8766D") + 
      geom_point(aes(y=CCB), size=3, color="#F8766D") + 
      geom_line(aes(y=CCK/13), linewidth=1.5, color="#00BFC4") + 
      geom_point(aes(y=CCK/13), size=3, color="#00BFC4") +
      theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta)) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) + 
      scale_y_continuous(name = "SDRL (CCB)", limits=c(0,5.5), sec.axis = sec_axis(~.*13, name="SDRL (CCK)") ) +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold"),
            axis.text.y.left = element_text(color = "#F8766D"),
            axis.text.y.right = element_text(color = "#00BFC4"),
            axis.title.y.left = element_text(color = "#F8766D"),
            axis.title.y.right = element_text(color = "#00BFC4"),
            plot.margin = margin(t = 5, b = 8, l = 5, r = 60))
    print(myplot)
    dev.off()
  }
}
# -------------------------------------------------------------------------------------------------


# Figure S10 ---------------------------------------------------------------------------------------
itv <- 2
Titleneams <- c("(a)", "(b)", "(c)", "(d)")
Figurenames <- paste0("FigS10",c("a","b","c","d"),".png")

for(sc in 1:4)
{
  if(sc!=4)
  {
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- ARL1_result %>% filter(Interval == itv, Scenario == sc) 
    myplot <-
      tmp_result_fig %>%
      ggplot( aes(x=Delta, y=ARL1, group=Method, color=Method)) +
      geom_line(linewidth=1.5) + geom_point(size=3) + theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta), y = "ARL1",
           color = "Method") +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold")) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) #+ 
    #scale_y_continuous(limits = c(0,38))
    print(myplot)
    dev.off()
  }else{
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- ARL1_result %>% filter(Interval == itv, Scenario == sc) %>% 
      pivot_wider(names_from = Method, values_from = ARL1)
    myplot <-tmp_result_fig %>%
      ggplot( aes(x=Delta, group = 1)) +
      geom_line(aes(y=CCB), linewidth=1.5, color="#F8766D") + 
      geom_point(aes(y=CCB), size=3, color="#F8766D") + 
      geom_line(aes(y=CCK/2.1), linewidth=1.5, color="#00BFC4") + 
      geom_point(aes(y=CCK/2.1), size=3, color="#00BFC4") +
      theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta)) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) + 
      scale_y_continuous(name = "ARL1 (CCB)", limits=c(0,50), sec.axis = sec_axis(~.*2.1, name="ARL1 (CCK)") ) +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold"),
            axis.text.y.left = element_text(color = "#F8766D"),
            axis.text.y.right = element_text(color = "#00BFC4"),
            axis.title.y.left = element_text(color = "#F8766D"),
            axis.title.y.right = element_text(color = "#00BFC4"),
            plot.margin = margin(t = 5, b = 8, l = 5, r = 60))
    print(myplot)
    dev.off()
  }
}
# -------------------------------------------------------------------------------------------------

# Figure S12 ---------------------------------------------------------------------------------------
itv <- 2
Titleneams <- c("(a)", "(b)", "(c)", "(d)")
Figurenames <- paste0("FigS12",c("a","b","c","d"),".png")

for(sc in 1:4)
{
  if(sc!=4)
  {
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- SDRL_result %>% filter(Interval == itv, Scenario == sc) 
    myplot <-
      tmp_result_fig %>%
      ggplot( aes(x=Delta, y=SDRL, group=Method, color=Method)) +
      geom_line(linewidth=1.5) + geom_point(size=3) + theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta), y = "SDRL",
           color = "Method") +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold")) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) #+ 
    #scale_y_continuous(limits = c(0,38))
    print(myplot)
    dev.off()
  }else{
    png(Figurenames[sc], width = 200, height = 150, units='mm', res = 900)
    tmp_result_fig <- SDRL_result %>% filter(Interval == itv, Scenario == sc) %>% 
      pivot_wider(names_from = Method, values_from = SDRL)
    myplot <-tmp_result_fig %>%
      ggplot( aes(x=Delta, group = 1)) +
      geom_line(aes(y=CCB), linewidth=1.5, color="#F8766D") + 
      geom_point(aes(y=CCB), size=3, color="#F8766D") + 
      geom_line(aes(y=CCK/7.3), linewidth=1.5, color="#00BFC4") + 
      geom_point(aes(y=CCK/7.3), size=3, color="#00BFC4") +
      theme_bw() +
      labs(title=Titleneams[sc], x =expression(delta)) +
      scale_x_discrete(breaks=c("d1","d2","d3","d4"), labels=c("0.05", "0.10", "0.15", "0.20")) + 
      scale_y_continuous(name = "SDRL (CCB)", limits=c(0,10), sec.axis = sec_axis(~.*7.3, name="SDRL (CCK)") ) +
      theme(plot.title  = element_text(size=40, face = "bold", hjust = 0.5),
            legend.title= element_text(size=27),
            legend.text = element_text(size=22),
            axis.text   = element_text(size=22),
            axis.title  = element_text(size=27,face="bold"),
            axis.text.y.left = element_text(color = "#F8766D"),
            axis.text.y.right = element_text(color = "#00BFC4"),
            axis.title.y.left = element_text(color = "#F8766D"),
            axis.title.y.right = element_text(color = "#00BFC4"),
            plot.margin = margin(t = 5, b = 8, l = 5, r = 60))
    print(myplot)
    dev.off()
  }
}
# -------------------------------------------------------------------------------------------------