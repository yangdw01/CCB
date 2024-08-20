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



###################################################################################################
# Bayesian approach based on SUGS algorithm
# -------------------------------------------------------------------------------------------------
calculate_parameters_NCT <- function(tmp_a,tmp_b,tmp_m,tmp_psi){
  
  hat_psi <- 1/(1/tmp_psi + 1)
  
  nu <- 2*tmp_a
  muy <- hat_psi * (1/tmp_psi) * tmp_m / (1-hat_psi) 
  sig2 <- (1/nu) * ((2*tmp_b + (tmp_m^2)*(1/tmp_psi - hat_psi/(tmp_psi^2)))/(1-hat_psi) - muy^2)
  sig <- sqrt(sig2)
  
  return(list(nu, muy, sig))
}

Noncentral_t <- function(y, others){
  
  nu <- others[[1]]
  mu <- others[[2]]
  sig <- others[[3]]
  
  tmp_val1 <- gamma(0.5)/beta(nu/2, 0.5) / (((pi*nu)^(1/2)) * sig)
  tmp_val2 <- (1 + 1/((sig^2)*nu) * ((y-mu)^2))^(-(nu+1)/2)
  result <- tmp_val1 * tmp_val2
  return(result)
}

SUGS_algorithm_sequential <- function(Y, a=1, m=0, psi=1, c=1, d=10, alpha_stars, eta_probs=NULL, first_iter=TRUE,
                                      prev_bhat_j=NULL, prev_Gamma_i=NULL, prev_phi_t=NULL, prev_psi_hj=NULL, prev_m_hj=NULL, prev_a_hj=NULL, prev_b_hj=NULL)
{
  # hyperparameters
  if(is.null(eta_probs)){ eta_probs <- rep(1/length(alpha_stars), length(alpha_stars))}
  
  # basic setting
  n <- nrow(Y)
  K <- ncol(Y)
  TT <- length(alpha_stars)
  
  # parameters
  bhat_j <- Gamma_i <- phi_t <- NULL
  psi_hj <- m_hj <- a_hj <- b_hj <- list()
  
  ## 0th value
  zero_psi_hj <- rep(psi, K)
  zero_m_hj <- rep(m, K)
  zero_a_hj <- rep(a, K)
  
  # algorithm
  
  if(first_iter == TRUE){
    
    previous_n <- 0
    
    ## 1st iteration
    i <- 1
    
    ### update b
    old_bhat_j <- rep(c/d, K)
    new_bhat_j <- rep(c/d, K)
    bhat_j <- rbind(bhat_j, new_bhat_j)
    
    ### choose gamma
    new_gamma_i <- 1
    Gamma_i <- c(Gamma_i, new_gamma_i)
    
    ### update alpha
    old_phi_t <- eta_probs
    new_phi_t <- eta_probs
    phi_t <- rbind(phi_t, new_phi_t)
    
    ### distribution of theta
    old_psi_hj <- zero_psi_hj
    new_psi_hj <- 1/(1/old_psi_hj + 1)
    
    old_m_hj <- zero_m_hj
    new_m_hj <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
    
    old_a_hj <- zero_a_hj
    new_a_hj <- old_a_hj + 1/2
    
    old_b_hj <- old_bhat_j
    new_b_hj <- old_b_hj + (Y[i,]^2 + (old_m_hj^2)/old_psi_hj - (new_m_hj^2)/new_psi_hj)/2 - old_bhat_j + new_bhat_j
    new_b0_hj <- old_b_hj - old_bhat_j + new_bhat_j
    
    psi_hj[[new_gamma_i]] <- t(new_psi_hj)
    m_hj[[new_gamma_i]] <- t(new_m_hj)
    a_hj[[new_gamma_i]] <- t(new_a_hj)
    b_hj[[new_gamma_i]] <- t(new_b_hj)
    
    ## next iterations
    for(i in 2:n)
    {
      cat(i, "   ")
      
      k_i1 <- length(unique(Gamma_i))
      
      ### update b
      old_bhat_j <- bhat_j[i-1,]
      new_bhat_j <- (c + a*k_i1) / (d + ((sapply(1:k_i1, function(x){ a_hj[[x]][i-1,] }) / sapply(1:k_i1, function(x){ b_hj[[x]][i-1,] })) %>% apply(1,sum)))
      bhat_j <- rbind(bhat_j, new_bhat_j)
      
      ### choose gamma
      tmp_m_hj <- matrix(zero_m_hj, nrow=i-1, ncol=K, byrow=T)
      tmp_psi_hj <- matrix(zero_psi_hj, nrow=i-1, ncol=K, byrow=T)
      tmp_a_hj <- matrix(zero_a_hj, nrow=i-1, ncol=K, byrow=T)
      tmp_b_hj <- matrix(old_bhat_j, nrow=i-1, ncol=K, byrow=T)
      
      old_phi_t <- phi_t[i-1,]
      pi_iht <- rbind(matrix(table(Gamma_i), nrow=k_i1, ncol=TT, byrow=F), alpha_stars) * matrix(1/((i-1) + alpha_stars), nrow=k_i1+1, ncol=TT, byrow=T)
      
      tmp_probs <- pi_iht * matrix(old_phi_t, nrow=k_i1+1, ncol=TT, byrow=T)
      Lih <- sapply(1:k_i1, function(x){ Noncentral_t(y=Y[i,], others=calculate_parameters_NCT(tmp_a=a_hj[[x]][i-1,], tmp_b=b_hj[[x]][i-1,], tmp_m=m_hj[[x]][i-1,], tmp_psi=psi_hj[[x]][i-1,])) }) %>% apply(2, prod)
      Lih <- c(Lih, Noncentral_t(y=Y[i,], others=calculate_parameters_NCT(tmp_a=tmp_a_hj[i-1,], tmp_b=tmp_b_hj[i-1,], tmp_m=tmp_m_hj[i-1,], tmp_psi=tmp_psi_hj[i-1,])) %>% prod)
      
      prob_gammai <- (tmp_probs*matrix(Lih, nrow=k_i1+1, ncol=TT)) %>% apply(1, sum)
      new_gamma_i <- prob_gammai %>% which.max()
      Gamma_i <- c(Gamma_i, new_gamma_i)
      
      ### update alpha
      old_phi_t <- phi_t[i-1,]
      new_phi_t <- old_phi_t * pi_iht[new_gamma_i,]
      new_phi_t <- new_phi_t/sum(new_phi_t)
      phi_t <- rbind(phi_t, new_phi_t)
      
      ### update theta
      if(sum(Gamma_i == new_gamma_i)==1)
      {
        for(ii in 1:k_i1){
          
          psi_hj[[ii]] <- rbind(psi_hj[[ii]], psi_hj[[ii]][i-1,])
          m_hj[[ii]]   <- rbind(m_hj[[ii]],   m_hj[[ii]][i-1,])
          a_hj[[ii]]   <- rbind(a_hj[[ii]],   a_hj[[ii]][i-1,])
          b_hj[[ii]]   <- rbind(b_hj[[ii]],   b_hj[[ii]][i-1,] - old_bhat_j + new_bhat_j)
        }
        
        old_psi_hj <- zero_psi_hj
        new_psi_hj <- 1/(1/old_psi_hj + 1)
        
        old_m_hj <- zero_m_hj
        new_m_hj <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
        
        old_a_hj <- zero_a_hj
        new_a_hj <- old_a_hj + 1/2
        
        old_b_hj <- old_bhat_j
        new_b_hj <- old_b_hj + (Y[i,]^2 + (old_m_hj^2)/old_psi_hj - (new_m_hj^2)/new_psi_hj)/2 - old_bhat_j + new_bhat_j
        
        psi_hj[[new_gamma_i]] <- rbind(matrix(zero_psi_hj, nrow=i-1, ncol=K, byrow=T), new_psi_hj)
        m_hj[[new_gamma_i]]   <- rbind(matrix(zero_m_hj, nrow=i-1, ncol=K, byrow=T), new_m_hj)
        a_hj[[new_gamma_i]]   <- rbind(matrix(zero_a_hj, nrow=i-1, ncol=K, byrow=T), new_a_hj)
        b_hj[[new_gamma_i]]   <- rbind(bhat_j[1:(i-1),], new_b_hj)
        
      }else{
        
        for(ii in 1:k_i1){
          
          if(ii == new_gamma_i){
            
            old_psi_hj <- psi_hj[[ii]][i-1,]
            new_psi_hj <- 1/(1/old_psi_hj + 1)
            
            old_m_hj <- m_hj[[ii]][i-1,]
            new_m_hj <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
            
            old_a_hj <- a_hj[[ii]][i-1,]
            new_a_hj <- old_a_hj + 1/2
            
            old_b_hj <- b_hj[[ii]][i-1,]
            new_b_hj <- old_b_hj + (Y[i,]^2 + (old_m_hj^2)/old_psi_hj - (new_m_hj^2)/new_psi_hj)/2 - old_bhat_j + new_bhat_j
            
            psi_hj[[new_gamma_i]] <- rbind(psi_hj[[new_gamma_i]], new_psi_hj)
            m_hj[[new_gamma_i]]   <- rbind(m_hj[[new_gamma_i]],   new_m_hj)
            a_hj[[new_gamma_i]]   <- rbind(a_hj[[new_gamma_i]],   new_a_hj)
            b_hj[[new_gamma_i]]   <- rbind(b_hj[[new_gamma_i]],   new_b_hj)
            
          }else{
            
            psi_hj[[ii]] <- rbind(psi_hj[[ii]], psi_hj[[ii]][i-1,])
            m_hj[[ii]]   <- rbind(m_hj[[ii]],   m_hj[[ii]][i-1,])
            a_hj[[ii]]   <- rbind(a_hj[[ii]],   a_hj[[ii]][i-1,])
            b_hj[[ii]]   <- rbind(b_hj[[ii]],   b_hj[[ii]][i-1,] - old_bhat_j + new_bhat_j)
          }
        }
      }
    }
    
  }else{
    
    # basic previous setting
    previous_n <- length(prev_Gamma_i)
    
    bhat_j  <- prev_bhat_j
    Gamma_i <- prev_Gamma_i
    phi_t   <- prev_phi_t
    psi_hj  <- prev_psi_hj
    m_hj    <- prev_m_hj
    a_hj    <- prev_a_hj
    b_hj    <- prev_b_hj
    
    for(i in 1:n)
    {
      cat(previous_n+i, "   ")
      
      k_i1 <- length(unique(Gamma_i))
      
      ### update b
      old_bhat_j <- bhat_j[previous_n+i-1,]
      new_bhat_j <- (c + a*k_i1) / (d + ((sapply(1:k_i1, function(x){ a_hj[[x]][previous_n+i-1,] }) / sapply(1:k_i1, function(x){ b_hj[[x]][previous_n+i-1,] })) %>% apply(1,sum)))
      bhat_j <- rbind(bhat_j, new_bhat_j)
      
      ### choose gamma
      tmp_m_hj   <- matrix(zero_m_hj,   nrow=previous_n+i-1, ncol=K, byrow=T)
      tmp_psi_hj <- matrix(zero_psi_hj, nrow=previous_n+i-1, ncol=K, byrow=T)
      tmp_a_hj   <- matrix(zero_a_hj,   nrow=previous_n+i-1, ncol=K, byrow=T)
      tmp_b_hj   <- matrix(old_bhat_j,  nrow=previous_n+i-1, ncol=K, byrow=T)
      
      old_phi_t <- phi_t[previous_n+i-1,]
      pi_iht <- rbind(matrix(table(Gamma_i), nrow=k_i1, ncol=TT, byrow=F), alpha_stars) * matrix(1/((previous_n+i-1) + alpha_stars), nrow=k_i1+1, ncol=TT, byrow=T)
      
      tmp_probs <- pi_iht * matrix(old_phi_t, nrow=k_i1+1, ncol=TT, byrow=T)
      Lih <- sapply(1:k_i1, function(x){ Noncentral_t(y=Y[i,], others=calculate_parameters_NCT(tmp_a=a_hj[[x]][previous_n+i-1,], tmp_b=b_hj[[x]][previous_n+i-1,], tmp_m=m_hj[[x]][previous_n+i-1,], tmp_psi=psi_hj[[x]][previous_n+i-1,])) }) %>% apply(2, prod)
      Lih <- c(Lih, Noncentral_t(y=Y[i,], others=calculate_parameters_NCT(tmp_a=tmp_a_hj[previous_n+i-1,], tmp_b=tmp_b_hj[previous_n+i-1,], tmp_m=tmp_m_hj[previous_n+i-1,], tmp_psi=tmp_psi_hj[previous_n+i-1,])) %>% prod)
      
      prob_gammai <- (tmp_probs*matrix(Lih, nrow=k_i1+1, ncol=TT)) %>% apply(1, sum)
      new_gamma_i <- prob_gammai %>% which.max()
      Gamma_i <- c(Gamma_i, new_gamma_i)
      
      ### update alpha
      old_phi_t <- phi_t[previous_n+i-1,]
      new_phi_t <- old_phi_t * pi_iht[new_gamma_i,]
      new_phi_t <- new_phi_t/sum(new_phi_t)
      phi_t <- rbind(phi_t, new_phi_t)
      
      ### update theta
      if(sum(Gamma_i == new_gamma_i)==1)
      {
        for(ii in 1:k_i1){
          
          psi_hj[[ii]] <- rbind(psi_hj[[ii]], psi_hj[[ii]][previous_n+i-1,])
          m_hj[[ii]]   <- rbind(m_hj[[ii]],   m_hj[[ii]][previous_n+i-1,])
          a_hj[[ii]]   <- rbind(a_hj[[ii]],   a_hj[[ii]][previous_n+i-1,])
          b_hj[[ii]]   <- rbind(b_hj[[ii]],   b_hj[[ii]][previous_n+i-1,] - old_bhat_j + new_bhat_j)
        }
        
        old_psi_hj <- zero_psi_hj
        new_psi_hj <- 1/(1/old_psi_hj + 1)
        
        old_m_hj   <- zero_m_hj
        new_m_hj   <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
        
        old_a_hj   <- zero_a_hj
        new_a_hj   <- old_a_hj + 1/2
        
        old_b_hj   <- old_bhat_j
        new_b_hj   <- old_b_hj + (Y[i,]^2 + (old_m_hj^2)/old_psi_hj - (new_m_hj^2)/new_psi_hj)/2 - old_bhat_j + new_bhat_j
        
        psi_hj[[new_gamma_i]] <- rbind(matrix(zero_psi_hj, nrow=previous_n+i-1, ncol=K, byrow=T), new_psi_hj)
        m_hj[[new_gamma_i]]   <- rbind(matrix(zero_m_hj,   nrow=previous_n+i-1, ncol=K, byrow=T), new_m_hj)
        a_hj[[new_gamma_i]]   <- rbind(matrix(zero_a_hj,   nrow=previous_n+i-1, ncol=K, byrow=T), new_a_hj)
        b_hj[[new_gamma_i]]   <- rbind(bhat_j[1:(previous_n+i-1),],                               new_b_hj)
        
      }else{
        
        for(ii in 1:k_i1){
          
          if(ii == new_gamma_i){
            
            old_psi_hj <- psi_hj[[ii]][previous_n+i-1,]
            new_psi_hj <- 1/(1/old_psi_hj + 1)
            
            old_m_hj <- m_hj[[ii]][previous_n+i-1,]
            new_m_hj <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
            
            old_a_hj <- a_hj[[ii]][previous_n+i-1,]
            new_a_hj <- old_a_hj + 1/2
            
            old_b_hj <- b_hj[[ii]][previous_n+i-1,]
            new_b_hj <- old_b_hj + (Y[i,]^2 + (old_m_hj^2)/old_psi_hj - (new_m_hj^2)/new_psi_hj)/2 - old_bhat_j + new_bhat_j
            
            psi_hj[[new_gamma_i]] <- rbind(psi_hj[[new_gamma_i]], new_psi_hj)
            m_hj[[new_gamma_i]]   <- rbind(m_hj[[new_gamma_i]],   new_m_hj)
            a_hj[[new_gamma_i]]   <- rbind(a_hj[[new_gamma_i]],   new_a_hj)
            b_hj[[new_gamma_i]]   <- rbind(b_hj[[new_gamma_i]],   new_b_hj)
            
          }else{
            
            psi_hj[[ii]] <- rbind(psi_hj[[ii]], psi_hj[[ii]][previous_n+i-1,])
            m_hj[[ii]]   <- rbind(m_hj[[ii]],   m_hj[[ii]][previous_n+i-1,])
            a_hj[[ii]]   <- rbind(a_hj[[ii]],   a_hj[[ii]][previous_n+i-1,])
            b_hj[[ii]]   <- rbind(b_hj[[ii]],   b_hj[[ii]][previous_n+i-1,] - old_bhat_j + new_bhat_j)
          }
        }
      }
    }
  }
  
  # approximate PML
  kn <- Gamma_i %>% unique %>% length
  
  final_bhat_j <- bhat_j[previous_n+n,]
  final_phi_t <- phi_t[previous_n+n,]
  final_pi <- rbind(matrix(table(Gamma_i), nrow=kn, ncol=TT, byrow=F), alpha_stars) * matrix(1/(previous_n+n + alpha_stars), nrow=kn+1, ncol=TT, byrow=T)
  final_prob <- (final_pi * matrix(final_phi_t, nrow=kn+1, ncol=TT, byrow=T)) %>% apply(1,sum)
  final_Lih <- NULL
  for(x in 1:kn){
    
    final_Lih <- rbind(final_Lih, sapply(1:n, function(ii){ Noncentral_t(y=Y[ii,], others=calculate_parameters_NCT(tmp_a=a_hj[[x]][previous_n+n,], tmp_b=b_hj[[x]][previous_n+n,], tmp_m=m_hj[[x]][previous_n+n,], tmp_psi=psi_hj[[x]][previous_n+n,])) }) %>% apply(2, prod) )
  }
  final_Lih <- rbind(final_Lih, sapply(1:n, function(ii){ Noncentral_t(y=Y[ii,], others=calculate_parameters_NCT(tmp_a=zero_a_hj, tmp_b=final_bhat_j, tmp_m=zero_m_hj, tmp_psi=zero_psi_hj)) }) %>% apply(2, prod) )
  log_aPML <- (final_Lih * matrix( final_prob, nrow=kn+1, ncol=n )) %>% apply(2, sum) %>% log %>% sum()
  
  return(list(log_aPML, bhat_j, Gamma_i, phi_t, psi_hj, m_hj, a_hj, b_hj))
}

SUGS_algorithm_sequential_fast <- function(Y, a=1, m=0, psi=1, c=1, d=10, alpha_stars, eta_probs=NULL, first_iter=TRUE,
                                           prev_bhat_j=NULL, prev_Gamma_i=NULL, prev_phi_t=NULL, prev_psi_hj=NULL, prev_m_hj=NULL, prev_a_hj=NULL, prev_b_hj=NULL)
{
  # hyperparameters
  if(is.null(eta_probs)){ eta_probs <- rep(1/length(alpha_stars), length(alpha_stars))}
  
  # basic setting
  n <- nrow(Y)
  K <- ncol(Y)
  TT <- length(alpha_stars)
  
  # parameters
  bhat_j <- Gamma_i <- phi_t <- NULL
  psi_hj <- m_hj <- a_hj <- b_hj <- list()
  
  ## 0th value
  zero_psi_hj <- rep(psi, K)
  zero_m_hj <- rep(m, K)
  zero_a_hj <- rep(a, K)
  
  # algorithm
  
  if(first_iter == TRUE){
    
    previous_n <- 0
    previous_n0 <- 0
    
    ## 1st iteration
    i <- 1
    
    ### update b
    old_bhat_j <- rep(c/d, K)
    new_bhat_j <- rep(c/d, K)
    bhat_j <- rbind(bhat_j, new_bhat_j)
    
    ### choose gamma
    new_gamma_i <- 1
    Gamma_i <- c(Gamma_i, new_gamma_i)
    
    ### update alpha
    old_phi_t <- eta_probs
    new_phi_t <- eta_probs
    phi_t <- rbind(phi_t, new_phi_t)
    
    ### distribution of theta
    old_psi_hj <- zero_psi_hj
    new_psi_hj <- 1/(1/old_psi_hj + 1)
    
    old_m_hj <- zero_m_hj
    new_m_hj <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
    
    old_a_hj <- zero_a_hj
    new_a_hj <- old_a_hj + 1/2
    
    old_b_hj <- old_bhat_j
    new_b_hj <- old_b_hj + (Y[i,]^2 + (old_m_hj^2)/old_psi_hj - (new_m_hj^2)/new_psi_hj)/2 - old_bhat_j + new_bhat_j
    new_b0_hj <- old_b_hj - old_bhat_j + new_bhat_j
    
    psi_hj[[new_gamma_i]] <- t(new_psi_hj)
    m_hj[[new_gamma_i]] <- t(new_m_hj)
    a_hj[[new_gamma_i]] <- t(new_a_hj)
    b_hj[[new_gamma_i]] <- t(new_b_hj)
    
    ## next iterations
    for(i in 2:n)
    {
      k_i1 <- length(unique(Gamma_i))
      
      ### update b
      old_bhat_j <- bhat_j[i-1,]
      new_bhat_j <- (c + a*k_i1) / (d + ((sapply(1:k_i1, function(x){ a_hj[[x]][i-1,] }) / sapply(1:k_i1, function(x){ b_hj[[x]][i-1,] })) %>% apply(1,sum)))
      bhat_j <- rbind(bhat_j, new_bhat_j)
      
      ### choose gamma
      tmp_m_hj <- matrix(zero_m_hj, nrow=i-1, ncol=K, byrow=T)
      tmp_psi_hj <- matrix(zero_psi_hj, nrow=i-1, ncol=K, byrow=T)
      tmp_a_hj <- matrix(zero_a_hj, nrow=i-1, ncol=K, byrow=T)
      tmp_b_hj <- matrix(old_bhat_j, nrow=i-1, ncol=K, byrow=T)
      
      old_phi_t <- phi_t[i-1,]
      pi_iht <- rbind(matrix(table(Gamma_i), nrow=k_i1, ncol=TT, byrow=F), alpha_stars) * matrix(1/((i-1) + alpha_stars), nrow=k_i1+1, ncol=TT, byrow=T)
      
      tmp_probs <- pi_iht * matrix(old_phi_t, nrow=k_i1+1, ncol=TT, byrow=T)
      Lih <- sapply(1:k_i1, function(x){ Noncentral_t(y=Y[i,], others=calculate_parameters_NCT(tmp_a=a_hj[[x]][i-1,], tmp_b=b_hj[[x]][i-1,], tmp_m=m_hj[[x]][i-1,], tmp_psi=psi_hj[[x]][i-1,])) }) %>% apply(2, prod)
      Lih <- c(Lih, Noncentral_t(y=Y[i,], others=calculate_parameters_NCT(tmp_a=tmp_a_hj[i-1,], tmp_b=tmp_b_hj[i-1,], tmp_m=tmp_m_hj[i-1,], tmp_psi=tmp_psi_hj[i-1,])) %>% prod)
      
      prob_gammai <- (tmp_probs*matrix(Lih, nrow=k_i1+1, ncol=TT)) %>% apply(1, sum)
      new_gamma_i <- prob_gammai %>% which.max()
      Gamma_i <- c(Gamma_i, new_gamma_i)
      
      ### update alpha
      old_phi_t <- phi_t[i-1,]
      new_phi_t <- old_phi_t * pi_iht[new_gamma_i,]
      new_phi_t <- new_phi_t/sum(new_phi_t)
      phi_t <- rbind(phi_t, new_phi_t)
      
      ### update theta
      if(sum(Gamma_i == new_gamma_i)==1)
      {
        for(ii in 1:k_i1){
          
          psi_hj[[ii]] <- rbind(psi_hj[[ii]], psi_hj[[ii]][i-1,])
          m_hj[[ii]]   <- rbind(m_hj[[ii]],   m_hj[[ii]][i-1,])
          a_hj[[ii]]   <- rbind(a_hj[[ii]],   a_hj[[ii]][i-1,])
          b_hj[[ii]]   <- rbind(b_hj[[ii]],   b_hj[[ii]][i-1,] - old_bhat_j + new_bhat_j)
        }
        
        old_psi_hj <- zero_psi_hj
        new_psi_hj <- 1/(1/old_psi_hj + 1)
        
        old_m_hj <- zero_m_hj
        new_m_hj <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
        
        old_a_hj <- zero_a_hj
        new_a_hj <- old_a_hj + 1/2
        
        old_b_hj <- old_bhat_j
        new_b_hj <- old_b_hj + (Y[i,]^2 + (old_m_hj^2)/old_psi_hj - (new_m_hj^2)/new_psi_hj)/2 - old_bhat_j + new_bhat_j
        
        psi_hj[[new_gamma_i]] <- rbind(matrix(zero_psi_hj, nrow=i-1, ncol=K, byrow=T), new_psi_hj)
        m_hj[[new_gamma_i]]   <- rbind(matrix(zero_m_hj, nrow=i-1, ncol=K, byrow=T), new_m_hj)
        a_hj[[new_gamma_i]]   <- rbind(matrix(zero_a_hj, nrow=i-1, ncol=K, byrow=T), new_a_hj)
        b_hj[[new_gamma_i]]   <- rbind(bhat_j[1:(i-1),], new_b_hj)
        
      }else{
        
        for(ii in 1:k_i1){
          
          if(ii == new_gamma_i){
            
            old_psi_hj <- psi_hj[[ii]][i-1,]
            new_psi_hj <- 1/(1/old_psi_hj + 1)
            
            old_m_hj <- m_hj[[ii]][i-1,]
            new_m_hj <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
            
            old_a_hj <- a_hj[[ii]][i-1,]
            new_a_hj <- old_a_hj + 1/2
            
            old_b_hj <- b_hj[[ii]][i-1,]
            new_b_hj <- old_b_hj + (Y[i,]^2 + (old_m_hj^2)/old_psi_hj - (new_m_hj^2)/new_psi_hj)/2 - old_bhat_j + new_bhat_j
            
            psi_hj[[new_gamma_i]] <- rbind(psi_hj[[new_gamma_i]], new_psi_hj)
            m_hj[[new_gamma_i]]   <- rbind(m_hj[[new_gamma_i]],   new_m_hj)
            a_hj[[new_gamma_i]]   <- rbind(a_hj[[new_gamma_i]],   new_a_hj)
            b_hj[[new_gamma_i]]   <- rbind(b_hj[[new_gamma_i]],   new_b_hj)
            
          }else{
            
            psi_hj[[ii]] <- rbind(psi_hj[[ii]], psi_hj[[ii]][i-1,])
            m_hj[[ii]]   <- rbind(m_hj[[ii]],   m_hj[[ii]][i-1,])
            a_hj[[ii]]   <- rbind(a_hj[[ii]],   a_hj[[ii]][i-1,])
            b_hj[[ii]]   <- rbind(b_hj[[ii]],   b_hj[[ii]][i-1,] - old_bhat_j + new_bhat_j)
          }
        }
      }
    }
    
  }else{
    
    # basic previous setting
    previous_n <- length(prev_Gamma_i)
    previous_n0 <- nrow(prev_phi_t)
    
    bhat_j  <- prev_bhat_j
    Gamma_i <- prev_Gamma_i
    phi_t   <- prev_phi_t
    psi_hj  <- prev_psi_hj
    m_hj    <- prev_m_hj
    a_hj    <- prev_a_hj
    b_hj    <- prev_b_hj
    
    for(i in 1:n)
    {
      k_i1 <- length(unique(Gamma_i))
      
      ### update b
      old_bhat_j <- bhat_j[previous_n0+i-1,]
      new_bhat_j <- (c + a*k_i1) / (d + ((sapply(1:k_i1, function(x){ a_hj[[x]][previous_n0+i-1,] }) / sapply(1:k_i1, function(x){ b_hj[[x]][previous_n0+i-1,] })) %>% apply(1,sum)))
      bhat_j <- rbind(bhat_j, new_bhat_j)
      
      ### choose gamma
      tmp_m_hj   <- matrix(zero_m_hj,   nrow=previous_n0+i-1, ncol=K, byrow=T)
      tmp_psi_hj <- matrix(zero_psi_hj, nrow=previous_n0+i-1, ncol=K, byrow=T)
      tmp_a_hj   <- matrix(zero_a_hj,   nrow=previous_n0+i-1, ncol=K, byrow=T)
      tmp_b_hj   <- matrix(old_bhat_j,  nrow=previous_n0+i-1, ncol=K, byrow=T)
      
      old_phi_t <- phi_t[previous_n0+i-1,]
      pi_iht <- rbind(matrix(table(Gamma_i), nrow=k_i1, ncol=TT, byrow=F), alpha_stars) * matrix(1/((previous_n+i-1) + alpha_stars), nrow=k_i1+1, ncol=TT, byrow=T)
      
      tmp_probs <- pi_iht * matrix(old_phi_t, nrow=k_i1+1, ncol=TT, byrow=T)
      Lih <- sapply(1:k_i1, function(x){ Noncentral_t(y=Y[i,], others=calculate_parameters_NCT(tmp_a=a_hj[[x]][previous_n0+i-1,], tmp_b=b_hj[[x]][previous_n0+i-1,], tmp_m=m_hj[[x]][previous_n0+i-1,], tmp_psi=psi_hj[[x]][previous_n0+i-1,])) }) %>% apply(2, prod)
      Lih <- c(Lih, Noncentral_t(y=Y[i,], others=calculate_parameters_NCT(tmp_a=tmp_a_hj[previous_n0+i-1,], tmp_b=tmp_b_hj[previous_n0+i-1,], tmp_m=tmp_m_hj[previous_n0+i-1,], tmp_psi=tmp_psi_hj[previous_n0+i-1,])) %>% prod)
      
      prob_gammai <- (tmp_probs*matrix(Lih, nrow=k_i1+1, ncol=TT)) %>% apply(1, sum)
      new_gamma_i <- prob_gammai %>% which.max()
      Gamma_i <- c(Gamma_i, new_gamma_i)
      
      ### update alpha
      old_phi_t <- phi_t[previous_n0+i-1,]
      new_phi_t <- old_phi_t * pi_iht[new_gamma_i,]
      new_phi_t <- new_phi_t/sum(new_phi_t)
      phi_t <- rbind(phi_t, new_phi_t)
      
      ### update theta
      if(sum(Gamma_i == new_gamma_i)==1)
      {
        for(ii in 1:k_i1){
          
          psi_hj[[ii]] <- rbind(psi_hj[[ii]], psi_hj[[ii]][previous_n0+i-1,])
          m_hj[[ii]]   <- rbind(m_hj[[ii]],   m_hj[[ii]][previous_n0+i-1,])
          a_hj[[ii]]   <- rbind(a_hj[[ii]],   a_hj[[ii]][previous_n0+i-1,])
          b_hj[[ii]]   <- rbind(b_hj[[ii]],   b_hj[[ii]][previous_n0+i-1,] - old_bhat_j + new_bhat_j)
        }
        
        old_psi_hj <- zero_psi_hj
        new_psi_hj <- 1/(1/old_psi_hj + 1)
        
        old_m_hj   <- zero_m_hj
        new_m_hj   <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
        
        old_a_hj   <- zero_a_hj
        new_a_hj   <- old_a_hj + 1/2
        
        old_b_hj   <- old_bhat_j
        new_b_hj   <- old_b_hj + (Y[i,]^2 + (old_m_hj^2)/old_psi_hj - (new_m_hj^2)/new_psi_hj)/2 - old_bhat_j + new_bhat_j
        
        psi_hj[[new_gamma_i]] <- rbind(matrix(zero_psi_hj, nrow=previous_n0+i-1, ncol=K, byrow=T), new_psi_hj)
        m_hj[[new_gamma_i]]   <- rbind(matrix(zero_m_hj,   nrow=previous_n0+i-1, ncol=K, byrow=T), new_m_hj)
        a_hj[[new_gamma_i]]   <- rbind(matrix(zero_a_hj,   nrow=previous_n0+i-1, ncol=K, byrow=T), new_a_hj)
        b_hj[[new_gamma_i]]   <- rbind(bhat_j[1:(previous_n0+i-1),],                               new_b_hj)
        
      }else{
        
        for(ii in 1:k_i1){
          
          if(ii == new_gamma_i){
            
            old_psi_hj <- psi_hj[[ii]][previous_n0+i-1,]
            new_psi_hj <- 1/(1/old_psi_hj + 1)
            
            old_m_hj <- m_hj[[ii]][previous_n0+i-1,]
            new_m_hj <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
            
            old_a_hj <- a_hj[[ii]][previous_n0+i-1,]
            new_a_hj <- old_a_hj + 1/2
            
            old_b_hj <- b_hj[[ii]][previous_n0+i-1,]
            new_b_hj <- old_b_hj + (Y[i,]^2 + (old_m_hj^2)/old_psi_hj - (new_m_hj^2)/new_psi_hj)/2 - old_bhat_j + new_bhat_j
            
            psi_hj[[new_gamma_i]] <- rbind(psi_hj[[new_gamma_i]], new_psi_hj)
            m_hj[[new_gamma_i]]   <- rbind(m_hj[[new_gamma_i]],   new_m_hj)
            a_hj[[new_gamma_i]]   <- rbind(a_hj[[new_gamma_i]],   new_a_hj)
            b_hj[[new_gamma_i]]   <- rbind(b_hj[[new_gamma_i]],   new_b_hj)
            
          }else{
            
            psi_hj[[ii]] <- rbind(psi_hj[[ii]], psi_hj[[ii]][previous_n0+i-1,])
            m_hj[[ii]]   <- rbind(m_hj[[ii]],   m_hj[[ii]][previous_n0+i-1,])
            a_hj[[ii]]   <- rbind(a_hj[[ii]],   a_hj[[ii]][previous_n0+i-1,])
            b_hj[[ii]]   <- rbind(b_hj[[ii]],   b_hj[[ii]][previous_n0+i-1,] - old_bhat_j + new_bhat_j)
          }
        }
      }
    }
  }
  
  # approximate PML
  kn <- Gamma_i %>% unique %>% length
  
  final_bhat_j <- bhat_j[previous_n0+n,]
  final_phi_t <- phi_t[previous_n0+n,]
  final_pi <- rbind(matrix(table(Gamma_i), nrow=kn, ncol=TT, byrow=F), alpha_stars) * matrix(1/(previous_n+n + alpha_stars), nrow=kn+1, ncol=TT, byrow=T)
  final_prob <- (final_pi * matrix(final_phi_t, nrow=kn+1, ncol=TT, byrow=T)) %>% apply(1,sum)
  final_Lih <- NULL
  for(x in 1:kn){
    
    final_Lih <- rbind(final_Lih, sapply(1:n, function(ii){ Noncentral_t(y=Y[ii,], others=calculate_parameters_NCT(tmp_a=a_hj[[x]][previous_n0+n,], tmp_b=b_hj[[x]][previous_n0+n,], tmp_m=m_hj[[x]][previous_n0+n,], tmp_psi=psi_hj[[x]][previous_n0+n,])) }) %>% apply(2, prod) )
  }
  final_Lih <- rbind(final_Lih, sapply(1:n, function(ii){ Noncentral_t(y=Y[ii,], others=calculate_parameters_NCT(tmp_a=zero_a_hj, tmp_b=final_bhat_j, tmp_m=zero_m_hj, tmp_psi=zero_psi_hj)) }) %>% apply(2, prod) )
  log_aPML <- (final_Lih * matrix( final_prob, nrow=kn+1, ncol=n )) %>% apply(2, sum) %>% log %>% sum()
  
  return(list(log_aPML, bhat_j, Gamma_i, phi_t, psi_hj, m_hj, a_hj, b_hj))
}
# -------------------------------------------------------------------------------------------------
#
###################################################################################################