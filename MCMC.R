library(Rcpp)
library(gtools)
sourceCpp("C:/Users/Haako/OneDrive/Documents/UiO/Master/Cluster LOWBMM/lowBMM/Cpp/leapandshift_sourceCpp.cpp") # use the C++ leap-and-shift function



cluster_mcmc <- function(data, M, C, N,n_star, alpha0, rho0, A_star0, clusters,  prob_back,prob_forw, thinning = 10, leap_size=1, L = 2){
  ##### Initialization ######
  ## Proposal/old ##
  rho_old_matrix <- matrix(data = NA, nrow= C, ncol = n_star)
  A_star_old_matrix <- matrix(data= NA, nrow = C, ncol = n_star)
  current_cluster <- clusters[1,]
  A_star_prop <-rho_mcmc <- rho_prop <- A_mcmc <- ACC_rho <- RATIO_rho <- ACC_A <- RATIO_A <- list()
  for(c in 1:C){
    A_star_prop[[c]]  <- A_star_old_matrix[c,] <- A_star0[[c]]
    rho_prop[[c]] <- rho_old_matrix[c,] <- rho0[[c]]
    alpha_prop <- alpha_old <- alpha0
    rho_mcmc[[c]] <- matrix(NA, floor(M/thinning), n_star)
    rho_mcmc[[c]][1,] <- rho0[[c]]
    A_mcmc[[c]] <- matrix(0, floor(M/thinning) , n_star)
    A_mcmc[[c]][1,] <- A_star0[[c]]
    ACC_rho[[c]] <- RATIO_rho[[c]] <- ACC_A[[c]] <- RATIO_A[[c]] <- c(0)
  }
  ##  Tau ## 
  dir <- rep(C,psi)
  tau <- rdirichlet(C,dir)
  
  
  #### MCMC ####
  for(m in 2:M) {
    
    tau <- rdirichlet(C, dir+as.numeric(table(current_cluster)))
    include <- m%%thinning
    for(c in 1:C){
      ### MH step 1: update rho restricted on A* ###
      # 1a. Sample rank proposal through leap and shift (cpp function)
      assesors <- which(current_cluster==c)
      N_a <- length(assesors) #number of assesors in each cluster
      if(N_a>0){
        A_star_old <- A_star_old_matrix[c,]
        rho_old <-  rho_old_matrix[c,]
        tmp <- leap_and_shift(rho_proposal = rho_prop[[c]], indices = c(1:n_star), prob_backward = prob_back, prob_forward = prob_forw, rho = rho_old, leap_size = leap_size, reduce_indices = F)
        rho_prop[[c]] <- tmp$rho_proposal[,]
        
        
        #Selection of data using the updating clusters, with old A*
        selection <- data[assesors, A_star_old,drop = FALSE]
        for(i in 1:N_a)selection[i, sort.int(selection[i,], index.return = TRUE)$ix] <- 1:n_star
        
        # 1b. Compute distances to current and proposed ranks
        dist_new <- abs(scale(selection, rho_prop[[c]], scale = FALSE)) # footrule distance metric
        dist_old <- abs(scale(selection, rho_old, scale = FALSE))
        rank_dist_sum <- sum(dist_new-dist_old)
        
        # 1c. Compute MH ratio and accept/reject
        prob_back <- tmp$prob_backward # probability backwards
        prob_forw <- tmp$prob_forward # probability forward
        prior_rho_old <- 1 # uniform prior
        prior_rho_prop <- 1
        C_rho <- (prob_back*prior_rho_old)/(prob_forw*prior_rho_prop) # leap and shift backwards/forwards probability factor (needs to be edited!!)
        ratio_rho <- min(1, C_rho*exp((-alpha_old/n_star)*rank_dist_sum))
        
        ### Log version ###
        #ratio_rho <- log(prob_back)-log(prob_forw)-(alpha_old/n_star)*rank_dist_sum
        
        if(runif(1)<ratio_rho){
          rho_old <- rho_prop[[c]]
          ACC_rho[[c]] <- c(ACC_rho[[c]], 1)
        }else{
          ACC_rho[[c]] <- c(ACC_rho[[c]], 0)
        }
        
        ### MH step 2: update A* ###
        # 2a. Sample A* proposed by perturbing L items
        removed_items <- sample.int(n = n_star, size = L) # index of items to be removed in A*
        new_items <- sample(setdiff(1:n, A_star_old), size = L)
        A_star_prop[[c]] <- sort(c(A_star_old[-removed_items], new_items))
        
        # 2b. Update data according to new set (ranks should go from in 1,..,n*)
        
        prop_selection <- data[assesors,A_star_prop[[c]],drop = FALSE]
        for(i in 1:N_a)prop_selection[i, sort.int(prop_selection[i,], index.return = TRUE)$ix] <- 1:n_star
        
        # 2c. Compute new corresponding rho_prop based on the items selected
        rho_prop_star <- rho_prop[[c]]
        
        # Acceptance A*
        idx_match_old <- match(A_star_prop[[c]], A_star_old)[!is.na(match(A_star_prop[[c]], A_star_old))]
        idx_match_prop <- match(A_star_old, A_star_prop[[c]])[!is.na(match(A_star_old, A_star_prop[[c]]))]
        rho_prop_star[idx_match_prop] <- rho_old[idx_match_old]
        idx <- match(new_items, A_star_prop[[c]]) # index of the new item(s) in A*
        
        # Alternative: if we want to order the new items based on their avg data rankings
        #idx2 <- A_star_prop[idx]
        #idx3 <- sort(avg_ranks[idx2],index.return = TRUE)$ix # sorting new items based on ranking
        #rho_pirop_star[idx[idx3]] <- rho_old[sort(removed_items)]
        
        # Alternative: randomly assign new items their ranking, not taking data into account
        
        
        rho_prop_star[idx] <- rho_old[removed_items]
        
        # 2d. Compute MH ratio and accept/reject
        ratio_A <- min(1, exp(-alpha_old/n_star*(sum(abs(scale(prop_selection, rho_prop_star, scale = FALSE)))-
                                                   sum(abs(scale(selection, rho_old, scale = FALSE))))))
        if(runif(1)<ratio_A){
          A_star_old <- A_star_prop[[c]]
          rho_old <- rho_prop_star
          ACC_A[[c]] <- c(ACC_A[[c]], 1)
        }else{
          ACC_A[[c]] <- c(ACC_A[[c]], 0)
        }
      }
      
      rho_old_matrix[c,] <- rho_old
      A_star_old_matrix[c,] <- A_star_old
      # Save current values
      if(!include){
        rho_mcmc[[c]][floor(m/thinning),] <- rho_old
        A_mcmc[[c]][floor(m/thinning),] <- A_star_old
        RATIO_rho[[c]] <- c(RATIO_rho, ratio_rho)
        RATIO_A[[c]] <- c(RATIO_A, ratio_A)
      }
      
    }
    #Update clusters Gibbs
    for(j in 1:N){
      p_cj <- c()
      for(c in 1:C){
        #Z(alpha) må finne
        A_star_old <- A_star_old_matrix[c,]
        rho_old <- rho_old_matrix[c,]
        relabel_data <- rank(data[j,A_star_old], ties.method = "min")
        p_cj <- c(p_cj, tau[c]* exp(-alpha_old/n_star*(sum(abs(relabel_data- rho_old)))))
      }
      z_N <- which(rmultinom(1,1,p_cj)==1)
      current_cluster[j] <- z_N
      if(!include){
        clusters[floor(m/thinning),j] <- z_N
        print(m)
      }
    }
    
  }
  return(list(rho_mcmc = rho_mcmc, clusters = clusters, A_mcmc = A_mcmc, ACC_A = ACC_A, ACC_rho = ACC_rho ))
}