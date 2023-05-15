#### If there are 0 assessors 
#### Two three measures of the results 
#### Trace plots 
#### Tempering 
#### Acceptance rate of M-H, around 30%, count 
#### functions in other files 
library(ggplot2)
library(prodlim)
source("Functions.R")


#### old method for selecting the items for A^*
old_map <- function(cluster, A_star, rho, n, burnin, threshold){
  C <- length(A_star)
  clust_map <- burnin_mat(cluster, burnin)
  n_star <- ncol(A_star[[1]])
  A_map <- matrix(NA, ncol = n_star, nrow = C )
  for(c in 1:C){
    A_mat <- burnin_mat(A_star[[c]], burnin)
    rho_mat <- burnin_mat(rho[[c]], burnin)
    freq_max <- nrow(A_mat)
    df_tmp <- NULL
    for(x in 1:n){
      freq <- length(which(A_mat == x))
      if(freq > freq_max*threshold){
        data <- table(rho_mat[which(A_mat==x, arr.ind = TRUE)])
        df_tmp <- rbind(df_tmp, data.frame(x =x , freq = freq, y= sum(rho_mat[which(A_mat==x, arr.ind = TRUE)])/freq))
      }
    }
    df_tmp <- df_tmp[order(df_tmp$y,decreasing = FALSE),]
    A_map[c,] <-  head(df_tmp$x, n=n_star)
  }
  return(A_map)
}




##### Function for running the data for a set of number of cluster ####
cluster_comparions <- function(data, cluster_set, iterations, n_star, alpha0, prob_back=0.5, prob_forw=0.5, burnin ,leap_size=4, L=15 ){
  N <- nrow(data)
  n <- ncol(data)
  sample_list <- list()
  t=1
  for(k in cluster_set){
    init <- generate_random_init(iterations,N,n,k,n_star)
    mcmc <- cluster_mcmc(data, iterations, k, N, n_star, alpha0, init$rho0, init$A_star0, init$clusters, prob_back, prob_forw,leap_size = leap_size, L=L)
    sample_list[[1]] <- list(A = mcmc$A_mcmc, R = mcmc$rho_mcmc, C = mcmc$clusters)
  }
  return(sample_list)
}



#### New alternative for selecting items for posterior of A^* ####
new_map <- function( A, rho, burnin, n,treshold_freq = 0, xi = 0){
  C <- length(A)
  n_star <- ncol(A[[1]])
  A_map <- NULL
  for(c in 1:C){
    df_tmp <- NULL
    A_mat <- burnin_mat(A[[c]], burnin)
    rho_mat <- burnin_mat(rho[[c]], burnin)
    freq_max <- nrow(A_mat)
    for(x in 1:n){
      freq <- length(which(A_mat == x))
      if(freq > treshold_freq*freq_max) {
        max_value <- 0
        data <- table(rho_mat[which(A_mat==x, arr.ind = TRUE)])
        for(xx in 1:(n_star-xi)){
          vec_tmp <- c(xx:(xx+xi))
          value <- sum(data[as.character(vec_tmp)],na.rm = TRUE)/freq
          if(value>max_value) max_value<-value
        }
        df_tmp <- rbind(df_tmp, data.frame(x =x , freq = freq, d= max_value, y= sum(rho_mat[which(A_mat==x, arr.ind = TRUE)])/freq))
      }
    }
    df_tmp <- head(df_tmp[order(-df_tmp$d),], n=n_star)
    A_map <- rbind(A_map, df_tmp[order(df_tmp$y),]$x)
    #print(head(df_tmp[order(-df_tmp$d),], n=100))
  }
  return(A_map)
  
}

#### Function for the matrix of the results #####
result_df <- function(map, true_ranks, cluster, true_clus, burnin){
  df <- data.frame(matrix(NA, nrow= nrow(map), ncol = 8))
  colnames(df) <- c("True Cluster", "Simulation Cluster","Cluster size","Nr. same assessors" , "Proportion of attached assessors", "Nr. Correct Items in MAP", "Proportion of Correct Items", "Relative distance")
  clus_mat <- burnin_mat(cluster, burnin)
  true_c <- NULL
  for(i in 1:max(true_clus)){
    true_c <- c(true_c, which(true_clus== i, arr.ind = TRUE)[1])
  }
  clus_pos <- cbind(true_c, apply(clus_mat, 2, get_cluster))
  clus_avg <- NULL
  clus_tmp <- NULL
  for(c in unique(clus_pos[,1])){
      df_tmp <- clus_pos[which(clus_pos[,1]==c),]
      y <- table(df_tmp[,2])
      sim_clus <- as.numeric(names(y[which.max(y)]))
      x <- intersect(true_ranks[c,], map[sim_clus,])
      df[c, "True Cluster"] <- c
      df[c, "Simulation Cluster"] <- sim_clus
      df[c, "Nr. same assessors"] <- max(y)
      df[c, "Cluster size"] <- nrow(clus_pos[which(clus_pos[,2]==sim_clus),])
      df[c, "Proportion of attached assessors"] <- max(y)/nrow(df_tmp)
      df[c, "Nr. Correct Items in MAP"] <- length(x)
      df[c, "Proportion of Correct Items"] <- length(x)/ncol(true_ranks)
      df[c, "Relative distance"] <- sum(abs(rank(match(x, true_ranks[c,]))- rank(match(x, map[sim_clus,]))))
  }
  return(df)
}




#### Function for running multiple chains on the same data ####
multichain <- function(data, A, R, chains = 5, n, N , n_star, M, burnin , L ,leap_size, prob_back= 0.5, prob_forw = 0.5,  psi = 50, C, thinning = 10, alpha0 = 3){
  output <-  list()
  for(k in 1:chains){
    init <- generate_random_init(M,N,n,C,n_star)
    mcmc <- cluster_mcmc(data, M, C, N, n_star, alpha0, init$rho0, init$A_star0, init$clusters, prob_back, prob_forw, leap_size = leap_size, L=L)
    output[[k]] <- list(A = mcmc$A_mcmc, R = mcmc$rho_mcmc, C = mcmc$clusters)
  }
  return(output)
}


#### Function to continue the run off a MCMC with starting point being the end of the last run ####
continue_run <-function(data, A, R, Clus, n , N , n_star, M  , burnin , L , leap_size , prob_back= 0.5, prob_forw = 0.5,  psi = 50,  thinning = 10, alpha0 = 3){
  C <- length(A)
  A_0 <- R_0 <- C_0 <- list()
  for(i in 1:C){
    A_0[[i]] <- tail(A[[i]],1)
    R_0[[i]] <- tail(R[[i]],1)
  }
  clusters <- matrix(NA, floor(M/thinning), N)
  clusters[1,] <- tail(Clus,1)
  mcmc <- cluster_mcmc(data, M, C, N, n_star, alpha0, R_0, A_0, clusters, prob_back, prob_forw, leap_size = leap_size, L=L)
  return(list(A = mcmc$A_mcmc, R = mcmc$rho_mcmc, C=mcmc$clusters))
}





