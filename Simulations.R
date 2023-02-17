library(BayesMallows)
sim_topK <- function(n, N, C, n_star, true_rho = NULL, true_A_star = NULL){
  true_ranks <- matrix(NA, nrow = C, ncol= n_star)
  rho_sim <- is.null(true_rho)
  if(rho_sim) true_rho <- matrix(NA, nrow = C, ncol= n_star)
  A_sim <- is.null(true_A_star)
  if(A_sim)true_A_star <- matrix(NA, nrow = C, ncol= n_star)
  true_index_cluster <-matrix(NA, nrow = C, ncol=N/C)
  data <- NULL
  for(i in 1:C){
    N_C = N/C
    if(rho_sim) rho_true <- sample.int(n = n_star, size = n_star) else rho_true <- true_rho[i,] 
    data_red <- sample_mallows(rho_true, alpha_true, N_C, burnin = 1e4, thinning = 100) # ranks/data in A*
    if(A_sim) A_star_true <- sort(sample.int(n = n, size = n_star)) else  A_star_true <- true_A_star[i,]# items in A*
    data_c <- matrix(NA, N_C, n) 
    data_c[,A_star_true] <- data_red
    true_rank <- A_star_true[sort(rho_true, index.return=TRUE)$ix]
    true_ranks[i,] <- true_rank
    true_A_star[i,] <- A_star_true
    true_rho[i,] <- rho_true
    # Complete the data by randomly assigning to the other n-n* items the ranks from n*+1 to n
    for(j in 1:(N/C))data_c[j,setdiff(1:n, A_star_true)] <- sample((n_star + 1):n, size = n-n_star)
    true_index_cluster[i,] <- c(1:N_C)+(N_C*(i-1)) 
      
    data <- rbind(data,data_c)
  }
  return(list(data = data, true_rho =true_rho, true_A_star=true_A_star, true_ranks = true_ranks, true_index_cluster= true_index_cluster ))
}


sim_rank_consistency <- function(n, N, C, n_star, true_rho = NULL, true_A_star = NULL){
  true_ranks <- matrix(NA, nrow = C, ncol= n_star)
  rho_sim <- is.null(true_rho)
  if(rho_sim) true_rho <- matrix(NA, nrow = C, ncol= n_star)
  A_sim <- is.null(true_A_star)
  if(A_sim)true_A_star <- matrix(NA, nrow = C, ncol= n_star)
  true_index_cluster <-matrix(NA, nrow = C, ncol=N/C)
  data <- NULL
  for(i in 1:C){
    N_C = N/C
    if(rho_sim) rho_true <- sample.int(n = n_star, size = n_star) else rho_true <- true_rho[i,]
    if(A_sim) A_star_true <- sort(sample.int(n = n, size = n_star)) else  A_star_true <- true_A_star[i,]
    data_red <- sample_mallows(rho_true, alpha_true, N_C, burnin = 1e4, thinning = 100)
    true_rank <- A_star_true[sort(rho_true, index.return=TRUE)$ix]
    rho_true_ordered <- sort(rho_true, index.return =TRUE)$ix
    for(j in 1:N_C)data_red[j,rho_true_ordered] <- sort(sample(1:n, size = n_star))
    data_c <- matrix(NA, N_C, n) 
    data_c[,A_star_true] <- data_red
    true_ranks[i,] <- true_rank
    true_A_star[i,] <- A_star_true
    true_rho[i,] <- rho_true
    for(j in 1:N_C)data_c[j,setdiff(1:n, A_star_true)] <- sample(setdiff(1:n, data_red[j,]), size = n-n_star)
    true_index_cluster[i,] <- c(1:N_C)+(N_C*(i-1))
    data <- rbind(data,data_c)
  }
  return(list(data = data, true_rho =true_rho, true_A_star=true_A_star, true_ranks = true_ranks, true_index_cluster= true_index_cluster ))
}

generate_random_init <- function(M,N,C,n_star, thinning = 10){
  #Clusters initializer 
  clusters <- matrix(NA, floor(M/thinning), N)
  clusters[1,] <- sample(1:C,N, replace = T)
  # Initial data in dimension (N,n*)
  init_data <- rho0 <- A_star0 <- list()
  for(j in 1:C){
    assesors <- which(clusters[1,]==j)
    N_c <- length(assesors)
    rho0[[j]] <- sample.int(n=n_star, size=n_star)
    A_star0[[j]] <- sort(sample.int(n=n, size=n_star))
  }
  return(list(clusters = clusters, rho0 = rho0, A_star0 =A_star0))
}





