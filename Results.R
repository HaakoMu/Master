#### If there are 0 assessors 
#### Two three measures of the results 
#### Trace plots 
#### Tempering 
#### Acceptance rate of M-H, around 30%, count 
#### functions in other files 
library(ggplot2)
source("Functions.R")


result_matrix <- function(burnin, cluster_mcmc, C, true_index_cluster,rho_mcmc,A_mcmc, n_star,true_rank){
  ##### Initalization for result matrix #####
  res_names <- c("Cluster", "True Cluster","No. Assessors","No. Correct Assessors","No. Correct Items", "No. Missclassification", "Prop. Correct Items")
  results <- matrix(NA, nrow = C, ncol = length(res_names))
  colnames(results) <- res_names
  MAP <- NULL
  
  ##### Creating matrix for the results #####
  clus_mat <- burnin_mat(cluster_mcmc,burnin)
  assessors <- ncol(cluster_mcmc)
  ass_cluster <- apply(clus_mat, 2, get_cluster)
  
  
  for(c in 1:C){
    index_cluster <- which(ass_cluster==c)
    results[c,"No. Assessors"] <- length(index_cluster)
    results[c, "Cluster"] <- c 
    correct_assesssors <- apply(true_index_cluster, 1, function(row) length(na.omit(match(index_cluster, row))))
    results[c, "No. Correct Assessors"] <- max(correct_assesssors)
    tc <- which.max(correct_assesssors)
    results[c, "True Cluster"] <- tc
    
    rho_mat <- burnin_mat(rho_mcmc[[c]],burnin)
    A_mat <- burnin_mat(A_mcmc[[c]],burnin)
    
    df <- Avg_ranks(rho_mat, A_mat,n)
    df <- df[df$y != 0,]
    df <- df[order(df$y,decreasing = FALSE),]
    rankedItems <- as.integer(rownames(head(df, n_star)))
    
    MAP <- rbind(MAP, rankedItems)
    
    # Proportion of correct items in the set
    tmp <- match(rankedItems, true_rank[tc,])
    
    results[c, "No. Correct Items"] <- length(tmp[!is.na(tmp)])
    results[c, "No. Missclassification"] <- n_star - results[c,"No. Correct Items"]
    results[c, "Prop. Correct Items"] <- results[c,"No. Correct Items"]/n_star 
  }
  return(list(matrix = results, MAP = MAP))
}



MAP <- function(cluster, A_star, rho, n, burnin){
  C <- length(A_star)
  clust_map <- burnin_mat(cluster, burnin)
  n_star <- ncol(A_star[[1]])
  A_map <- matrix(NA, ncol = n_star, nrow = C )
  for(c in 1:C){
    A_mat <- burnin_mat(A_star[[c]], burnin)
    rho_mat <- burnin_mat(rho[[c]], burnin)
    df <- Avg_ranks(rho_matrix = rho_mat,A_matrix = A_mat, n = n)
    df <- df[df$y != 0,]
    df <- df[order(df$y,decreasing = FALSE),]
    items <- as.integer(rownames(head(df, n_star)))
    A_map[c,] <- items
  }
  return(A_map)
}


cluster_comparions <- function(data, cluster_set, iterations, n_star, alpha0, prob_back=0.5, prob_forw=0.5, burnin ){
  N <- nrow(data)
  n <- ncol(data)
  result_list <- list()
  MAP_list <- list()
  cluster_list <- list()
  for(k in cluster_set){
    init <- generate_random_init(iterations,N,k,n_star)
    mcmc <- cluster_mcmc(data, iterations, k, N, n_star, alpha0, init$rho0, init$A_star0, init$clusters, prob_back, prob_forw)
    cluster_list[[k]] <- mcmc$clusters
    MAP_list[[k]] <- MAP(cluster =mcmc$clusters, A_star = mcmc$A_mcmc, rho = mcmc$rho_mcmc, n = n, burnin = burnin)
  }
  return(elbow_plot(MAP = MAP_list, data = data, cluster = cluster_list, burnin = burnin))
}




new_map <- function(cluster, A, rho, burnin, n){
  C <- length(A)
  clust_map <- burnin_mat(cluster, burnin)
  n_star <- ncol(A[[1]])
  A_map <- NULL
  for(c in 1:C){
    A_mat <- burnin_mat(A[[c]], burnin)
    rho_mat <- burnin_mat(rho[[c]], burnin)
    df <- Avg_ranks(rho_matrix = rho_mat,A_matrix = A_mat, set = unique(c(A_mat)))
    df$y <- floor(df$y)
    freq <- NULL
    newx <- NULL
    for(i in df$x){
      ii <- as.integer(strsplit(i, " ")[[1]][3])
      newx <- c(newx,ii)
      freq <- c(freq, length(which(rho_mat[which(A_mat==ii,arr.ind = TRUE),] <= df[which(df$x == i), "y"])) )
    }
    df$freq <-freq 
    df$x <- newx
    A_map <- rbind(A_map, head(df[order(df$y, -df$freq),], n= n_star)$x)
  }
  return(A_map)
}


