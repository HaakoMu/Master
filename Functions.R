Avg_ranks <- function(rho_matrix, A_matrix, n, set = NULL){
  df = NULL
  if(is.null(set)) set <- c(1:n)
  for(i in set){
    idx <- which(A_matrix ==  i)
    df <- rbind(df, data.frame(x = paste("item ", i), y = sum(rho_matrix[idx]/length(idx))))
  }
  return(df)
}

get_df <- function(rho_matrix, A_matrix, item_set){
  df <- NULL
  for(i in item_set){
    idx <- which(A_matrix == i)
    new_row <- data.frame(
      name = i,
      sum = sum(rho_matrix[idx]),
      freq = length(idx)
    )
    df<- rbind(df, new_row)
  }
  return(df)
}


burnin_mat <- function(df, burnin){
  first <- burnin*nrow(df)
  last <- nrow(df)
  return(df[first:last,])
} 

get_cluster <- function(x) {
  tab <- table(x)
  as.integer(names(tab)[which.max(tab)])
}


find_match <- function(vec1, vec2){
  length(intersect(vec1,vec2))
}




save_file <- function(data, clusters_true, true_ranks, A_mcmc, rho_mcmc, cluster_mcmc, filename ){
  value_list <- list(data = data, clusters_true = clusters_true, A_mcmc = A_mcmc, rho_mcmc = rho_mcmc, cluster_mcmc = cluster_mcmc)
  save(value_list, file = paste("Data/",filename, ".Rdata", sep = ""))
}