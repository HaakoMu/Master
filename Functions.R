Avg_ranks <- function(rho_matrix, A_matrix, n, set = NULL){
  df = NULL
  if(is.null(set)) set <- c(1:n)
  for(i in set){
    idx <- which(A_matrix ==  i)
    df <- rbind(df, data.frame(x = paste("item ", i), y = sum(rho_matrix[idx]/length(idx))))
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
