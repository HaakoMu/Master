Avg_ranks <- function(rho_matrix, A_matrix, n){
  df = NULL
  for(i in 1:n){
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
  sample(as.integer(names(tab)[tab == max(tab)]), size=1)
}
