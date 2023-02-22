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
    # Proportion of correct items in the set
    tmp <- match(rankedItems, true_rank[tc,])
    
    results[c, "No. Correct Items"] <- length(tmp[!is.na(tmp)])
    results[c, "No. Missclassification"] <- n_star - results[c,"No. Correct Items"]
    results[c, "Prop. Correct Items"] <- results[c,"No. Correct Items"]/n_star 
  }
  return(list(matrix = results, rankedItems = rankedItems))
}







