#### If there are 0 assessors 
#### Two three measures of the results 
#### Trace plots 
#### Tempering 
#### Acceptance rate of M-H, around 30%, count 
#### functions in other files 
library(ggplot2)

Avg_ranks <- function(rho_matrix, A_matrix){
  avg <- c()
  for(i in 1:n){
    idx <- which(A_matrix==i)
    avg <- c(avg, sum(rho_matrix[idx]/length(idx)))
    
  }
  return(sort(avg, index.return=T, na.last = T)$ix)
}



result_matrix <- function(burnin_percentage, cluster_mcmc, C, true_index_cluster,rho_mcmc,A_mcmc, n_star,true_rank,thinning=10){
  ##### Initalization for result matrix #####
  M <- floor(M/thinning)
  burnin = M*burnin_percentage 
  res_names <- c("Cluster", "True Cluster","No. Assessors","No. Correct Assessors","No. Correct Items", "No. Missclassification", "Prop. Correct Items")
  results <- matrix(NA, nrow = C, ncol = length(res_names))
  colnames(results) <- res_names
  avg_rho_ranks <- rankedItems <- list()
  
  ##### Creating matrix for the results #####
  clus <- c()
  for(k in 1:ncol(cluster_mcmc)){
    clus <- c(clus,as.integer(names(which.max(table(cluster_mcmc[burnin:M,k])))))
  }
  
  
  for(c in 1:C){
    index_cluster <- which(clus==c)
    results[c,"No. Assessors"] <- length(index_cluster)
    results[c, "Cluster"] <- c 
    correct_assesssors <- apply(true_index_cluster, 1, function(row) sum(!is.na(match(row, index_cluster))))
    results[c, "No. Correct Assessors"] <- max(correct_assesssors)
    tc <- which.max(correct_assesssors)
    results[c, "True Cluster"] <- tc
    
    rho_use <- rho_mcmc[[c]][burnin:M,]
    A_use <- A_mcmc[[c]][burnin:M,]
    
    rankedItems <- Avg_ranks(rho_use, A_use)
    # Proportion of correct items in the set
    tmp <- match(rankedItems[1:n_star], true_rank)
    results[c, "No. Correct Items"] <- length(tmp[!is.na(tmp)])
    results[c, "No. Missclassification"] <- n_star - results[c,"No. Correct Items"]
    results[c, "Prop. Correct Items"] <- results[c,"No. Correct Items"]/n_star 
  }
  return(list(matrix = results, rankedItems = rankedItems))
}




trace_plot <- function(clusters, number = 5, spread = 100){
  assessors <- sample.int(n = ncol(clusters), size = number)
  matrix_data <- clusters[seq(1, nrow(clusters), spread), assessors ]
  x <- rep(1:nrow(matrix_data))
  name <- rep(assessors, each = nrow(matrix_data))
  df <- data.frame(cbind(x, matrix_data))
  y_cols <- names(df) %>%
    setdiff(c("x"))
  
  # Create the trace plot
  
  y_col <- y_cols[1]
  fig <- plot_ly(df, x = ~x) 
  for (i in seq_along(y_cols)) {
    y_col <- y_cols[i]
    fig <- fig %>% add_trace(y = as.formula(paste0("~`", y_col, "`")), name = y_col,type = 'scatter', mode = 'lines') 
  }
  return(fig)
}



boxplot_cluster <- function(clusters){
  C <-  max(clusters, na.rm = TRUE)
  box <- vector(mode = "list", length = C)
  max_integers <- apply(clusters, 2, function(x) {
    max(table(x))/nrow(clusters)
  })
  for(i in 1:N){
    cluster_nr <- which.max(table(clusters[,i]))
    box[[cluster_nr]] <- c(box[[cluster_nr]], max_integers[i])
  }
  boxplot(box[[1]])
  print(box)
  my_matrix <- do.call(rbind, box)
  return(my_matrix)
}