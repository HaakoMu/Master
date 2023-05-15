library(ggplot2)
library(reshape2)
library(hrbrthemes)
library(forcats)
library(gridExtra)
library(stringr)
library(egg)
library(dplyr)
library(scales)
library(lattice)


source("Functions.R")

heatplot <- function(A, rho, burnin, n, n_star, C , items = NULL, true_rank, cut_off_frequency = NULL, results_df){
  plot_list <- list()
  for(c in 1:C){
    tc <- results_df[which(results_df$"Simulation Cluster" == c), "True Cluster"]
    df <- NULL
    A_mat <- burnin_mat(A[[c]], burnin)
    Rho_mat <- burnin_mat(rho[[c]], burnin)
    if(is.null(items)) item_set <- c(1:n) else item_set <- items
    df <- get_df(rho_matrix = Rho_mat, A_matrix = A_mat, item_set = item_set)
    df$true_rank <- 0
    for(j in 1:ncol(true_rank)){
      item <- true_rank[tc, j]
      df$true_rank[df$name == item] <- j
    }
    if(!is.null(cut_off_frequency)) df <- df[df$freq >= max(df$freq)*cut_off_frequency,] 
    df$rank <- df$sum/df$freq
    df <- df[order(df$rank),]
    g1 <- ggplot(df, aes(x=fct_inorder(paste(name)), y=freq)) +
      geom_bar(color = "black", stat = "identity") +
      ylab("Proportion") +
      theme_minimal() +
      ggtitle(paste("Simulation cluster", c))+
      #scale_y_continuous(breaks= pretty_breaks()) +
      theme(axis.text.x=element_blank(), 
            axis.title.x = element_blank(), 
            axis.ticks.x=element_blank(), 
            axis.title = element_text(size = 14),
            plot.margin=unit(c(1,2.85,-1,0.025), "cm")) 
    
    
    col_fun <- c("#CCCCCC",rainbow(ncol(true_rank)))
    
    g2 <- ggplot(df, aes(fct_inorder(paste(name)), y = 1, fill=true_rank)) + 
      geom_tile(color = "black")+
      theme_minimal() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.x = element_blank(), 
            axis.text.y=element_blank(), 
            axis.ticks.y=element_blank(), 
            axis.title.y = element_blank(), 
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12), 
            legend.position = "right",
            plot.margin=unit(c(0.8,0.6,-0.6,1.28), "cm"))+
      #scale_colour_manual(name = "grp",values = myColors)
      scale_fill_gradientn(colours=col_fun, name="True \n ranking")
    
    
    df_heat <- NULL
    for(i in df$name){
      data <- table(Rho_mat[which(A_mat == i, arr.ind = TRUE)])
      z <- rep(0,n_star)
      z[as.integer(names(data))] <- as.vector(data)/sum(as.vector(data))
      df_heat <- rbind(df_heat, data.frame(
        x = rep(paste(i), n_star),
        y = c(1:n_star),
        z = z,
        value = rep(sum(as.vector(data), n_star))))
    }
    g3 <- ggplot(df_heat, aes(fct_inorder(x), y, fill= z) , de) + 
      geom_tile(color = "grey")+
      #theme_ipsum() +
      xlab("Item") +
      ylab("Ranking") +
      theme_minimal() +
      scale_y_continuous(breaks= pretty_breaks()) +
      #scale_x_discrete(breaks= pretty_breaks()) +
      theme(axis.text=element_text(size=11),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title=element_text(size=14),
            text = element_text(family = "sans"), 
            axis.ticks.x=element_blank(),
            #axis.text.x=element_blank(), 
            legend.title = element_text(size = 11),
            legend.position = "right",
            legend.text = element_text(size = 11)) +
      #scale_fill_gradient(low="blue", high="red",name="Posterior \n probability")
      
      scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5, name="Posterior \n probability")
    gg <- ggarrange(g1, g2, g3, nrow = 3, ncol = 1, widths = 1, heights = c(0.4, 0.2,3)) 
    plot_list[[c]] <- gg
  }
  return(plot_list)
}



heatplot_fix <- function(A, rho, burnin, n, n_star,C, ranks = NULL, result = NULL){
  plot_list <- list()
  for(i in 1:C){
    A_mat <- burnin_mat(A[[c]], burnin)
    Rho_mat <- burnin_mat(rho[[c]], burnin)
    
  }
}





barplot_item <- function(A, rho, burnin, n, n_star = NULL){
  plots <- list()
  for(c in 1:length(A)){
    A_mat <- burnin_mat(A[[c]], burnin)
    rho_mat <- burnin_mat(rho[[c]], burnin)
    df <- Avg_ranks(rho_mat, A_mat,n )
    df <-df[df$y !=0,]
    if(is.null(n_star)) nr <- nrow(df) else nr <- n_star
    df <- df[order(df$y,decreasing = FALSE),]
    p <- ggplot(data = df[1:nr,], aes(x = fct_inorder(x), y = y)) +
      geom_bar(stat = "identity", fill = "blue") +
      labs(title = "Barplot Example", x = "Item", y = "Average Rank")
    plots[[c]] <- p
  }
  return(plots)
}





boxplot_cluster <- function(clusters, burnin= NULL){
  C <-  max(clusters, na.rm = TRUE)
  if(is.null(burnin)) clusters <- clusters else clusters <- burnin_mat(clusters, burnin)
  df <- NULL
  for(i in 1:ncol(clusters)){
    data <- table(clusters[,i])
    x <- names(data)[which.max(data)]
    y <- max(data)/nrow(clusters)
    df <- rbind(df, data.frame(x = x, y = y))
  }
  ggplot(df, aes(x = x, y = y)) +
    geom_boxplot()
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


elbow_plot <- function(MAP, data, items, burnin){
  df <- NULL
  for(t in 1:length(MAP)){
    C <- nrow(MAP[[t]])
    clus_mat <- burnin_mat(items[[t]]$C, burnin)
    ass_cluster <- apply(clus_mat, 2, get_cluster)
    for(c in 1:C){
      ass_idx <- which(ass_cluster==c)
      if(length(ass_idx) !=0){
        data_tmp <- NULL
        for(i in ass_idx){
          data_tmp <- rbind(data_tmp, rank(data[i, MAP[[t]][c,]], ties.method =  "min"))
        }
        df <- rbind(df, data.frame(x= rep(paste(C),nrow(data_tmp)), y= rowSums(abs(scale(data_tmp, c(1:ncol(MAP[[1]])), scale = FALSE)))))
      }
    }
  }
  p <- ggplot(df, aes(x = x, y = y)) +
      geom_boxplot()+labs(title="Incluster Distance")+xlab("Number of clusters(C)") + ylab(expression(paste("Distance between ranking ",R[j]," and consensus ", rho[ z[j] ])))
  return(p)
}

selected_items <- function(true_rankings,n){
  C <- nrow(true_rankings)
  df <- NULL
  for(i in 1:C){
    mat_tmp <- cbind(c(1:n), rep(0,n))
    mat_tmp[true_rankings[i,],2] <- c(1:ncol(true_rankings))
    df <- rbind(df, data.frame(x = mat_tmp[,1], y = rep(i, n), z = mat_tmp[,2] ))
  }
  col_fun <- c("#CCCCCC",rainbow(ncol(true_rankings)*2))
  g3 <- ggplot(df, aes(fct_inorder(paste(x)), y, fill= z) , de) + 
    geom_tile(color = "grey")+
    #theme_ipsum() +
    xlab("Item") +
    ylab("Cluster") +
    theme_minimal() +
    scale_y_continuous(breaks= pretty_breaks()) +
    #scale_x_discrete(breaks= pretty_breaks()) +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=14),
          text = element_text(family = "sans"), 
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(), 
          legend.title = element_text(size = 11),
          legend.position = "right",
          legend.text = element_text(size = 11)) +
    #scale_fill_gradient(low="blue", high="red",name="Posterior \n probability")
    
    scale_fill_gradientn(colours=col_fun, name="True \n ranking")
}
