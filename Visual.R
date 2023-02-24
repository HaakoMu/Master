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

heatplot_rho <- function(A, rho, burnin, n, n_star,C, ranks = NULL, result = NULL){
  plots <- list()
  for(c in 1:C){
    A_mat <- burnin_mat(A[[c]], burnin)
    Rho_mat <- burnin_mat(rho[[c]], burnin)
    df <- NULL
    if(is.null(ranks)) items <- c(1:n) else items <- ranks[result[c,2],]
    for(i in items){
      data <- table(Rho_mat[which(A_mat == i, arr.ind = TRUE)])
      z <- rep(0,n_star)
      z[as.integer(names(data))] <- as.vector(data)/sum(as.vector(data))
      df <- rbind(df, data.frame(
        x = rep(paste(i), n_star),
        y = c(1:n_star),
        z = z,
        value = rep(sum(as.vector(data)/n_star, n_star))))
    }
    df <- df[order(df$z ==0, df$y, -ifelse(df$z == 0, NA, df$z)),]
    g1 <- ggplot(df, aes(x=fct_inorder(x), y=value)) +
      geom_bar(color = "black", stat = "identity") +
      ylab("Proportion") +
      theme_minimal() +
      theme(axis.text.x=element_blank(), 
            axis.title.x = element_blank(), 
            axis.ticks.x=element_blank(), 
            axis.title = element_text(size = 14),
            plot.margin=unit(c(1,1,-1,1), "cm")) 
    
    
    
    g3 <- ggplot(df, aes(fct_inorder(x), y, fill= z) , de) + 
      geom_tile(color = "grey")+
      #theme_ipsum() +
      xlab("Item") +
      ylab("Ranking") +
      theme_minimal() +
      scale_y_continuous(breaks= pretty_breaks()) +
      #scale_x_discrete(breaks= pretty_breaks()) +
      theme(axis.text=element_text(size=11),
            axis.title=element_text(size=14),
            text = element_text(family = "sans"), 
            axis.ticks.x=element_blank(),
            #axis.text.x=element_blank(), 
            legend.title = element_text(size = 11),
            legend.position = c(1.14, 0.3),
            legend.text = element_text(size = 11)) +
      #scale_fill_gradient(low="blue", high="red",name="Posterior \n probability")
      
      scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5, name="Posterior \n probability")
    
    gg <- ggarrange(g1, g3, nrow = 2, ncol = 1, widths = 1, heights = c(0.4, 3))
    plots[[c]] <- gg
  }
  return(plots)
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



