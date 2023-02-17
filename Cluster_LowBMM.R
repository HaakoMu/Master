library(BayesMallows)
library(Rcpp)
library(dplyr)
library(gtools)
library(plotly)
sourceCpp("C:/Users/Haako/OneDrive/Documents/UiO/Master/Cluster LOWBMM/Cluster_LowBMM/Master/Cpp/leapandshift_sourceCpp.cpp") # use the C++ leap-and-shift function
source("MCMC.R")
source("Simulations.R")
source("Results.R")
library(ggplot2)


##### Fixing parametes #####
N <- 20 # number of assessors 
n <- 40 # total number of items
n_star <- 20 # number of items selected to be "relevant", i.e. will have the highest ranks
n_star_true <- n_star
alpha_true <- alpha0 <- 2
C <- 2 #number of clusters
thinning = 10 

# Tuning parmameters for the MCMC
psi <- 10
M <- 2e4
leap_size = round(n_star/5)
L <- 2
prob_back <- prob_forw <- 0.5
burnin <- 0.3
A_star <- matrix(c(1,20,2,19,3,18,4,17,5,16,6,15,7,14,8,13),nrow = 2)

simulation <- sim_rank_consistency(n,N,C,n_star)

#simulation <- sim_topK(n,N,C,n_star)

init <- generate_random_init(M,N,C,n_star)

#cluster_mcmc <- function(data, M, N,n_star, alpha0, rho0, A_star0, clusters )
  
mcmc <- cluster_mcmc(simulation$data, M, C, N, n_star, alpha0, init$rho0, init$A_star0, init$clusters, prob_back, prob_forw)

results <- result_matrix(burnin_percentage = burnin, cluster_mcmc = mcmc$clusters, C, true_index_cluster =simulation$true_index_cluster,rho_mcmc= mcmc$rho_mcmc,A_mcmc = mcmc$A_mcmc, n_star,true_rank = simulation$true_ranks)


trace <- trace_plot(mcmc$clusters)

boxplot_mat <- boxplot_cluster(mcmc$clusters)
boxplot(t(boxplot_mat))
boxplot_mat





fix_burnin <- function(x){
  m <- nrow(x)
  b <<- burnin
  print(x)
  print(m)
  return(x[m*b:m,])
}
asd <- mcmc$A_mcmc
fff <- nrow(asd[[1]])*burnin
A_post <- lapply(asd, function(x) x[-c(1:fff), ])


my_matrix <- A_post[[1]]
paste(my_matrix[1,], collapse = ",")


count_unique_strings <- function(x) {
  paste(my_matrix[1,], collapse = ",")
}

# Apply the custom function to each row of the matrix
string_counts <- apply(my_matrix, 1, count_unique_strings)
MAP <- as.integer(strsplit((dimnames(my_table)[[1]][which.max(table(string_counts))]), ",")[[1]])


out <- matrix(nrow=1, ncol=0)

for(i in 1:n){
  col <- matrix(sum(my_matrix==i)/nrow(my_matrix), ncol = 1)
  colnames(col) <- as.character(i)
  out <- cbind(out, col,deparse.level = 0)
}

barplot(out[, rev(order(out[1,]))[1:n_star]], main = "Probability of items from MAP", xlab = "Column Index", ylab = "Value")





