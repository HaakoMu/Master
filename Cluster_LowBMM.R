library(BayesMallows)
library(Rcpp)
library(dplyr)
library(gtools)
library(plotly)
sourceCpp(file.path(".","/Cpp/leapandshift_sourceCpp.cpp")) # use the C++ leap-and-shift function
source("MCMC.R")
source("Simulations.R")
source("Results.R")
library(ggplot2)
source("Functions.R")
source("Visual.R")


##### Fixing parametes #####
N <- 80 # number of assessors 
n <- 40 # total number of items
n_star <- 12 # number of items selected to be "relevant", i.e. will have the highest ranks
n_star_true <- n_star
alpha_true <- alpha0 <- 5
C <- 4#number of clusters
thinning = 5

# Tuning parmameters for the MCMC
psi <- 10
M <- 4e4
leap_size = round(n_star/5)
L <- 2
prob_back <- prob_forw <- 0.5
burnin <- 0.3
A_star <- matrix(c(1,20,2,19,3,18,4,17,5,16,6,15,7,14,8,13),nrow = 2)

#simulation <- sim_rank_consistency(n,N,C,n_star)

simulation <- sim_topK(n,N,C,n_star)

init <- generate_random_init(M,N,C,n_star)

#cluster_mcmc <- function(data, M, N,n_star, alpha0, rho0, A_star0, clusters )

mcmc <- cluster_mcmc(simulation$data, M, C, N, n_star, alpha0, init$rho0, init$A_star0, init$clusters, prob_back, prob_forw)


#function(burnin, cluster_mcmc, C, true_index_cluster,rho_mcmc,A_mcmc, n_star,true_rank){
results <- result_matrix(burnin = burnin, cluster_mcmc = mcmc$clusters, C = C, true_index_cluster =simulation$true_index_cluster,rho_mcmc= mcmc$rho_mcmc,A_mcmc = mcmc$A_mcmc, n_star = n_star,true_rank = simulation$true_ranks)


trace <- trace_plot(mcmc$clusters)

cluster_boxplot <- boxplot_cluster(mcmc$clusters, burnin)




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
############### VISUAL ###############
#### HEATMAP for items  ####

heat <- heatplot_rho(A = mcmc$A_mcmc, rho = mcmc$rho_mcmc, burnin, n, n_star, C, simulation$true_ranks, results$matrix)

heat_1 <- heatplot_rho(mcmc$A_mcmc, mcmc$rho_mcmc, burnin, n, n_star, C)


#### Barplot for MAP of A* items and rho ####

bar <- barplot_item(mcmc$A_mcmc, mcmc$rho_mcmc, burnin, n = n, n_star = n_star)




