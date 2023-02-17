library(BayesMallows)
library(Rcpp)
library(dplyr)
library(gtools)
library(plotly)
sourceCpp("C:/Users/Haako/OneDrive/Documents/UiO/Master/Cluster LOWBMM/lowBMM/Cpp/leapandshift_sourceCpp.cpp") # use the C++ leap-and-shift function
source("MCMC.R")
source("Simulations.R")
source("Results.R")
library(ggplot2)


##### Fixing parametes #####
N <- 80 # number of assessors 
n <- 40 # total number of items
n_star <- 20 # number of items selected to be "relevant", i.e. will have the highest ranks
n_star_true <- n_star
alpha_true <- alpha0 <- 2
C <- 4 #number of clusters
thinning = 10 

# Tuning parmameters for the MCMC
psi <- 10
M <- 2e4
leap_size = round(n_star/5)
L <- 2
prob_back <- prob_forw <- 0.5
burnin <- 0.3
A_star <- matrix(c(1,20,2,19,3,18,4,17,5,16,6,15,7,14,8,13),nrow = 2)


simulation <- sim_topK(n,N,C,n_star)

init <- generate_random_init(M,N,C,n_star)

#cluster_mcmc <- function(data, M, N,n_star, alpha0, rho0, A_star0, clusters )
  
mcmc <- cluster_mcmc(simulation$data, M, C, N, n_star, alpha0, init$rho0, init$A_star0, init$clusters, prob_back, prob_forw)

results <- result_matrix(burnin_percentage = burnin, cluster_mcmc = mcmc$clusters, C, true_index_cluster =simulation$true_index_cluster,rho_mcmc= mcmc$rho_mcmc,A_mcmc = mcmc$A_mcmc, n_star,true_rank = simulation$true_ranks)





