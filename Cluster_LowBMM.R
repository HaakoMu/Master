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
library(gridExtra)
library(ggplot2)
library(matrixStats)


set.seed(1555)

##### Fixing parametes #####
N <- 180 # number of assessors 
n <- 300 # total number of items
n_star <- 18 # number of items selected to be "relevant", i.e. will have the highest ranks
n_star_true <- 15
alpha_true <- 3
alpha0 <- 3
C <- 6#number of clusters
thinning = 10

save = FALSE

# Tuning parmameters for the MCMC
psi <- 50
M <- 4e4
leap_size = round(n_star/5)
L <- 3
prob_back <- prob_forw <- 0.5
burnin <- 0.5
bruA_star <- matrix(c(1,20,2,19,3,18,4,17,5,16,6,15,7,14,8,13),nrow = 2)
A_star_1 <- matrix(c(1,1,40,40,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12), nrow=2)  
A <- NULL
A <- rbind(A, c(1:15))
A <- rbind(A, c(9:23))
A <- rbind(A, c(20:34))
A <- rbind(A, c(1:15))
A <- rbind(A, c(9:23))
A <- rbind(A, c(20:34))


R <- NULL
R <- rbind(R, c(1:n_star_true))
R <- rbind(R, c(1:n_star_true))
R <- rbind(R, c(1:n_star_true))
R <- rbind(R, c(1:n_star_true))
R <- rbind(R, c(1:n_star_true))
R <- rbind(R, c(1:n_star_true))

simulation <- sim_data(n,N,C,n_star_true ,top_rank = 3, consistency = 3, alpha_true = alpha_true, true_A_star = A, true_rho = R)
init <- generate_random_init(M,N,n,C,n_star, thinning = thinning)
mcmc <- cluster_mcmc(simulation$data, M, C, N, n_star, alpha0, init$rho0, init$A_star0, init$clusters, prob_back, prob_forw, leap_size = leap_size, L=L)
test <- new_map( A = mcmc$A_mcmc, rho = mcmc$rho_mcmc,  burnin = burnin, n=n, treshold_freq = 0.1)
r_df <- result_df(test, simulation$true_ranks, mcmc$clusters, simulation$true_index_cluster, burnin)
r_df

t <- old_map(cluster =mcmc$clusters, A = mcmc$A_mcmc, rho = mcmc$rho_mcmc,  burnin = burnin, n=n, threshold = 0.1)
r_s <- result_df(t, simulation$true_ranks, mcmc$clusters, simulation$true_index_cluster, burnin)
r_s

heat <- heatplot(A = mcmc$A_mcmc, rho = mcmc$rho_mcmc, burnin = burnin, n = n, n_star = n_star, C =C, items = NULL, true_rank = simulation$true_ranks, cut_off_frequency = 0.1, results_df = r_df)

samples <- continue_run(data=simulation$data, A =  samples$A, R= samples$R, Clus= samples$C, n=n , N=N , n_star=n_star, M=2e4 , burnin=burnin , L=1, leap_size=leap_size, prob_back= 0.5, prob_forw = 0.5,  psi = 50,  thinning = 10, alpha0 = 3)
test1 <- new_map( A = samples$A, rho = samples$R,  burnin = burnin, n=n, treshold_freq = 0.1, xi=2)
r_df <- result_df(test1, simulation$true_ranks, mcmc$clusters, simulation$true_index_cluster, burnin)
r_df
heat1 <- heatplot(A = samples$A, rho =samples$R, burnin = burnin, n = n, n_star = n_star, C =C, items = NULL, true_rank = simulation$true_ranks, cut_off_frequency = 0.1, results_df = r_df)

#SHOWING WHICH ITEMS CHOSEN#
true_rank_visual <- selected_items(simulation$true_ranks, n)

run_top <- list()
run_cons <- list()
for(i in 1:10){
  print(i)
  N <- 180 # number of assessors 
  n <- 300 # total number of items
  n_star <- 18 # number of items selected to be "relevant", i.e. will have the highest ranks
  n_star_true <- 15
  alpha_true <- 3
  alpha0 <- 3
  C <- 6#number of clusters
  thinning = 10
  
  save = FALSE
  
  # Tuning parmameters for the MCMC
  psi <- 50
  M <- 5e4
  leap_size = round(n_star/5)
  L <- 4
  prob_back <- prob_forw <- 0.5
  burnin <- 0.5
  bruA_star <- matrix(c(1,20,2,19,3,18,4,17,5,16,6,15,7,14,8,13),nrow = 2)
  A_star_1 <- matrix(c(1,1,40,40,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12), nrow=2)  
  A <- NULL
  A <- rbind(A, c(1:15))
  A <- rbind(A, c(9:23))
  A <- rbind(A, c(20:34))
  A <- rbind(A, c(30:44))
  A <- rbind(A, c(39:53))
  A <- rbind(A, c(49:63))
  simulation <- sim_data(n,N,C,n_star_true ,top_rank = 6, alpha_true = alpha_true, true_A_star = A)
  init <- generate_random_init(M,N,n,C,n_star, thinning = thinning)
  mcmc <- cluster_mcmc(simulation$data, M, C, N, n_star, alpha0, init$rho0, init$A_star0, init$clusters, prob_back, prob_forw, leap_size = leap_size, L=L)
  test <- new_map( A = mcmc$A_mcmc, rho = mcmc$rho_mcmc,  burnin = burnin, n=n, treshold_freq = 0.1, xi=1)
  r_df <- result_df(test, simulation$true_ranks, mcmc$clusters, simulation$true_index_cluster, burnin)
  run_cons[[i]] <- r_df
  
}


