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
N <- 200 # number of assessors 
n <- 40 # total number of items
n_star <- 12 # number of items selected to be "relevant", i.e. will have the highest ranks
n_star_true <- n_star
alpha_true <- alpha0 <- 10
C <- 4#number of clusters
thinning = 10

# Tuning parmameters for the MCMC
psi <- 50
M <- 2e4
leap_size = round(n_star/5)
L <- 2
leap_size = 2
prob_back <- prob_forw <- 0.5
burnin <- 0.8
A_star <- matrix(c(1,20,2,19,3,18,4,17,5,16,6,15,7,14,8,13),nrow = 2)
A_star_1 <- matrix(c(1,1,40,40,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12), nrow=2)  
A <- NULL
A <- rbind(A, c(1:12))
A <- rbind(A, c(1:12))
A <- rbind(A, c(9:20))
A <- rbind(A, c(16:27))

#simulation <- sim_rank_consistency(n,N,C,n_star, true_A_star = A)

simulation <- sim_topK(n,N,C,n_star, true_A_star = A)

init <- generate_random_init(M,N,C,n_star)

#cluster_mcmc <- function(data, M, N,n_star, alpha0, rho0, A_star0, clusters )

mcmc <- cluster_mcmc(simulation$data, M, C, N, n_star, alpha0, init$rho0, init$A_star0, init$clusters, prob_back, prob_forw, leap_size = leap_size, L=L)


#function(burnin, cluster_mcmc, C, true_index_cluster,rho_mcmc,A_mcmc, n_star,true_rank){
results <- result_matrix(burnin = burnin, cluster_mcmc = mcmc$clusters, C = C, true_index_cluster =simulation$true_index_cluster,rho_mcmc= mcmc$rho_mcmc,A_mcmc = mcmc$A_mcmc, n_star = n_star,true_rank = simulation$true_ranks)


trace <- trace_plot(mcmc$clusters)

cluster_boxplot <- boxplot_cluster(mcmc$clusters, burnin)


barplot(out[, rev(order(out[1,]))[1:n_star]], main = "Probability of items from MAP", xlab = "Column Index", ylab = "Value")
############### VISUAL ###############
#### HEATMAP for items  ####

heat <- heatplot_rho(A = mcmc$A_mcmc, rho = mcmc$rho_mcmc, burnin, n, n_star, C, simulation$true_ranks, results$matrix)

heat_1 <- heatplot_rho(mcmc$A_mcmc, mcmc$rho_mcmc, burnin, n, n_star, C)

#### Barplot for MAP of A* items and rho ####

bar <- barplot_item(mcmc$A_mcmc, mcmc$rho_mcmc, burnin, n = n, n_star = n_star)



s <- MAP(cluster =mcmc$clusters, A_star = mcmc$A_mcmc, rho = mcmc$rho_mcmc, n =n, burnin = burnin)


t <- new_map(cluster =mcmc$clusters, A = mcmc$A_mcmc, rho = mcmc$rho_mcmc,  burnin = burnin, n=n)



dis_check_all <- NULL
for(j in 1:4){
  dis_check <- NULL
  for(i in 101:150){
    ndata <- rank(simulation$data[i,t[j,]], ties.method = "min")
    dis_check <- rbind(dis_check, sum(abs(ndata- c(1:n_star))))
  }
  dis_check_all <- cbind(dis_check_all, dis_check)
}


table(apply(tail(mcmc$clusters), 2, get_cluster))
ga <- NULL
for(i in 1:200){
  dat <- table(tail(mcmc$clusters[,i]))
  ga <- c(ga,as.integer(names(dat)[which.max(dat)]))
}
