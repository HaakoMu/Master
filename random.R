source(MCMC)
library(ggplot2)
test <- matrix(c(1,1,3,3,2,2,2,3,1), ncol= 3)
test2 <- matrix(c(5,5,5,6,6,6,7,7,7), ncol =3 )
sad <- replicate(8,0)
e2 <- apply(test, 1, function(row) which(row==1))

for(i in 1:nrow(test2)){
  sad[test2[i,e2[i]]] <- sad[test2[i,e2[i]]] + 1 
}

which.max(sad)
fsa <- c(4,3,NA,NA,NA)
asa <- c(1,4,2,5,6,7,8)
asa[fsa] <- 0


which(test == 1)
test[which(test==1)]


mat <- matrix(c(1,2,3,1,2,1), nrow=3, ncol=2)
table(mat[,1])



df <- data.frame(
  x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2 ,2 ,2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
        4, 4, 4, 4, 4, 4, 4, 4),
  y = c(1,  3,  4,  5,  6,  7,  8,  9, 10, 11,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
        4,  5,  6,  7,  8,  9, 10, 11, 12, 13,  1,  2,  4,  5,  7,  8,  9, 10, 11),
  z = c(0.0005, 0.0015, 0.0140, 0.0985, 0.6785, 0.1340, 0.0350, 0.0135, 0.0075,
         0.0030, 0.0020, 0.0025, 0.0095, 0.0290, 0.1650, 0.0380, 0.0140, 0.0020,
         0.0010, 0.0025, 0.0005, 0.0045, 0.0070, 0.0300, 0.1090, 0.6825, 0.1195,
         0.0305, 0.0045, 0.0035, 0.0005, 0.0015, 0.0005, 0.0005, 0.0005, 0.0015,
         0.0005, 0.0005, 0.0010, 0.0050)
)

df <- data.frame(
  x = c(1, 2, 3, 4, 5),
  y = c(1, 2, 3, 4, 5),
  z = c(10, 20, 30, 40, 50)
)

# Create a heat map using ggplot2
ggplot(df, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")




ga <- tail(mcmc$rho_mcmc[[1]],20)
fa <- tail(mcmc$A_mcmc[[1]], 20)
which(ga == 2)
table(ga[which(fa == 21)])
