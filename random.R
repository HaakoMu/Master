source(MCMC)

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
