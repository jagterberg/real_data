library(igraph)
library(nonparGraphTesting)


dat <- read.csv("stormofswords.csv")
g <- graph.data.frame(dat)
A <- as_adjacency_matrix(g)



A_eigen <- eigen(A)
ds <- A_eigen$values
