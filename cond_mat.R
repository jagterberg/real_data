getElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE, main="", ...) {
  ## Given a decreasingly sorted vector, return the given number of elbows
  ##
  ## Args:
  ##   dat: a input vector (e.g. a vector of standard deviations), or a input feature matrix.
  ##   n: the number of returned elbows.
  ##   threshold: either FALSE or a number. If threshold is a number, then all
  ##   the elements in d that are not larger than the threshold will be ignored.
  ##   plot: logical. When T, it depicts a scree plot with highlighted elbows.
  ##
  ## Return:
  ##   q: a vector of length n.
  ##
  ## Reference:
  ##   Zhu, Mu and Ghodsi, Ali (2006), "Automatic dimensionality selection from
  ##   the scree plot via the use of profile likelihood", Computational
  ##   Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006.
  
  #  if (is.unsorted(-d))
  
  if (is.matrix(dat)) {
    d <- sort(apply(dat,2,sd), decreasing=TRUE)
  } else {
    d <- sort(dat,decreasing=TRUE)
  }
  
  if (!is.logical(threshold))
    d <- d[d > threshold]
  
  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")
  
  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  if (n > 1 && q < (p-1)) {
    q <- c(q, q + getElbows(d[(q+1):p], n-1, plot=FALSE))
  }
  
  if (plot==TRUE) {
    if (is.matrix(dat)) {
      sdv <- d # apply(dat,2,sd)
      plot(sdv,type="b",xlab="dim",ylab="stdev",main=main,...)
      points(q,sdv[q],col=2,pch=19)
    } else {
      plot(dat, type="b",main=main,...)
      points(q,dat[q],col=2,pch=19)
    }
  }
  
  return(q)
}
library(nonparGraphTesting)
library(igraph)
library(irlba)

A1 <- read_graph("cond-mat.gml",format= "gml")
A2 <- read_graph("cond-mat-2003.gml",format= "gml")

#A1 <- graph_from_adjacency_matrix(A1)
#A2 <- graph_from_adjacency_matrix(A2)

dg1 <- decompose.graph(A1)[[1]]
dg2 <- decompose.graph(A2)[[1]]
#lcc1 <- components(A1,"strong")
#lcc2 <- components(A2,"strong")
A1 <- as_adjacency_matrix(dg1)
A2 <- as_adjacency_matrix(dg2)

rm(dg1,dg2)#,blogosphere,dg,cons,A)

print("computing eigen1")
eigen_1 <- eigen(A1)
#eigen_1_2 <- partial(-A1,120)
#svd_1 <- svd(A1)
#plot(eigen_1$values)
#plot(svd_1$d)
print("computing eigen2")
eigen_2 <- eigen(A2)
#svd_2 <- svd(A2)
#plot(eigen_2$values)
#plot(svd_2$d)


f1 <- getElbows(sort(eigen_1$values[eigen_1$values > 0],decreasing=TRUE))#,threshold = sqrt(622))#eigen_1$values > 0)
f2 <- getElbows(sort(eigen_2$values[eigen_2$values > 0],decreasing=TRUE))#,threshold = sqrt(622))#eigen_1$values > 0)
p <- min(f1,f2)

f1<- getElbows(sort(abs(eigen_1$values[eigen_1$values < 0]),decreasing=TRUE))#eigen_1$values > 0)
f2 <- getElbows(sort(abs(eigen_2$values[eigen_2$values < 0]),decreasing=TRUE))#eigen_1$values > 0)
q <-  min(f1,f2)


#d <- 103

eigen_1_ngs <- eigen_1$values[eigen_1$values < 0]
eigen_2_ngs <- eigen_2$values[eigen_2$values < 0]

#get the ASEs:
Xhat <- cbind(eigen_1$vectors[,c(1:p)] %*% diag(eigen_1$values[c(1:p)]^(1/2)),
              eigen_1$vectors[,c( (length(eigen_1_ngs)+ 1-q): length(eigen_1_ngs))]
              %*% diag(abs(eigen_1_ngs[c(1:q)])^(1/2)))

Yhat <- cbind(eigen_2$vectors[,c(1:p)] %*% diag(eigen_2$values[c(1:p)]^(1/2)),
              eigen_2$vectors[,c( (length(eigen_2_ngs)+ 1-q): length(eigen_2_ngs))]
              %*% diag(abs(eigen_2_ngs[c(1:q)])^(1/2)))


print("aligning...")
#find the alignment without negative and positive eigenvalues
get_matched <- iterative_optimal_transport(Xhat,Yhat
                                           #,lambda_init = .5
                                           #,alpha = .5
                                           #,lambda_final = .27
                                           #, Q = bdiag(1,signs[[l]]),numReps = 10
)#,p=p,q=q)

#project to block diagonal
final_Q <- procrustes(X=diag(1,p+q),Y=get_matched$Q,Pi = diag(1,p+q),p=p,q=q)

Xnew <- Xhat%*% final_Q

print("performing nonpar test...")
res <- nonpar.test(Xnew,Yhat,nsims = 500)
#res$`estimated p-value`

final_results <- list(res,Xnew,Xhat,Yhat,p,q,final_Q)
names(final_results) <- c("nonpar.test","Xnew","xhat","yhat","p","q","final_Q")

save(final_results,file = "res_12-22.Rdata")
print("finished.")




