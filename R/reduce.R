reduce <-
function(beta, indices, assured.intercept)
{
C <- diag(1, length(beta))
dimnames(C) <- list(rep(1:length(indices[[1]]),times=indices[[1]]), rownames(beta))

i <- 1
while (i <= dim(C)[2]) {
  if ((indices[[2]]+indices[[3]])[eval(parse(text=rownames(C)[which(C[,i]==1)]))]!=0 ){ # exclude only penalized coefficients
      # exclude zero coefficients
      if (sum(beta*C[,i])==0 && ifelse(assured.intercept, sum(as.matrix(C[,-i])[c(1:indices[[1]][1]),])!=0, TRUE)) { # keep assured.intercept, even if zero
          namen <- colnames(C)[-i]
          C <- as.matrix(C[,-i])
          colnames(C) <- namen
      } else {
      # exclude equal coefficients
          equal <- which((beta - rep(sum(beta*C[,i]), length(beta)))==0)
          married <- which(rownames(C)==eval(parse(text=rownames(C)[which(C[,i]==1)])) )
          null <- equal[which(equal %in% married)]
          A <- C
          A[null,i] <- 1
          if (length(null)>1 && ifelse(assured.intercept, sum(as.matrix(C[,-which(colSums(A * A[,i])==1)])[c(1:indices[[1]][1]),])!=0, TRUE)){
             namen <- colnames(C)[-which(colSums(A * A[,i])==1)]
             C <- as.matrix(C[,-which(colSums(A * A[,i])==1)])
             colnames(C) <- namen
             }
          i <- i+1
          }
  } else {i <- i+1}
}

return(list(C=C, beta=t(C) %*% beta))

}
