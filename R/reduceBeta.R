reduceBeta <-
function(beta, index1, index2, index3) # 
{ 
C <- diag(1, length(beta))
dimnames(C) <- list(rep(1:length(index1),times=index1), rownames(beta)) 

i <- 1
while (i <= dim(C)[2]) {
  if ((index2+index3)[eval(parse(text=rownames(C)[which(C[,i]==1)]))]!=0 ){ #
  if (sum(beta*C[,i])==0) {C <- as.matrix(C[,-i])} else { 
  equal <- which((beta - rep(sum(beta*C[,i]), length(beta)))==0) 
  married <- which(rownames(C)==eval(parse(text=rownames(C)[which(C[,i]==1)])) ) 
  null <- equal[which(equal %in% married)] 
  A <- C
  A[null,i] <- 1
  if (length(null)>1) {C <- as.matrix(C[,-which(colSums(A * A[,i])==1)])}
  i <- i+1 
  }
} else {i <- i+1} #
}

return(t(C) %*% beta)
      
}

