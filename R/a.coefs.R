a.coefs <-
function(index1, index2, index3, assured.intercept=TRUE, p.ord.abs=TRUE) {

  a.coefs.nominal <- function (p = NULL, ...)   
      { if (p < 2)
          stop ("There must be at least two levels, if the variable varies! \n")

        a.coefs.mat <- diag(p)
        for (i in 1:(p-1)){
        a.coefs.mat <- cbind(a.coefs.mat, 
                rbind(matrix(0,ncol=p-i,nrow=i-1), -1, diag(1,p-i))) 
        }

        return (a.coefs.mat) 
      }

  a.coefs.ordinal <- function (p = NULL, ...)  
      { if (p < 2)
          stop ("There must be at least two levels, if the variable varies! \n")

        if (p > 2)
          { h1 <- cbind (-diag (p-1), 0)
            h2 <- cbind (0, diag (p-1))
            mat1 <- h1 + h2
            mat2 <- diag (p)
            a.coefs.mat <- cbind (mat2, t (mat1)) }
        else { a.coefs.mat <- cbind (diag (2), c(-1,1)) }

        return (a.coefs.mat) 
      }

  A <- matrix(0,ncol=0,nrow=0)
  if (assured.intercept==TRUE){
  if(index2[1]<0){A <- bdiag(A,a.coefs.nominal(p=index1[1])[,-(1:(index1[1]))])}  
  if(index2[1]>0){A <- bdiag(A,a.coefs.ordinal(p=index1[1])[,-(1:(index1[1]))])}
  if(index2[1]==0){A <- bdiag(A,matrix(0,ncol=0,nrow=index1[1]))}
  b <- 2
  } else {b <- 1}

  if (b <= length(index1)) {
  for (i in b:length(index1)) {
  if(index2[i]<0){A <- bdiag(A,a.coefs.nominal(p=index1[i]))}
  if(index2[i]>0){A <- bdiag(A,a.coefs.ordinal(p=index1[i]))}
  
  if(index3[i]<0){A <- bdiag(A,a.coefs.nominal(p=index1[i]))}
  if(index3[i]>0){
     if(p.ord.abs == TRUE){
     A <- bdiag(A,a.coefs.ordinal(p=index1[i]))
     } else {
     A <- bdiag(A,a.coefs.ordinal(p=index1[i])[,-(2:(index1[i]))])}
  } 
  if(index2[i]==0 && index3[i]==0){A <- bdiag(A,matrix(0,ncol=0,nrow=index1[i]))}
  }
  }

  a.coefs.mat<-as.matrix(A) 
  return(a.coefs.mat)
}

