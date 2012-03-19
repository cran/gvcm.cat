lambda.max <-
function (X, y, method, family, weight, weights, control, ...)
{
# definitions
  max.accuracy <- 1
  accuracy.option <- 7
  p <- dim(X)[2] 
  path <- matrix(nrow = p, ncol=0)
  rownames(path) <- dimnames(as.matrix(control$oml))[[1]]

# search upper bound
  option.option <- control$lambda.upper/10/accuracy.option
  option <- round(cumsum(rep(c(option.option, option.option*9), rep(accuracy.option,2))), digits=2)
  lambdas <- option[1]/2
  
  opt <- optimierung(X, y, method, family, lambda=lambdas, weight, weights, control)

  beta <- round(opt$beta.i,digits=control$accuracy)
  path <- cbind(path, matrix(beta, ncol=1)) 
  number.removed.parameters <- p-length(reduceBeta(beta,control$index1, control$index2, control$index3)) 
  
  i <- 1
  if (number.removed.parameters < control$number.selectable.parameters) {
    while (number.removed.parameters < control$number.selectable.parameters && i <= length(option)) { 
      opt <- optimierung(X, y, method, family, lambda=option[i], weight, weights, control)

      beta <- round(opt$beta.i,digits=control$accuracy)
      path <- cbind(path, matrix(beta, ncol=1)) 
      lambdas <- c(lambdas, option[i])
      number.removed.parameters <- p-length(reduceBeta(beta,control$index1, control$index2, control$index3))
      i <- i + 1  
    }  
  }
  
# reduce upper bound if possible
  if (number.removed.parameters >= control$number.selectable.parameters) {
      lambda.upper <- lambdas[i]
      lambda.lower <- ifelse(i==1,0,lambdas[i-1])
      
      while ((lambda.upper-lambda.lower) > max.accuracy) {
          lambda <- round((lambda.lower+lambda.upper)/2, digits=2)
          
          opt <- optimierung(X, y, method, family, lambda, weight, weights, control)
          beta <- round(opt$beta.i,digits=control$accuracy)
          path <- cbind(path, matrix(beta, ncol=1))       
          lambdas <- c(lambdas, lambda)
          
          number.removed.parameters <- p-length(reduceBeta(beta,control$index1, control$index2, control$index3))
          if (number.removed.parameters >= control$number.selectable.parameters) {lambda.upper <- lambda}
          if (number.removed.parameters < control$number.selectable.parameters) {lambda.lower <- lambda}
          }    
      } else {lambda.upper <- control$lambda.upper}       
  colnames(path) <- as.character(lambdas)

# order path matrix
  path <- path[,order(as.numeric(colnames(path)))]

# check upper bound
  new.lambdas <- round((1-sin(seq(pi/2,0,length.out=control$steps))[-c(1,control$steps)])*lambda.upper, digits=2)
  new.path <- matrix(nrow=dim(path)[1],ncol=length(new.lambdas))  
  colnames(new.path) <- as.character(new.lambdas)
  for (i in 1:(length(new.lambdas))){
    opt <- optimierung(X, y, method, family, lambda=new.lambdas[i], weight, weights, control)
    new.path[,i] <- round(opt$beta.i, digits=control$accuracy)
    }
  
  path <- cbind(new.path,path)  
  path <- path[,order(as.numeric(colnames(path)))]
  
  i <- 2  
  while(i<dim(path)[2] && dim(path)[1]-length(reduceBeta(path[,i],control$index1, control$index2, control$index3)) < 
         control$number.selectable.parameters) {i<-i+1}
  if (as.numeric(colnames(path)[i])<lambda.upper) {lambda.upper <- as.numeric(colnames(path)[i])}
  path <- path[,1:which(as.numeric(colnames(path))==lambda.upper)]
  
  
# upper bound too large?!  
if (number.removed.parameters < control$number.selectable.parameters && lambda.upper==10000) {
    warning("The minimal lambda setting all coefficients to zero is assumed to be 10 000. \n", call. = FALSE)
    } 

# return
  return(list(lambda.upper=lambda.upper, path=path))

}

