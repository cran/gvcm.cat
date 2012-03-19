cv.lambda.phi <-
function(X, y, method, family, lambda, phi, L.index, T.index, 
dev, weights, control,...)
{

# lambda vektor
if(!is.logical(lambda)){
  cross <- cv.vectors(X, y, method, family, lambda=closer, phi, L.index, T.index, dev, weights, control)      
  lambda <- cross$lambda
  lambdas <- cross$lambdas
  phi <- cross$phi
  score <- cross$score
  colnames(score) <- as.character(lambdas)
}

# lambda==TRUE
if(is.logical(lambda)){

  # definitions
  p <- dim(X)[2]
  score <- matrix(nrow=length(phi), ncol=0)
  lambdas <- c()
  basis <- exp(1/4*log(control$lambda.upper))
  closer <- round(basis^(0:4), digits=2)  # lambda
  i <- 1

  while (max(closer)-min(closer)>=control$lambda.accuracy && i <11){
  
      cross <- cv.vectors(X, y, method, family, lambda=closer, phi, L.index, T.index, dev, weights, control)      

      score <- cbind(score, cross$score)
      lambdas <- c(lambdas, cross$lambdas)
      lambda <- max(lambdas[(which(score==min(score))-1)%/%dim(score)[1]+1])[1]
      
      i <- i+1
      j <- log(lambda)/log(basis)
      closer <- round(basis^(c((j-2^(-i)), (j+2^(-i)))),digits=2)
      closer <- closer[which(closer<control$lambda.upper)]
  
  }
  
  phi <- phi[(which(score==min(score))-1)%%(dim(score)[1])+1]
  if(length(phi)!=1) {phi <- phi[which.min((phi-0.5)^2)]}
  if(length(phi)!=1) {phi <- phi[1]}
  
  score <- score[,order(as.numeric(colnames(score)))]
  
}
    
return(list(lambda=lambda, phi=phi, score=score, lambdas=lambdas))
}

