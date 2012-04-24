cv.lambda.phi <-
function(X, y, method, family, lambda, phi, L.index, T.index, 
weights, control, dev, d, l, oml, indices, ...)           
{

# lambda vektor
if(!is.logical(lambda)){
  cross <- cv.vectors(X, y, method, family, lambda=lambda, phi, L.index, T.index, weights, control, dev, d, l, oml, indices)
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
  basis <- exp(1/4*log(control$lambda.upper))
  opt.lambdas <- matrix(nrow=2,ncol=length(phi))
  colnames(opt.lambdas) <- as.character(phi)
  scores <- list()

  for (r in 1:length(phi)){
  score <- matrix(nrow=1, ncol=0)
  lambdas <- c()
  if (basis >1) {closer <- round(basis^(0:4), digits=2)} else {closer <- seq(0,control$lambda.upper,by=control$lambda.accuracy)}
  i <- 1

  while (max(closer)-min(closer)>=control$lambda.accuracy && i <11){

      closer <- closer[which(closer<=control$lambda.upper)]
      closer <- unique(closer)
      if (any(closer %in% lambdas)) closer <- closer[-which(closer %in% lambdas)]
      if (length(closer)>0){
          cross <- cv.vectors(X, y, method, family, lambda=closer, phi[r], L.index, T.index, weights, control, dev, d, l, oml, indices)

          score <- cbind(score, cross$score)
          lambdas <- c(lambdas, cross$lambdas)
          lambda <- max(lambdas[which(score==min(score))])[1]
      }

      i <- i+1
      j <- log(lambda)/log(basis)
      closer <- round(basis^(c((j-2^(-i)), (j+2^(-i)))),digits=2)

  }

  opt.lambdas[1,r] <- lambda
  opt.lambdas[2,r] <- round((score[which(score==min(score))])[1], digits=2)
  score <- score[,order(as.numeric(colnames(score)))]
  scores[[r]] <- score
  }
  
  which.phi <- min(which(opt.lambdas[2,]==min(opt.lambdas[2,])))[1]
  colnames(opt.lambdas) <- NULL
  lambda <- opt.lambdas[1,which.phi]
  phi <- phi[which.phi]
  score <- scores[[which.phi]]
  lambdas <- as.numeric(colnames(score))
  
}
    
return(list(lambda=lambda, phi=phi, score=score, lambdas=lambdas))
}

