path.matrix <-
function(X, y, method, family, lambda, coefficients, path, 
weight, weights, control, l, oml, indices, phi, ...)
{   
new.lambdas <- round(c(lambda-lambda/10*(7:9), (lambda+lambda/10*(7:9))[which((lambda+lambda/10*(7:9))<control$lambda.upper)]),digits=2) 

new.path <- matrix(nrow=dim(path)[1],ncol=length(new.lambdas)) 
for (i in 1:(length(new.lambdas))){
  opt <- optimierung(X, y, method, family, lambda=new.lambdas[i], weight, weights, control, l, oml, indices, phi, FALSE)
  new.path[,i] <- round(opt$beta.i, digits=control$accuracy)
  }
 
# add oml, coefficients
new.path <- cbind(oml, coefficients, new.path)
colnames(new.path) <- as.character(c(0, lambda, new.lambdas))
path <- cbind(new.path,path)  

# sortieren nach lambdas
path <- path[,order(as.numeric(colnames(path)))]

return(path) 
}

