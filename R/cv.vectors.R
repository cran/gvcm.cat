cv.vectors <-
function(X, y, method, family, lambda, phi, L.index, T.index, 
weights, control, dev, d, l, oml, indices, ...)             
{
n <- dim(X)[1]

if(control$tuning.criterion == "deviance"){crit <- dev}
if(control$tuning.criterion == "SSE"){crit <- function(y,mudach,weights){sum((y-mudach)^2*weights)}}
if(control$cv.refit == FALSE){  ###
   evalcv <- function(X.tr, X.te, beta, control, training.y, training.weights) {
        output <- list(beta=beta, test.X=X.te)
        return(output)
        }   
   } else {
   evalcv <- function(X.tr, X.te, beta, control, training.y, training.weights){
        coefficients <- round(beta,digits=control$accuracy)
        reductionC <- reduce(coefficients, indices, control$assured.intercept)$C
        X.reduced.tr <- as.matrix(X.tr %*% reductionC)     
        X.reduced.te <- as.matrix(X.te %*% reductionC) 
        try(beta.refitted <- as.matrix(round(glm.fit(X.reduced.tr,training.y, training.weights, family=family,
            intercept = FALSE)$coefficients, digits=control$accuracy)))
        output <- list(beta=beta.refitted, test.X=X.reduced.te)
        return(output)
       }
   }

losses <- matrix(ncol=length(lambda), nrow=length(phi))
rownames(losses) <- as.character(phi)
colnames(losses) <- as.character(lambda)

for (r in 1: length(phi)) {
     weight <- weight.function(phi[r], indices, oml, control)

for (i in 1:length(lambda)) {
    loss <- 0
    for (j in 1:control$K){
        training.X <- X[L.index[[j]],]  
        training.y <- y[L.index[[j]]]
        training.weights <- weights[L.index[[j]]]
        test.X <- X[T.index[[j]],]
        test.y <- y[T.index[[j]]]
        test.weights <- weights[T.index[[j]]]
        
        model <- optimierung(X=training.X, y=training.y, method, family, lambda=lambda[i], weight, training.weights, control, l, oml, indices, phi, FALSE)
        eval.cv <- evalcv(training.X, test.X, model$beta.i, control, training.y, training.weights)#
        test.mudach <- family$linkinv(eval.cv$test.X %*% eval.cv$beta)                            #
#        test.mudach <- family$linkinv(test.X %*% model$beta.i)
        loss <- loss + crit(y=test.y, mudach=test.mudach, weights=test.weights)
    }
    losses[r,i] <- loss 
}

}

opt <- max(lambda[(which(losses==min(losses))-1)%/%dim(losses)[1]+1])
if (length(opt)!=1) {opt <- opt[1]}
phi <- phi[(which(losses==min(losses))-1)%%(dim(losses)[1])+1]
if(length(phi)!=1) {phi <-phi[which.min((phi-0.5)^2)]}
if(length(phi)!=1) {phi <- phi[1]}

return(list(lambda=opt, phi=phi, score=losses, lambdas=lambda))
}

