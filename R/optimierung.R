optimierung <-
function (X, y, method, family, lambda, weight, weights, control, l, oml, indices, phi, rank = TRUE, ...)

{

if (!is.null(control$initials) && length(oml)!=length(control$initials)) { 
    warning ("control$initials do not fit, instead the default value is employed. \n ")
    }
if (!is.null(control$initials) && length(oml)==length(control$initials)) { 
    initials <- matrix(control$initials, ncol=1)
    dimnames(initials) <- dimnames(oml) 
    } else {
    initials <- oml
    }
if (any(is.na(initials))) { initials[is.na(initials)] <- 0 }

acoefs <- a.coefs(indices, control, oml)
n <- nrow(X)


if (method == "nlm") {

    # functions nlm
      link <- family$linkinv
      abs.aprox <- function(x){x2<-x^2; x2 / (sqrt(x2+control$c))}
      f <- function(beta) {- l(y=y,mudach=link(X %*% beta),weights) + lambda * sum(weight * abs.aprox(t(acoefs)%*%beta))}

    # nlm
      try(lsg <- (nlm(f, initials)) ) #suppressWarnings

    # compute return / forward to method 'lqa'
      if (exists("lsg")==TRUE) {
          beta.i <- as.matrix(lsg$estimate)
            dimnames(beta.i) <- dimnames(as.matrix(oml))
          reduction <- reduce(round(beta.i, control$accuracy), indices, control$assured.intercept)
          df.residual <- dim(X)[1] - length(reduction$beta) # df(error)
          rank <- length(reduction$beta) # df(model)
          iter <- lsg$iterations
          converged <- (lsg$code %in% c(1,2))
          } else {
          method <- "lqa"
          warning ("Method 'nlm' failed, instead method 'lqa' was used. \n")
          }
    }
    
if (method == "lqa") { # penalty lqa

      appro <- function(betak, acoefs, weight){
        diagonale <- as.vector( first.derivative(betak, lambda, weight, acoefs, indices, control)/
        ( sqrt((t(acoefs)%*%betak)^2 + control$c) ))
        output <- diag(diagonale, nrow=length(diagonale))
        return(output)
        }
}   

if (method %in% c("lqa")) {

    # functions
      A <- function(betak) {acoefs %*% appro(betak, acoefs, weight) %*% t(acoefs)}  

    # initialization
      beta.i <- initials
      eta.i <- X %*% as.vector(initials)
      beta.m <- matrix(0, nrow = control$maxi, ncol = length(initials))
      stop.at <- control$maxi
      converged <- FALSE


   # P-IRLS
     for (i in 1:control$maxi) {
        beta.m[i, ] <- beta.i
        mu.i <- family$linkinv(eta.i)
        d.i <- family$mu.eta(eta.i)
        v.i <- family$variance(mu.i)/weights
        w.wurzel.i <- as.vector(d.i/sqrt(v.i))
        X.star <- w.wurzel.i * X
        y.schlange.star <- w.wurzel.i * (eta.i + (y - mu.i)/d.i)
        A.lambda <- A(beta.i)
        p.imat.new <- crossprod(X.star) + A.lambda
        chol.pimat.new <- chol(p.imat.new)
        inv.pimat.new <- chol2inv(chol.pimat.new)
        beta.new <- control$g * drop(inv.pimat.new %*% t(X.star) %*%
            y.schlange.star) + (1 - control$g) * beta.m[i, ]
        if ((sum(abs(beta.new - beta.i))/sum(abs(beta.i)) <= control$epsilon)){
            converged <- TRUE
            stop.at <- i
            if (i < control$maxi)
                break
        } else {
            beta.i <- beta.new
            eta.i <- drop(X %*% beta.new)
        }
      }

    beta.i <- matrix(beta.i, ncol=1)
    rownames(beta.i) <- rownames(oml)
    if (rank){
      H.i.1 <- Matrix(X.star) 
      suppressWarnings(try(H.i.2 <- H.i.1 %*% inv.pimat.new))
      if (exists("H.i.2")) rank <- sum(H.i.2 * H.i.1) else rank <- 0
      } else {rank <- 0}
    df.residual <- dim(as.matrix(X))[1]- rank # df(error)
    iter <- stop.at

    if (!converged && (stop.at == control$maxi))
        cat("Convergence warning for lambda = ", lambda, "\n")

    }
    

return(list(
     beta.i = beta.i,
     df.residual = df.residual, # df(error)
     rank = rank, # df(model)
     iter = iter,
     converged = converged
     ))

}

