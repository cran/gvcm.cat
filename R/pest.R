pest <-    
function(
X,
y,
ind,
family = gaussian,
method = c("lqa", "nlm"),
tuning = list(lambda=TRUE, phi=0.5),
weights,
control = cat_control(),
plot=FALSE,
...
)
{
# checks
if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
if (is.function(family))
    family <- family()
if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
}

method <- match.arg(method)
if (!(method %in% c("lqa", "nlm", "jump.n.nmk", "jump.n.PIRLS", "jump.root", "jump.weighted.ridge", "jump.log")))
     stop ("method is incorrect. \n")

if (!is.list(tuning) || length(tuning)!=2)
     stop ("tuning must be a list. \n")
lambda <- tuning[[1]]
phi <- tuning[[2]]
if (!is.numeric(lambda) && !is.logical(lambda))
     stop ("lambda must be numeric or 'TRUE'. \n")
if (is.numeric(lambda)){
     if (lambda<0 || !is.vector(lambda) || is.matrix(lambda) || is.list(lambda)
     || is.array(lambda))
     stop ("lambda must be numeric or 'TRUE'. \n")
     }
if (is.logical(lambda) && lambda!=TRUE)
     stop ("lambda must be numeric or 'TRUE'. \n")
if (!is.numeric(phi) && !is.logical(phi))
     stop ("phi must be numeric or 'TRUE'. \n")
if (is.numeric(phi)){
     if (phi<0 || !is.vector(phi) || is.matrix(phi) || is.list(phi)
     || is.array(phi) || phi>1)
     stop ("phi must be numeric and out of ]0;1[. \n")
     }

# model
# plot
# weights
    
# definitions
control$index1 <- ind[[1]]
control$index2 <- ind[[2]]
control$index3 <- ind[[3]] 

n <- nrow(X)
if (control$K > n)
    ("K must be a single integer < dim(data)[1]. \n")
control$number.selectable.parameters <- sum(abs(control$index1*(control$index2+control$index3)))- 
  (as.integer(control$assured.intercept))*abs(control$index2[1])

# functions
  control$psi <- 1
if(family$family == "binomial"){
  control$l <- function(y,mudach){sum((y*log(mudach) + (1-y)*log(1-mudach))*weights + lchoose(weights,y*weights))}
  control$d <- function(y,mudach){-2*log(1-abs(y-mudach))}
  }
if(family$family == "gaussian"){
  control$l <- function(y,mudach){-1/2 * sum((y - mudach)^2)/control$psi - log(sqrt(2*pi*control$psi))}
  control$d <- function(y,mudach){(y-mudach)^2}
  }
if(family$family == "poisson"){
  control$l <- function(y,mudach){sum((y*log(mudach)) - mudach - lgamma(y+1))}
  control$d <- function(y,mudach){e <- matrix(0, nrow=length(y), ncol=1)
   e[which(y==0),] <-(2*mudach)[which(y==0)]
   e[which(y!=0),] <- (-2*((y*log(mudach))-mudach-(y*log(y))+y))[which(y!=0)]
   return(e)}
  }
dev.res <- function(y,mudach,weights) {(y-mudach)/abs(y-mudach) * sqrt(control$d(y, mudach)*weights) }
dev <- function(y,mudach,weights) {sum(control$d(y=y,mudach=mudach)*weights)}

# oml
suppressWarnings(try(oml.model <- glm.fit(X, y, weights = weights,
    family = family, intercept = FALSE), silent = TRUE))
if(exists("oml.model")==FALSE) {
    oml.model <- list(coefficients=rep(NA, times=dim(X)[2]), rank=NA, aic=NA,
      iter=NA, df.residual=NA, converged=FALSE)
    warning("Ordinary maximum likelihood estimate does not exsist. \n")
    control$adapted.weights <- FALSE
    }
control$oml <- as.matrix(round(oml.model$coefficients, control$digits))
if(sum(as.integer(is.na(control$oml)))>0) {
    control$adapted.weights <- FALSE
    warning("Ordinary maximum likelihood estimate contains NAs. \n")
    }
if(length(which(abs(control$oml)<.0001))>0) {
    control$adapted.weights <- FALSE
    warning("control$adapted.weights set to FALSE as at least one ML-estimate is too close to zero. \n")
    }

# model selection
if (control$number.selectable.parameters==0){# || lambda==0){
    coefficients <- control$oml
    tuning <- list(lambda=0, phi=.5)
    rank <- oml.model$rank
    iter <- oml.model$iter
    df.residual <- oml.model$df.residual
    converged <- oml.model$converged
    plot <- list(NA,NA)
    }

if (control$number.selectable.parameters>0){# && lambda>0){

    # definitions
    if ( is.logical(phi) ) {phi <- seq(from=0.1, to=0.9, by=0.1)}
    weight <- weight.function(phi=0.5, control)
    if (is.logical(lambda) || length(lambda)>1 || length(phi)>1)
      {cross <- TRUE} else {cross <- FALSE}

    # upper boundary lambda
    if ( is.logical(lambda) || plot==TRUE ) {
    highest.lambda <- lambda.max(X, y, method, family, weight, weights, control)
    control$lambda.upper <- highest.lambda$lambda.upper
    path <- highest.lambda$path
    } else {control$lambda.upper <- NULL }

    # cross-validation
    if (cross==TRUE) {
        if (!is.list(control$K)){
        T.index <- split(sample(1:n), rep(1:control$K, length = n))
        } else {T.index <- control$K}
        L.index <- lapply(T.index, function(i) setdiff(1:n, i))

        cross <- cv.lambda.phi(X, y, method, family, lambda, phi, L.index, T.index, dev, weights, control)
        lambda <- cross$lambda
        lambdas <- cross$lambdas
        phi <- cross$phi
        score <- cross$score

        if ( is.na(lambda * phi)==TRUE ) stop("Error in cross validation. \n")
    } else {score <- NA}

    # model
    opt <- optimierung(X, y, method, family, lambda, weight, weights, control, rank=TRUE)
    coefficients <- round(opt$beta.i, digits=control$accuracy)
    tuning <- list(lambda=lambda, phi=phi)
    rank <- opt$rank
    iter <- opt$iter
    df.residual <- opt$df.residual
    converged <- opt$converged

    # plot
    if (plot==TRUE){
    path <- path.matrix(X, y, method, family, lambda, coefficients, path, weight, weights, control)
    } else {path <- NA}
    plot <- list(path=path, score=score)
    

}

# prepare output
linear.predictors <- X %*% coefficients
mudach <- family$linkinv(linear.predictors)

X.reduced <- as.matrix(reduceX(coefficients,X,control$index1, control$index2, control$index3))
X.reduction <- reductionX(coefficients,X,control$index1, control$index2, control$index3)
beta.reduced <- as.matrix(reduceBeta(coefficients,control$index1, control$index2, control$index3))
try(beta.refitted <- suppressWarnings(as.matrix(round(glm.fit(X.reduced,y, weights, family=family,
    intercept = FALSE)$coefficients, digits=control$accuracy))), silent=TRUE)
if(!exists("beta.refitted")){beta.refitted <- as.matrix(rep(NA,
    times=length(beta.reduced)))}
if (rank == 0) rank <- dim(X.reduced)[2]

if (family$family=="gaussian")
    control$psi <- 1/(n-length(beta.reduced))*sum((y-mudach)^2)
aic <- -2*control$l(y,mudach) + 2*(length(beta.reduced)+(family$family=="gaussian"))

null <- glm(y~1, weights=weights, family=family, x=TRUE)
beta.null <- null$coefficients
X.null <- null$x

unless.null <- function(x, if.null) if (is.null(x))
    if.null else x
valideta <- unless.null(family$valideta, function(eta) TRUE)
validmu <- unless.null(family$validmu, function(mu) TRUE)
if (!(valideta(linear.predictors) && validmu(mudach))) 
    boundary <- TRUE  
    else boundary <- FALSE

# output
output <- list(
    coefficients = coefficients,
    coefficients.reduced = beta.reduced,
    coefficients.refitted = beta.refitted,
    coefficients.oml = control$oml,

    residuals = dev.res(y = y, mudach = mudach, weights = weights),
    fitted.values = mudach,
    rank = rank,
    family = family,
    linear.predictors = linear.predictors,
    deviance = round(dev(y=y,mudach=mudach,weights=weights),digits=2),
    aic = aic,
    null.deviance = round(dev(y=y, mudach=family$linkinv(X.null%*%beta.null),weights=weights),2),
    iter = iter,
    weights = weights, prior.weights = NULL,
    df.residual = df.residual,
    df.null = n-1,
    converged = converged,
    boundary = boundary,
    offset = NULL,
    control = control,
    method = method,
    contrasts = options("contrasts"),
    na.action  = "na.omit",
    plot = plot,
    tuning = tuning,
    index1 = control$index1,
    index2 = control$index2,
    index3 = control$index3, 
    number.selectable.parameters = control$number.selectable.parameters,
    number.removed.parameters = dim(X)[2]-dim(X.reduced)[2],
    x.reduction = X.reduction
    )

return(output)

}

