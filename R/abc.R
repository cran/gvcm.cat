abc <-
function(
X,
y,
indices,
family = gaussian,
method = c("AIC", "BIC"),
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
    stop("'family' not recognized. \n")
}

method <- match.arg(method)
if (!(method %in% c("AIC", "BIC")))
     stop ("method is incorrect. \n")

if (missing(control))
    control <- cat_control(...)
if (missing(weights))
    weights <- rep(1, times=dim(X)[1])
if (length(weights)!=nrow(X) || !is.vector(weights) || !is.numeric(weights))
    stop("Error in input weights. ") 
if (!is.logical(plot) || !is.matrix(X) || !is.numeric(X) || !is.vector(y) || !is.numeric(y) || nrow(X)!=length(y) || !is.list(indices))
     stop ("Error in input arguments. \n")
    
# definitions
n <- nrow(X)
if (!is.list(control$K) && control$K > n)
    ("K must be a single integer < dim(data)[1]. \n")
n.sp <- sum(abs(indices[[1]]*(indices[[2]]+indices[[3]])))- 
  (as.integer(control$assured.intercept))*abs(indices[[2]][1])

# functions
if(family$family == "binomial"){
  l <- function(y,mudach,weights=weights){sum((y*log(mudach) + (1-y)*log(1-mudach))*weights + lchoose(weights,y*weights))}
  d <- function(y,mudach,weights=weights){e <- matrix(0, nrow=length(y), ncol=1)
   rein.binaere <- if (any(y==0) || any(y==1)) c(which(y==0),which(y==1)) else -(1:length(y))
   e[rein.binaere,] <- (-2*log(1-abs(y-mudach))*weights)[rein.binaere]
   e[-rein.binaere,] <- (weights*(y*log(y/mudach)+(1-y)*log((1-y)/(1-mudach))))[-rein.binaere]
   return(e)}
  }
if(family$family == "gaussian"){
  l <- function(y,mudach,weights=weights){-1/2 * sum(weights*(y - mudach)^2) -     log(sqrt(2*pi))}
# l <- function(y,mudach){-1/2 * sum(weights*(y - mudach)^2)/psi - log(sqrt(2*pi*psi))}
  d <- function(y,mudach,weights=weights){weights*(y-mudach)^2}   # 
  }
if(family$family == "poisson"){
  l <- function(y,mudach,weights=weights){sum((y*log(mudach)) - mudach - lgamma(y+1))}
  d <- function(y,mudach,weights=weights){e <- matrix(0, nrow=length(y), ncol=1)
   e[which(y==0),] <-(2*mudach*weights)[which(y==0)]
   e[which(y!=0),] <- (2*weights*((y*log(y/mudach))+mudach-y))[which(y!=0)]
   return(e)}
  }
dev.res <- function(y,mudach,weights) {(y-mudach)/abs(y-mudach) * sqrt(d(y, mudach, weights)) }
dev <- function(y,mudach,weights) {sum(d(y=y,mudach=mudach, weights))}

# oml
suppressWarnings(try(oml.model <- glm.fit(X, y, weights = weights,
    family = family, intercept = FALSE), silent = TRUE))
if(exists("oml.model")==FALSE) {
    oml.model <- list(coefficients=rep(NA, times=dim(X)[2]), rank=NA, aic=NA,
      iter=NA, df.residual=NA, converged=FALSE)
    warning("Ordinary maximum likelihood estimate does not exsist. \n")
    control$adapted.weights <- FALSE
    }
oml <- as.matrix(round(oml.model$coefficients, control$digits))
if(sum(as.integer(is.na(oml)))>0) {
    control$adapted.weights <- FALSE
    warning("Ordinary maximum likelihood estimate contains NAs. \n")
    }
if(length(which(abs(oml)<.0001))>0) {
    control$adapted.weights <- FALSE
    warning("control$adapted.weights set to FALSE as at least one ML-estimate is too close to zero. \n")
    }

# model selection
if (method == "AIC") {criterion <- function(beta){0} } else
                     {criterion <- function(beta){(log(n)-2)*length(beta)} } 

if (n.sp==0){
    coefficients <- round(oml, control$accuracy)
    tuning <- round(oml.model$aic + criterion(coefficients), digits=2)
    rank <- oml.model$rank
    iter <- oml.model$iter
    df.residual <- oml.model$df.residual
    converged <- oml.model$converged
    plot <- list(NA,NA)
    }

if (n.sp>0){
    
    suppressWarnings(try(opt <- abc.optimierung(method, X, y, family, weights, control, indices, oml)))
    if(exists("opt")==FALSE) {stop("Variable selection via AIC/BIC failed. \n")}
    best <- (opt$model)
    beta.reduced <- matrix(round(best$coefficients,digits=control$accuracy),ncol=1)
    rownames(beta.reduced) <- colnames(opt$A.model)
    coefficients <- round(opt$A.model %*% beta.reduced, digits=control$accuracy)
    
    tuning <- round(opt$abc.model, digits=2)
    rank <- best$rank # numeric rank of the fitted linear model, sum(diag(H.i))
    iter <- opt$iter
    df.residual <- best$df.residual
    converged <- best$converged

    # plot
    if (plot==TRUE) {plot <- opt$A.models} else {plot <- NA}
    
    # aic 
    aic <- best$aic
    
    }


# prepare output
linear.predictors <- X %*% coefficients
mudach <- family$linkinv(linear.predictors)
residuals <- dev.res(y = y, mudach = mudach, weights = weights)
deviance <- round(dev(y=y,mudach=mudach,weights=weights),digits=2)

reduction <- reduce(coefficients, indices, control$assured.intercept)
X.reduced <- as.matrix(X %*% reduction$C)
X.reduction <- reduction$C
beta.reduced <- as.matrix(reduction$beta)
beta.refitted <- beta.reduced

aic <- -2*l(y,mudach,weights) + 2*( rank + (family$family=="gaussian"))

null <- glm(y~1, weights=weights, family=family, x=TRUE)
beta.null <- null$coefficients
X.null <- null$x
null.deviance <- round(dev(y=y, mudach=family$linkinv(X.null%*%beta.null),weights=weights),2)

unless.null <- function(x, if.null) if (is.null(x))
    if.null else x
valideta <- unless.null(family$valideta, function(eta) TRUE)
validmu <- unless.null(family$validmu, function(mu) TRUE)
if (!(valideta(linear.predictors) && validmu(mudach))) 
    boundary <- TRUE else 
    boundary <- FALSE


# output
output <- list(
    coefficients = coefficients,
    coefficients.reduced = beta.reduced,
    coefficients.refitted = beta.refitted,
    coefficients.oml = oml,

    residuals = residuals,
    fitted.values = mudach,
    rank = rank,
    family = family,
    linear.predictors = linear.predictors,
    deviance = deviance,
    aic = aic,
    null.deviance = null.deviance,
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
    indices = indices,
    number.selectable.parameters = n.sp,
    number.removed.parameters = dim(X)[2]-dim(X.reduced)[2],
    x.reduction = X.reduction
    )

return(output)
   
}

