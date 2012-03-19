abc <-
function(
X,
y,
ind,
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
    stop("'family' not recognized")
}

method <- match.arg(method)
if (!(method %in% c("AIC","BIC")))
     stop ("method is incorrect. \n")

# model
# plot

# definitions
control$index1 <- ind[[1]]
control$index2 <- ind[[2]]
control$index3 <- ind[[3]]
n <- nrow(X)
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
suppressWarnings(try(oml.model <- glm.fit(X, y, weights, offset = rep(0,n),
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

# model selection
if (control$number.selectable.parameters==0){
    coefficients <- control$oml
    tuning <- oml.model$aic 
    rank <- oml.model$rank
    iter <- oml.model$iter
    df.residual <- oml.model$df.residual
    converged <- oml.model$converged
    plot <- NA
    aic <- oml.model$aic
    }

if (control$number.selectable.parameters>0){
    
    suppressWarnings(try(opt <- abc.optimierung(method, X, y, family, weights, control)))
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

X.reduced <- as.matrix(reduceX(coefficients, X, control$index1, control$index2, control$index3))
X.reduction <- reductionX(coefficients, X, control$index1, control$index2, control$index3)
beta.reduced <- as.matrix(reduceBeta(coefficients,control$index1, control$index2, control$index3))
beta.refitted <- beta.reduced

null <- glm(y~1, weights=weights, family=family, x=TRUE)
beta.null <- null$coefficients
X.null <- null$x

unless.null <- function(x, if.null) if (is.null(x))
    if.null
    else x
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

    residuals = dev.res(y = y, mudach = mudach, weights=weights),
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

