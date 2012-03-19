simulation <-
function(
n,
covariates, 
correlation = NULL,
formula, 
coefficients, 
family, 
sd=1,
seed=rpois(1,2348)*rnorm(1)
)

{
check.simulation(n, covariates, correlation = NULL, formula, 
coefficients, family, sd, seed)

# definitions
set.seed(seed)
m <- length(covariates)
data <- data.frame(ncol=0) 
namen <- c("ncol")
mu <- c()
sigma <- c()

# variables
for (i in 1:m) {  
    distribution <- covariates[[i]][[1]]
    p <- covariates[[i]][[2]]

    if (distribution %in% c("exp", "gamma", "pois")) 
      {x <- eval (parse( text = paste("r",distribution,"(",n,",",p[1],")", sep="") ))}
    if (distribution %in% c("lnorm", "beta", "unif")) 
      {x <- eval (parse( text = paste("r",distribution,"(",n,",",p[1], ",",p[2], ")", sep="") ))}
    if ( distribution %in% c("norm") ) { 
       mu <- c(mu, p[1])
       names(mu)[length(mu)] <- names(covariates)[i]
       sigma <- c(sigma, p[2])
       names(sigma)[length(sigma)] <- names(covariates)[i]
      } 
    if (distribution == "multinom")
      {x <- as.factor( t( t(1:length(p)) %*% rmultinom(n, 1, as.vector(p)) ))
      if (covariates[[i]][[3]] == "ordinal") {x <- as.ordered(x)}  }
                                                                                                                                                      
    if (distribution %in% c("exp", "gamma", "pois", "lnorm", "beta", "multinom", "unif")) 
      {data <- cbind(data, x)
       namen <- c(namen, names(covariates)[i])
       names(data) <- namen
       } 
    }

if (length(mu)>0) {
    if (is.null(correlation)==TRUE ) Sigma <- diag(sigma)
    if (is.matrix(correlation)==TRUE) Sigma <- diag(sigma) %*% correlation %*% diag(sigma)
    x <- mvrnorm(n=n, mu=mu, Sigma=Sigma, tol = 1e-6, empirical = FALSE)
    data <- cbind(data, x) 
    namen <- c(namen, names(mu))
    names(data) <- namen
    }

if (ncol(data)!=2){
    data <- data[,-1]
    r <- c()
    for (i in 1:length(names(covariates))){
        r <- c(r, which(names(data)==names(covariates)[i]))
        }
    data <- data[,r] 
    }else{
    data <- data.frame(data[,-1])
    colnames(data) <- namen[2]
    }

# model.matrix
special <- c("v", "p")
int <- if (grepl("v\\(", strsplit(deparse(formula[3]), "\\+")[[1]][1]) )
    0 else 1
m <- model.frame(formula=terms(formula[c(1,3)], specials=special, data=data),data)
Terms <- attr(m, "terms")
attr(Terms, "intercept") <- 1
options(contrasts = c("contr.effect", "contr.effect"))
X <- model.matrix(Terms, m)
if (int==0) X <- X[,-1]
namen <- colnames(X)
namen <- gsub(" ", "", namen, fixed=TRUE)
namen <- sub("v(", "", namen, fixed=TRUE)
namen <- sub("p(", "", namen, fixed=TRUE)
namen <- gsub("(", "", namen, fixed=TRUE)
namen <- gsub(")", "", namen, fixed=TRUE)
namen <- sub(",", ".", namen, fixed=TRUE)
colnames(X) <- namen

# response
link <- family$linkinv
E <- as.vector(link(as.matrix(X) %*% coefficients ))
if (family$family == "gaussian") {y <- rnorm(n,mean=E, sd= sd)}
if (family$family == "binomial") {y <- (rbinom(n, 1, E))}
if (family$family == "poisson")  {y <- rpois(n, E)} 
if (sum(as.integer(is.na(y))) > 0)
    stop ("Response contains NAs. \n")

#return
daten <- cbind(y,data)
names(daten) <- c("y",names(data))
return(daten)
}

