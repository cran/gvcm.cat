gvcm.cat.default <-
function(
formula,
data,
family = gaussian,
method = "lqa", 
tuning = list(lambda=TRUE, phi=0.5),
weights,
control,
model = FALSE,
x = FALSE,
y = FALSE,
plot=FALSE,
...
)

{

# check
    Call <- match.call()
    indx <- match(c("formula", "data"),
        names(Call), nomatch = 0)
    if (indx[1] == 0)
        stop("A formula argument is required")
    if (indx[2] == 0)
        stop("A data argument is required")
    if (missing(control))
        control <- cat_control(...)
    
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }

    if (!is.logical(model) || !is.logical(x) || !is.logical(y) || !is.logical(plot))
         stop ("Error in input arguments. \n")

# standardize + na remove
    data <- na.omit(data)

    if(control$standardize){
       no <- which(names(data)==formula[[2]])
       for (i in 1:dim(data)[2]) {
          if (is.factor(data[,i])) {no <- c(no,i)}
       }
       data[,-no] <- scale(data[,-no], center = FALSE, scale = TRUE)
    }
   
# model.frame
    special <- c("v","p")
    int <- if (grepl("v\\(", strsplit(deparse(formula[3]), "\\+")[[1]][1]) )
        0 else 1
    m <- model.frame(formula=terms(formula, specials=special, data=data), data)
    if (nrow(m) == 0)
        stop("No (non-missing) observations")

# model.matrix
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 1
#    options(contrasts = c("contr.effect", "contr.effect"))     
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
    Y <- y
    y <- model.extract(m, "response")
    if (is.factor(y)==TRUE){y <- as.numeric(y)-1}
    if (missing(weights))
        weights <- rep(1, times=dim(X)[1])
    if (!is.null(dim(y)[2]) && family$family=="binomial") {
        weights <- (y[,1]+y[,2])*weights 
        y <- y[,1]/(y[,1]+y[,2])
        } else {y <- weights*y}
    if (family$family=="binomial" && (sum(y>1) || sum(y<0))) 
        stop("No binomial response") 

# definitions
    ind <- index(formula,data)

# default method        
    if (method %in% c("lqa", "nlm"))
    output <- pest(X, y, ind, family, method, tuning, weights, control, plot)
    if (method %in% c("AIC", "BIC"))
    output <- abc(X, y, ind, family, method, weights, control, plot)
    if (!exists("output"))
        stop("Error in argument 'method'")     
                                
    output$call <- Call
    output$formula <- formula
    output$terms <- Terms
    output$data <- data 
    output$x <- if(x==TRUE) X else NULL
    output$y <- if(Y==TRUE) y else NULL
    output$model <- if(model==TRUE) m else NULL
    output$xlevels <- .getXlevels(Terms, m)
    class(output) <- c("gvcm.cat", "glm", "lm")
    output

}

