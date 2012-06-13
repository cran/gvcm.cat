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
        stop("A formula argument is required. \n")
    if (indx[2] == 0)
        stop("A data argument is required. \n")
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
#       data[,-no] <- scale(data[,-no], center = FALSE, scale = TRUE)
       data[,-no] <- scale(data[,-no], center = FALSE, scale = apply(data[,-no],2,sd,na.rm=TRUE))
    }
   
# model.matrix
    dsgn <- design(formula,data)
    X <- dsgn$X
    
# response
    Y <- y
    y <- model.extract(dsgn$m, "response")
    if (is.factor(y)==TRUE){y <- as.numeric(y)-1}
    if (missing(weights))
        weights <- rep(1, times=dim(X)[1])
    if (length(weights)!=nrow(X) || !is.vector(weights) || !is.numeric(weights))
        stop("Error in input weights. ") 
    if (!is.null(dim(y)[2]) && family$family=="binomial") {
        weights <- (y[,1]+y[,2])*weights 
        y <- y[,1]/(y[,1]+y[,2])
        } 
    if (family$family=="binomial" && (sum(y>1) || sum(y<0))) 
        stop("No binomial response. \n") 

# definitions
    indices <- index(dsgn,data)

# default method        
    if (method %in% c("AIC", "BIC")) {
        output <-  abc(X, y, indices, family, method, weights, control, plot)
        } else {
        output <- pest(X, y, indices, family, method, tuning, weights, control, plot)
        }
    if (!exists("output"))
        stop("Error in argument 'method'. \n")     
                                
# output
    output$call <- Call
    output$formula <- dsgn$formula
    output$terms <- dsgn$Terms
    output$data <- data 
    output$x <- if(x==TRUE) X else NULL
    output$y <- if(Y==TRUE) y else NULL
    output$model <- if(model==TRUE) dsgn$m else NULL
    output$xlevels <- .getXlevels(dsgn$Terms, dsgn$m)
    class(output) <- c("gvcm.cat", "glm", "lm")
    output

}

