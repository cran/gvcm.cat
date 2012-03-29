predict.gvcm.cat <-
function(
object,   
newdata, 
type="link", # "link", "response"
...
)
{
# check input
  if (("gvcm.cat" %in% is(object))==FALSE )
       stop ("object must be a 'gvcm.cat' object. \n")     
  if (("data.frame" %in% is(newdata))== FALSE)
       stop ("newdata must be a dataframe. \n")
  if ( (type %in% c("link", "response"))==FALSE )
       stop ("type must be one out of 'link', 'response'. \n")

# definitions out of given object
formula <- object$formula
control <- object$control

# na action
data <- na.omit(newdata)

# standardize??
if(control$standardize){
   no <- which(names(data)==formula[[2]])
   for (i in 1:dim(data)[2]) {
      if (is.factor(data[,i])) {no <- c(no,i)}
   }
   data[,-no] <- scale(data[,-no], center = FALSE, scale = TRUE)
}

# define X.new
special <- c("v", "p")
int <- if (grepl("v\\(", strsplit(deparse(formula[3]), "\\+")[[1]][1]) )
    0 else 1
m <- model.frame(formula=terms(formula, specials=special, data=data), data)
if (nrow(m) == 0)
    stop("No (non-missing) observations")
Terms <- attr(m, "terms")
attr(Terms, "intercept") <- 1
#options(contrasts = c("contr.effect", "contr.effect"))
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
X.reduced <- X %*% object$x.reduction

# response
y <- model.extract(m, "response")
if (is.factor(y)==TRUE){y <- as.numeric(y)-1}

# type
if (type == "link") {link <- function(x){x}} else # predictor level
                    {link <- object$family$linkinv}        # response level

# return
return(list (
fit = link(X %*% object$coefficients), 
fit.refitted = link(X.reduced %*% object$coefficients.refitted),
fit.oml = link(X %*% object$coefficients.oml),
na.action = "na.omit"))  

}

