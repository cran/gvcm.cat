summary.gvcm.cat <-
function (
object, ...
)
{
# check input
if (("gvcm.cat" %in% is(object))==FALSE )
     stop ("object must be a 'gvcm.cat' object. \n")
       
# defintions
dev.res <- strsplit(summary(round(object$residuals,3)), ":")
index.reduced <- c()
for (i in 1:dim(object$x.reduction)[2]){
   index.reduced <- c(index.reduced,min(which(object$x.reduction[,i]==1))) 
   }

coefs <- data.frame(object$coefficients)
coefs[,1]<- object$coefficients.oml
coefs[,2]<- object$coefficients
coefs[index.reduced,3]<- object$coefficients.reduced
coefs[index.reduced,4]<- object$coefficients.refitted
colnames(coefs) <- c("coefficients.oml", "coefficients", "coefficients.reduced", 
   "coefficients.refitted" )

# summary
cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
    "\n\n", sep = "")

cat("Deviance Residuals: \n")
cat("    Min        1Q    Median        3Q       Max \n")
cat(dev.res[[1]][2],dev.res[[2]][2],dev.res[[3]][2],dev.res[[5]][2],
    dev.res[[6]][2],"\n \n")

cat("Coefficients: \n")
print(coefs)
cat("\n")

cat("(Dispersion parameter for ", object$family$family," family taken to be ", 
    object$control$psi, ") \n", sep="")

cat("    Null deviance: ", object$null.deviance," on ", object$df.null,
    " degrees of freedom \n", sep="")
cat("Residual deviance: ", object$deviance," on ", round(object$df.residual, 2),
    " degrees of freedom \n \n", sep="")

cat("Removed parameters: ", object$number.removed.parameters, " out of ", 
    object$number.selectable.parameters, "\n", sep="")
if(object$method %in% c("nlm","lqa")){
cat("Penalization parameter lambda = ", object$tuning[[1]], "\n", sep="")
cat("Weighting parameters: phi = ", object$tuning[[2]], ", adapted.weights = ", 
    object$control$adapted.weights, ", assured.intercept = ", 
    object$control$assured.intercept, "\n", sep="")
cat("Minimal lambda causing maximal penalization (", object$control$accuracy,
    " digits): ",object$control$lambda.upper, "\n \n", sep="")
}
if(object$method %in% c("AIC")){
cat("AIC of chosen model: ", object$tuning, "\n", sep="")
}
if(object$method %in% c("BIC")){
cat("BIC of chosen model: ", object$tuning, "\n", sep="")
}
cat("Number of iterations: ", object$iter, "\n", sep="")
if (object$converged==TRUE) {cat("The model converged. \n", sep="")} else 
   {cat("The model did not converge. \n", sep="")}

}

