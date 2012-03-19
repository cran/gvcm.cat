first.derivative <-
function(beta, lambda, weight, acoefs, control) 
{ 
if (is.null (weight) && sum(control$index2 + control$index3)!=0) 
{stop ("'weight' must be the current weight vector \n")}

p <- length (weight)

if (p < 2) { stop ("There must be at least two regressors! \n") }

return (rep (lambda, p) * weight * as.integer (drop(t(acoefs)%*%beta) != 0))
 
}

