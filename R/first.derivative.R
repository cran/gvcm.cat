first.derivative <-
function(beta, lambda, weight, acoefs, indices, control) 
{ 
if (is.null (weight) && sum(indices[[2]] + indices[[3]])!=0) 
{stop ("'weight' must be the current weight vector \n")}
p <- length (weight)
return (rep (lambda, p) * weight * as.integer (drop(t(acoefs)%*%beta) != 0))
 
}

