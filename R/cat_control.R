cat_control <-
function (
standardize = FALSE, 
accuracy = 2, digits = 4,
initials = NULL,
g = 0.5, epsilon = 10^(-5), maxi = 250, steps = 25, c = 10^(-5), 
tuning.criterion = "deviance", K = 5,
cv.refit = FALSE, 
lambda.upper=10000, lambda.accuracy=.01,
adapted.weights=TRUE,
assured.intercept=TRUE, 
...)

{
if (!is.logical(standardize))
     stop ("standardize must be logical. \n")

if (!is.numeric(accuracy) || accuracy < 0)
    stop("'accuracy' must be >= 0")
accuracy <- as.integer(accuracy)

if (!is.numeric(digits) || digits < 0)
    stop("'digits' must be >= 0")
digits <- as.integer(digits)

if (!is.null(initials)){
if (!is.vector(initials) || !is.numeric(initials))     
     stop ("initials must be a numeric vector of proper length. \n")
}

if (!is.numeric(g) || length(as.vector(g))!=1 || g>1 || g<0)
     stop ("g must be a singl, numeric value out of ]0;1[. \n")

if (!is.numeric(epsilon) || epsilon <= 0 || epsilon>1 || length(epsilon)!=1)
   stop ("epsilon is ment to be a single, small, positive and numeric value. \n")

if (!is.numeric(maxi) || length(maxi)!=1 || maxi < 1)
     stop ("maxi must be a sinlge integer > 0 \n")
maxi <- as.integer(maxi)

if (!is.numeric(steps) || length(steps)!=1 || steps < 1)
     stop ("steps must be a sinlge integer > 0 \n")
steps <- as.integer(steps)

if (!is.numeric(c) || c<0 || length(c)!=1)
    stop ("c is ment to be a single, small, positive and numeric value. \n")


if (!(tuning.criterion %in% c("deviance","SSE")))
     stop ("tuning.criterion must be one out of 'deviance' or 'SSE'. \n")

if (!is.vector(K))
         stop ("K must be a single integer > 1. \n")
if (!is.list(K)) {
    if (!is.numeric(K) || K<2)
         stop ("K must be a single integer > 1. \n")
    }

if (!is.logical(cv.refit))
     stop ("cv.refit must be logical. \n")

if (!is.numeric(lambda.upper) || length(lambda.upper)!=1)
     stop ("lambda.upper must be numeric and positive \n")

if (!is.numeric(lambda.accuracy) || lambda.accuracy <= 0 || length(lambda.accuracy)!=1)
   stop ("lambda.accuracy is ment to be a single, positive and numeric value. \n")

if (!is.logical(adapted.weights))
     stop ("adapted.weights must be logical. \n")

if (!is.logical(assured.intercept))
     stop ("assured.intercept must be logical. \n")


list(
standardize = standardize,
accuracy = accuracy, digits = digits,
initials = initials,
g = g, epsilon = epsilon, maxi = maxi, steps = steps, c = c,
tuning.criterion = tuning.criterion, K = K, cv.refit = cv.refit,
lambda.upper = lambda.upper, lambda.accuracy = lambda.accuracy , 
adapted.weights=adapted.weights, assured.intercept=assured.intercept
)

}

