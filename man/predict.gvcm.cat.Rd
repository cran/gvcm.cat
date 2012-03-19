\name{predict.gvcm.cat}
\alias{predict.gvcm.cat}
\title{Predict Method for gvcm.cat Fits}
\description{
Obtains predictions from a fitted \code{gvcm.cat} object. 
}
\usage{
\method{predict}{gvcm.cat}(object, newdata, type = "link", ...)
}
\arguments{
  \item{object}{a fitted object of class \code{\link{gvcm.cat}}}
  \item{newdata}{a data frame in which to look for variables with which to predict}
  \item{type}{the type of prediction required. The default is on the scale of the linear predictors; the alternative \code{"response"} is on the scale of the response variable. Thus for a binomial model the default predictions are of log-odds (probabilities on logit scale) and \code{type = "response"} gives the predicted probabilities}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{Observations containing \code{\link{NA}}s are always omitted.}
\value{
  \item{fit}{predictions}
  \item{fit.refitted}{predictions assuming refitted coefficients}
  \item{fit.oml}{predictions assuming maximum likelihood estimates}    
  \item{na.action}{information returned by \code{\link{model.frame}} on the special handling of \code{\link{NA}}s; currently always \code{\link{na.omit}}}
}
\seealso{Function \code{\link{gvcm.cat}}}
\examples{
## continues example of function gvcm.cat 
\donttest{newdata <- simulation(200, covariates, NULL, true.f, true.coefs , binomial(), seed=789)
prediction <- predict.gvcm.cat(m1, newdata) }
}
\keyword{gvcm.cat}
