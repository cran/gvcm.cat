\name{simulation}
\alias{simulation}
\title{Simulates data with categorial covariates}
\description{Simulates data with categorial covariates/categorial effect modifiers}
\usage{
simulation(n, covariates, correlation = NULL, formula, coefficients, 
family, sd = 1, seed = rpois(1, 2348) * rnorm(1))
}
\arguments{
  \item{n}{number of observations; must be large enough, so that all categories of all factor variables exist and therefore vector \code{coefficients} fits}
  \item{covariates}{description of the covariates and effect modifiers included in the model; format: \code{list( name of variable 1 = list("distribution", c(parameters), "level of measurement"), name of variable 2 = ...)}}
  \item{correlation}{optional matrix, specifies the correlation of Gaussian covariates}
  \item{formula}{formula like in \code{\link{gvcm.cat}} (all variables contained in \code{\link{formula}} must be defined in \code{covariates})}
  \item{coefficients}{true parameter vector}
  \item{family}{a \code{\link{family}} object; currently only \code{gaussian}, \code{binomial}, \code{poisson}}
  \item{sd}{if needed, i.e. if \code{family = gaussian}, standard deviation of response}
  \item{seed}{specifies the to be used seed}
}
\details{
Remarks on covariates:
\itemize{
	\item all parameterizations like default in \code{\link{Distributions}}.
	\item possible distributions of covariates (required as characters), their parameters (required as vectors) and constraints (in parentheses):
	      \itemize{
	        \item  \code{beta}    : shape1 (>0), shape2 (>0)
	        \item  \code{exp}     : rate (>0)
	        \item  \code{gamma}   : shape (>0)
	        \item  \code{lnorm}   : mean , sd (>0)
	        \item  \code{multinom}: vector of the categories' probabilities (all elements must be >0, sum over all elements must be 1)
	        \item  \code{norm}    : mean, sd (>0)
	        \item  \code{pois}    : lambda (>0)
	        \item  \code{unif}    : min, max
	      }
	\item \code{level of measurement} is only needed for \code{distribution = "multinom"}, must be \code{"nominal"} or \code{"ordinal"}.
	\item If any, the covariates' correlation is specified by argument \code{correlation}. 
	Correlations are defined for Gaussian covariates only. 
	Matrix \code{correlation} refers to these covariates according to the order they are listed in \code{covariates}. So that the dimensions of \code{correlation} must fit to the number of normal distributed variables in \code{covariates}. 
} 

}
\value{A data frame containing all specified covariates (even if they are not included in \code{\link{formula}}) and the response (named \code{y})}
\seealso{Function \code{\link{gvcm.cat}}}
\examples{
## example function simulation
covariates <- list(x1=list("unif", c(0,2)),
                  x2=list("unif", c(0,2)),
                  x3=list("unif", c(0,2)),
                  u=list("multinom",c(0.3,0.4,0.3), "nominal")
                  )
true.f <- y ~ 1 + v(x1,u) + x2
true.coefs <- c(0.2,  0.3,.7,.7, -.5)
data <- simulation(400, covariates, NULL, true.f, true.coefs , binomial(), seed=456)
}
\keyword{gvcm.cat}
