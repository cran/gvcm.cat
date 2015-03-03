\name{index}
\alias{index}
\alias{design}
\title{Functions to Build Design Matrices and Indices for Function gvcm.cat()}
\description{
\code{design()} builds design matrices for function \code{\link{gvcm.cat}}; \code{index()} computes indices with information about the terms of the formula. 
}
\usage{
design(formula, data)

index(dsgn, data = data, formula = formula)

}
\arguments{
  \item{formula}{an object of class \code{"formula"}; see \code{\link{gvcm.cat}}}
  \item{data}{a data frame; see \code{\link{gvcm.cat}}}
  \item{dsgn}{value of function \code{design()}}
}
\value{
  \item{X}{the model matrix}
  \item{Terms}{the according \code{\link{terms.object}}}
  \item{m}{the model frame}
  \item{int}{either 0, indicating that the intercept is varying, or 1 indicating that the intercept is constant}
  \item{formula}{sorted version of the given formula, \code{index} vectors will refer to this formula}
  \item{a matrix}{value of function \code{index}}
}
\details{
Function \code{index} returns a matrix with one indicator vector per row. The columns refer to the elements of the formula (same order).
The indicator/indices are:
	      \itemize{
	        \item  \code{index1} : gives the number of coefficients belonging to each term in the formula. 
                  An entry is 1 if the according term is metric, it equals the number of the coded variable's categories, 
                  if the variable is a factor. If a continuous variable is modified by a factor \code{u} the entry equals 
                  the number of \code{u}'s categories
	        \item  \code{index2} : indicates varying coefficients. An entry is 0 if the according coefficient is not varying, 
                     it is -1 if the according coefficient is nominal, 1 if it is ordinal
	        \item  \code{index2b}: conforms to indicator \eqn{b_j} in Oelker et. al. 2012
	        \item  \code{index3} : indicates penalized covariates \code{p(u)}. An entry is 0 if the according covariate is 
                     not penalized, it is -1 if the according covariate is nominal, 1 if it is ordinal or metric
	        \item  \code{index4} : indicates penalized covariates \code{grouped(u)}. An entry is 0 if the according covariate is 
                     not penalized, it is -1 if the according covariate is nominal, 1 if it is ordinal or metric
	        \item  \code{index5} : experimental
	        \item  \code{index6} : indicates penalized covariates \code{sp}
	        \item  \code{index7} : indicates penalized covariates \code{SCAD}. An entry is 0 if the according covariate is 
                     not penalized, it is -1 if the according covariate is nominal, 1 if it is ordinal or metric
	        \item  \code{index8} : indicates penalized covariates \code{elastic}. An entry is 0 if the according covariate is 
                     not penalized, it is -1 if the according covariate is nominal, 1 if it is ordinal or metric
	        \item  \code{index9} : experimental
	      }

}
\references{
Oelker, M.-R., J. Gertheiss and G. Tutz (2012). Regularization and model melection with categorial predictors and effect modifiers in generalized linear models. \emph{Department of Statistics at the University of Munich: Technical Report 122}.\cr
}
\seealso{
Functions \code{\link{pest}}, \code{\link{abc}}
}
\examples{
## example for function simulation()
covariates <- list(x1=list("unif", c(0,2)),
                  x2=list("unif", c(0,2)),
                  x3=list("unif", c(0,2)),
                  u=list("multinom",c(0.3,0.4,0.3), "nominal")
                  )
true.f <- y ~ 1 + v(x1,u) + x2
true.coefs <- c(0.2,  0.3,.7,.7, -.5)
data <- simulation(400, covariates, NULL, true.f, true.coefs , binomial(), seed=456)
## example for function index()
f <- y ~ v(1,u) + v(x1,u) + v(x2,u)
dsgn <- design(f, data)
index(dsgn, data)
}
\keyword{gvcm.cat}
