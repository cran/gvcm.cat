\name{index}
\alias{index}
\alias{design}
\title{Functions to Build Design Matrices and Indices for Function gvcm.cat()}
\description{
\code{design()} builds design matrices considering special formula terms \code{v()} and \code{p()}. 
\code{index()} computes several indices with different information about each term of the formula. 
}
\usage{
design(formula, data)

index(dsgn, data)

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
  \item{index1}{a vector; its length equals the number of the formula's terms including the intercept, each entry gives the number of coefficients belonging to the according term. An entry is 1 if the according term is metric, it equals the number of the coded variable's categories, if the variable is a factor. If a continuous variable is modified by a factor \code{u} the entry equals the number of \code{u}'s categories}
  \item{index2}{a vector of the same length as \code{index1}; indicates varying coefficients. An entry is 0 if the according coefficient is not varying, it is -1 if the according coefficient is nominal, 1 if it is ordinal}
  \item{index3}{a vector of the same length as \code{index1}; indicates penalized covariates \code{p(u)}. An entry is 0 if the according covariate is not penalized, it is -1 if the according covariate is nominal, 1 if it is ordinal or metric}
}
\seealso{
Functions \code{\link{pest}}, \code{\link{abc}}
}
\examples{
## continues example of function simulation
\donttest{f <- y ~ v(1,u) + v(x1,u) + v(x2,u)
dsgn <- design(f, data)
index(dsgn, data)}
}
\keyword{gvcm.cat}
