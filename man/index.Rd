\name{index}
\alias{index}
\title{Index Function for Function gvcm.cat()}
\description{
\code{index()} computes three statements per variable; one indicating how many coefficients belong to this variable (\code{index1}), the second indicating whether this variable varies or not (\code{index2}), the third indicating wheter the variable is categorial and penalized or not.
}
\usage{
index(formula, data)
}
\arguments{
  \item{formula}{an object of class \code{"formula"}; see \code{\link{gvcm.cat}}}
  \item{data}{a data frame; see \code{\link{gvcm.cat}}}
}
\value{
  \item{index1}{a vector; its length equals the number of the model's variables including the intercept, each entry gives the number of coefficients belonging to the according variable. An entry is 1 if the according variable is metric, it equals the number of the coded variable's categories, if the variable is a factor. If a continuous variable is modified by a factor \code{u} the entry equals the number of \code{u}'s categories}
  \item{index2}{a vector of the same length as \code{index1}; indicates varying coefficients. An entry is 0 if the according coefficient is not varying, it is -1 if the according coefficient is nominal, 1 if it is ordinal}
  \item{index3}{a vector of the same length as \code{index1}; indicates penalized categorial covariates. An entry is 0 if the according covariate is not categorial, it is -1 if the according covariate is nominal, 1 if it is ordinal}
}
\seealso{
Functions \code{\link{pest}}, \code{\link{abc}}
}
\examples{
## continues example of function simulation
\donttest{f <- y ~ v(1,u) + v(x1,u) + v(x2,u)
index(f, data)}
}
\keyword{gvcm.cat}
