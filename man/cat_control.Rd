\name{cat_control}
\alias{cat_control}
\title{Auxiliary Function for gvcm.cat Fitting}
\description{
Auxiliary function for \code{\link{gvcm.cat}} fitting. Used to modify the algorithm's internal parameters.
}
\usage{
cat_control(standardize = FALSE, accuracy = 2, digits = 4, initials = NULL, 
g = 0.5, epsilon = 10^(-5), maxi = 250, steps = 25, c = 10^(-5), 
tuning.criterion = "deviance", K = 5, cv.refit = FALSE, lambda.upper = 10000, 
lambda.accuracy = 0.01, adapted.weights = TRUE, assured.intercept = TRUE, 
...)
}
\arguments{
  \item{standardize}{logical; if \code{TRUE} all metric covariates are standardized by their empirical variances}
  \item{accuracy}{integer; number of digits being compared when setting coefficents equal/to zero; must be 4 or less}
  \item{digits}{integer; number of digits for estimates}
  \item{initials}{numeric vector; starting parameter values for penalized estimation}
  \item{g}{step length parameter for the \acronym{PIRLS}-algorithm, see Ulbricht (2010); out of )0,1(}
  \item{epsilon}{small, positive, real constant; termination criterion for the \acronym{PIRLS}-algorithm, see Ulbricht (2010)}
  \item{maxi}{integer; maximal number of iterations in the fitting algorithm}
  \item{steps}{integer; tuning parameter for path-plotting; minimal number of estimates employed for path-plotting}
  \item{c}{small, positive, real constant; needed for a approximation of the absolute value function in the \acronym{PIRLS}-algorithm, see Ulbricht (2010)}
  \item{tuning.criterion}{loss criterion for cross-validation; one out of \code{"deviance"}, \code{"SSE"} (for sum of squared errors)}
  \item{K}{integer; folds for cross-validation}
  \item{cv.refit}{logical; if \code{TRUE} cross-validation is based on a refit of selected coefficients}
  \item{lambda.upper}{integer; upper bound for cross-validation of \code{lambda}}
  \item{lambda.accuracy}{numeric; how accurate shall \code{lambda} be cross-validated?; minimal absolute difference between two candidates for \code{lambda}}
  \item{adapted.weights}{logical; if \code{TRUE}, penalty terms are weighted adaptively, that is by inverse ML-estimates; set to \code{FALSE}, if ML-estimates do not exist/are to close to zero}
  \item{assured.intercept}{logical; shall a constant intercept remain in the model in any case?}
  \item{\dots}{further arguments passed to or from other methods}  
}
\value{Returns a list containing the (checked) input arguments.}
\references{
Ulbricht, J. (2010). \emph{Variable Selection in Generalized Linear Models}. Dissertation an der Fakultaet fuer Mathematik, Informatik und Statistik an der Ludwig-Maximilians-Universitaet Muenchen: Verlag Dr. Hut.
}

\seealso{
Function \code{\link{gvcm.cat}}
}
\keyword{gvcm.cat}
