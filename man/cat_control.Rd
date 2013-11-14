\name{cat_control}
\alias{cat_control}
\title{Auxiliary Function for gvcm.cat}
\description{
Auxiliary function for \code{\link{gvcm.cat}}. Modifies the algorithm's internal parameters.
}
\usage{
cat_control(center = FALSE, standardize = FALSE, accuracy = 2, digits = 4,
g = 0.5, epsilon = 10^(-5), maxi = 250, c = 10^(-5), gama = 20, steps = 25, 
nu = 1, tuning.criterion = "GCV", K = 5, cv.refit = FALSE, 
lambda.upper=50, lambda.lower=0, lambda.accuracy=.01, scaled.lik=FALSE,
adapted.weights=FALSE, adapted.weights.adj = FALSE, adapted.weights.ridge =
FALSE, assured.intercept=TRUE, 
level.control = FALSE, case.control = FALSE, pairwise = TRUE, 
grouped.cat.diffs = FALSE, bootstrap = 0, start.ml = FALSE, L0.log = TRUE, 
subjspec.gr = FALSE, high = NULL, ...)
}
\arguments{
  \item{center}{logical; if \code{TRUE}, all metric covariates are centered by their empirical mean}
  \item{standardize}{logical; if \code{TRUE}, the design matrix is standardized by its (weighted) empirical variances}
  \item{accuracy}{integer; number of digits being compared when setting coefficents equal/to zero}
  \item{digits}{integer; number of digits for estimates}
  \item{g}{step length parameter for the \acronym{PIRLS}-algorithm; out of )0,1(}
  \item{epsilon}{small, positive, real constant; the \acronym{PIRLS}-algorithm is terminated when the (scaled, absolute) difference of the coefficients of the current iteration and the coefficients of the previous iteration is smaller than \code{epsilon}}
  \item{maxi}{integer; maximal number of iterations in the fitting algorithm}
  \item{c}{small, positive, real constant; needed for the approximation of the absolute value function in the \acronym{PIRLS}-algorithm}
  \item{gama}{positive number; tuning parameter for the approximation of the L0 norm}
  \item{steps}{integer; tuning parameter for path-plotting; minimal number of estimates employed for path-plotting}
  \item{nu}{optional weighting parameter}
  \item{tuning.criterion}{loss criterion for cross-validation; one out of \code{"GCV"} (generalized cross validation criterion), \code{"deviance"} (K-fold cross-validation with the predictive deviance as criterion)}
  \item{K}{integer; number of folds for cross-validation}
  \item{cv.refit}{logical; if \code{TRUE}, cross-validation is based on a refit of the selected coefficients} 
  \item{lambda.upper}{integer; upper bound for cross-validation of \code{lambda}}
  \item{lambda.lower}{integer; lower bound for cross-validation of \code{lambda}}
  \item{lambda.accuracy}{numeric; how accurate shall \code{lambda} be cross-validated?; minimal absolute difference between two candidates for \code{lambda}}
  \item{scaled.lik}{if \code{TRUE}, the likelihood in the objective function is scaled by 1/n}
  \item{adapted.weights}{logical; if \code{TRUE}, penalty terms are weighted adaptively, that is by inverse ML-estimates; set to \code{FALSE}, if ML-estimates do not exist/are to close to zero; only for specials v, p, grouped, SCAD, elastic}
  \item{adapted.weights.adj}{logical; if \code{TRUE}, adapted weights of several categorical covariates are scaled such that they are comparable}
  \item{adapted.weights.ridge}{logical; if \code{TRUE}, adapted weights are based on aa estimate that is slightly penalized by a Ridge penalty}
  \item{assured.intercept}{logical; shall a constant intercept remain in the model in any case?}
  \item{level.control}{logical; if \code{TRUE}, the penalty terms are adjusted for different number of penalty terms per covariate}
  \item{case.control}{logical; if \code{TRUE}, the penalty terms are adjusted for the number of observations on each level of a categorical covariate}
  \item{pairwise}{experimental option; disabled if \code{TRUE}}
  \item{grouped.cat.diffs}{experimental option; disabled if \code{FALSE}}
  \item{bootstrap}{experimental option; disabled if \code{0}}
  \item{start.ml}{logical; if \code{TRUE}, the initial value is the ML-estimate}
  \item{L0.log}{experimental option; disabled if \code{TRUE}}  
  \item{subjspec.gr}{experimental option; disabled if \code{FALSE}}
  \item{high}{experimental option; disabled if \code{NULL}}
  \item{\dots}{further arguments passed to or from other methods}  
}
\value{Returns a list containing the (checked) input arguments.}

\seealso{
Function \code{\link{gvcm.cat}}
}
\keyword{gvcm.cat}
