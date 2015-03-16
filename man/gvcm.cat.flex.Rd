\name{gvcm.cat.flex}
\alias{gvcm.cat.flex}
\title{Regularized Effects with Flexible Smoothing Parameters}
\description{
The function fits the same models with the same approximation as in \code{\link{gvcm.cat}} but the choice of the tuning parameter lambda for the penalty differs: instead of weighting the penalty terms and choosing on global tuning parmeter based on (generalized) cross-validation methods that again rely on the converged model, \code{\link{gvcm.cat.flex}} estimates several penalty parameteres lambda_i by linking the local quadratic approximation of \code{\link{gvcm.cat}} with the fantastic methods implemented in the package \code{\link{mgcv}}. This is why the arguments of \code{\link{gvcm.cat}} and \code{\link{gvcm.cat.flex}} differ. \code{\link{gvcm.cat.flex}} is not as well-developed as \code{\link{gvcm.cat}}.
}
\usage{
gvcm.cat.flex(whichCoefs, intercept = TRUE, data, family = gaussian(), method = "REML",
tuning = NULL, indexNrCoefs, indexPenNorm, indexPenA, indexPenWeight, 
control = list(c=1e-05, epsilon=1e-07, gama=35, maxi=1500, nu=.5))
}
\arguments{
  \item{whichCoefs}{vector with covariates (as characters)}
  \item{intercept}{logical}
  \item{data}{a data frame, with named and coded covariates}
  \item{family}{a \code{family} object describing the error distribution and link function to be used in the model; see \code{\link{family}} for details; everyl family that is compatible with \code{\link{gam}} is working}
  \item{method}{see \code{\link{gam}}}
  \item{tuning}{for function \code{\link{gam}}: argument \code{sp}}
  \item{indexNrCoefs}{vector with number of coefficients per covariate}
  \item{indexPenNorm}{vector with norm of the employed penalty (as.character)}
  \item{indexPenA}{list with the penalty matrices A_j for each covariate j}
  \item{indexPenWeight}{list, possible weights for the penalty terms (each entry is a vector)}
  \item{control}{a list of parameters for controlling the fitting process; must be \code{NULL} or contain all named elements}
}
\details{
The local quadratic approximation are linked to the methods of \code{\link{mgcv}} by alternating the update of the penalty and the update of the PIRLS algorithm/estimating the tuning parameters lambda_i via \code{\link{mgcv}}. Therefore, \code{\link{gvcm.cat.flex}} can be slow (but will be faster than \code{\link{gvcm.cat}} for the most part).
}
\value{
A \code{\link{gamObject}}.
}
\seealso{
Function \code{\link{gvcm.cat}}.
}
\examples{
\dontrun{
# compare gvcm.cat.flex and gvcm.cat for Lasso-type penalties:
n <- 100
ncov <- 7
set.seed(123)
X <- matrix(rnorm(n*ncov, sd=5), ncol=ncov)
coefs <- rpois(ncov + 1, 2)
y <- cbind(1, X) %*% coefs + rnorm(n, 1.5)
data <- as.data.frame(cbind(y, X))
names(data) <- c("y", paste("x", 1:ncov, sep=""))

m1 <- gvcm.cat.flex(
whichCoefs = paste("x", 1:ncov, sep=""),
data=data,
indexNrCoefs=rep(1, ncov),
indexPenNorm=rep("L1", ncov),
indexPenA=list(1,1,1,1,1,1,1),
indexPenWeight=list(1,1,1,1,1,1,1)
)

m2 <- gvcm.cat(y ~ 1 + p(x1) + p(x2) + p(x3) + p(x4) + p(x5) + p(x6) + p(x7),
              data=data, tuning=list(lambda=m1$sp, specific=TRUE), start=rep(1, 8))

rbind(m1$coefficients, m2$coefficients)

# Lasso-type fusion penalty with gvcm.cat.flex
n <- 100
ncat <- 8
set.seed(567)
X <- t(rmultinom(n, 1, rep(1/ncat, ncat)))[, -1]
coefs <- c(rpois(1, 2), sort(rpois(ncat-1, 1)))

y <- cbind(1, X) %*% coefs + rnorm(n, 1.5)
data <- as.data.frame(y)
data$x1 <- X
names(data) <- c("y", "x1")

A <- a(1:(ncat-1), ncat-2)

m3 <- gvcm.cat.flex(
whichCoefs = c("x1"),
data = data,
indexNrCoefs = c(ncat-1),
indexPenNorm = c("L1"),
indexPenA = list(A),
indexPenWeight = list(rep(1, ncol(A))),
tuning = 100 # fixed and large - in order to demonstrate the fusion of the coefficients
)
m3$coefficients
}
}

