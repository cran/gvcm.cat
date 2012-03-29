\name{gvcm.cat}
\alias{gvcm.cat}
\alias{gvcm.cat.default}
\alias{pest}
\alias{abc}
\alias{print.gvcm.cat}
\alias{summary.gvcm.cat}
\title{Regularized Categorial Effects/Categorial Effect Modifiers in GLMs}
\description{
The function fits generalized linear models with categorial effects and categorial effect modifiers. 
The model is specified by giving a symbolic description of the linear predictor and a description of the error distribution. 
Estimation employs regularization and model selection strategies to fuse and/or select a covariate's categories. 
These strategies are either a penalty combining the fused and the pure Lasso or a forward selection strategy employing AIC/BIC (see Oelker et. al. 2012).
}
\usage{
gvcm.cat(formula, data, family = gaussian, method = "lqa", 
tuning = list(lambda=TRUE, phi=0.5), weights, control, 
model = FALSE, x = FALSE, y = FALSE, plot=FALSE, ...)

\method{gvcm.cat}{default}(formula, data, family = gaussian, method = "lqa", 
tuning = list(lambda=TRUE, phi=0.5), weights, control, 
model = FALSE, x = FALSE, y = FALSE, plot=FALSE, ...)
 
pest(X, y, ind, family = gaussian, method = c("lqa", "nlm"), 
tuning = list(lambda=TRUE, phi=0.5), weights, 
control = cat_control(), plot=FALSE, ...)

abc(X, y, ind, family = gaussian, method = c("AIC", "BIC"), 
weights, control = cat_control(), plot=FALSE, ...)
}
\arguments{
  \item{formula}{an object of class \code{\link{formula}}: a symbolic description of the model to be fitted. See details}
  \item{data}{an optional data frame, containing the variables in the model.}
  \item{family}{a \code{family} object describing the error distribution and link function to be used in the model; this can be a character string naming a family function, a family function or the result of a call to a family function, see \code{\link{family}} for details; currently only \code{gaussian}, \code{binomial}, \code{poisson} are working}
  \item{method}{fitting method; one out of \code{"lqa"}, \code{"nlm"}, \code{"AIC"} or \code{"BIC"}; methods \code{"lqa"} and \code{"nlm"} induce penalized estimation; the default method \code{"lqa"} employs a \acronym{PIRLS}-algorithm (see Fan and Li, 2001, Ulbricht, 2010); method \code{"nlm"} employs an absolute value approximation similar to this of Ulbricht (2010), but fits the approximated optimization problem directly with function \code{\link{nlm}}. Methods \code{"AIC"} and \code{"BIC"} employ a forward selection strategy like described in Oelker et. al. (2012)}
  \item{tuning}{a list; tuning parameters for penalized estimation; \code{lambda} is the penalty parameter used; \code{phi} weights the penalty's terms, out of intervall )0,1(, the default \code{0.5} corresponds to equal weights; if \code{lambda} and/or \code{phi} are vectors, these values are cross-validated; if \code{lambda} and/or \code{phi} are \code{TRUE}, they will be cross-validated on their whole range}
  \item{weights}{an optional weight vector}
  \item{control}{a list of parameters for controlling the fitting process; if emtpy, set to \code{cat_control()}; see \code{\link{cat_control}}}
  \item{model}{a logical value indicating whether the employed model frame shall be returned or not}
  \item{x, y}{for \code{gvcm.cat}: logical values indicating whether the response vector and model matrix used in the fitting process shall be returned or not; for \code{pest} and \code{abc}: \code{y} must be a response vector}
  \item{X}{only for \code{pest} and \code{abc}: a proper coded design matrix}
  \item{plot}{logical; if TRUE, estimates for path-plotting are computed}
  \item{ind}{for \code{pest} and \code{abc} only: the to be used index arguments; see function \code{\link{index}}}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
A typical \code{\link{formula}} has the form \code{response ~ 1 + terms}; where \code{response} is the response vector and \code{terms} is a series of terms which specifies a linear predictor. 
Varying coefficients enter the \code{\link{formula}} as \code{v(x,u)} where \code{u} denotes the categorial effect modifier and \code{x} the modfied covariate. 
A varying intercept is denoted by \code{v(1,u)} and has stand at the formula's very beginning. If the \code{\link{formula}} contains no (varying) intercept, \code{gvcm.cat} assumes a constant intercept. 
There is no way to avoid an intercept. \cr
Ordinal/nominal covariates \code{u} given as \code{p(u)} are penalized as described in Gertheiss and Tutz (2010). 
For binomial and quasibinomial families the response can also be a success/failure rate or a two-column matrix with the columns giving the numbers of successes and failures. \cr
Function \code{pest} computes \strong{p}enalized \strong{est}imates, that is, it implements methods  \code{"lqa"} (\acronym{PIRLS}-algorithm) and \code{"nlm"}. \cr
Function \code{abc} implements the forward selection strategy employing \strong{A}IC/\strong{B}I\strong{C}. \cr
Categorial effect modifiers and penalized categorial covariates are dummy coded as required by the penalty. If \code{x} in \code{v(x,u)} is binary, it is effect coded (first category refers to -1). Other covariates are coded like given by \code{\link{getOption}}.
}

\value{
\code{gvcm.cat} returns an object of class \dQuote{\code{"gvcm.cat"}} which inherits from class \dQuote{\code{\link{glm}}} which inherits from class \dQuote{\code{\link{lm}}}.
An object of class \dQuote{\code{gvcm.cat}} contains:
  \item{coefficients}{named vector of coefficients}
  \item{coefficients.reduced}{reduced vector of coefficients; selected coefficients/differences of coefficients are set to zero}
  \item{coefficients.refitted}{refitted vector of coefficients; i.e. maximum likelihood estimate of that model containing selected covariates only; same length as \code{coefficients.reduced}}
  \item{coefficients.oml}{maximum likelihood estimate of the full model}
  \item{residuals}{deviance residuals}
  \item{fitted.values}{fitted mean values}
  \item{rank}{degrees of freedom model; for \code{method="lqa"} estimated by the trace of the generalized head matrix; for \code{method="nlm"} the estimate is the number of selected coefficients; for methods \code{"AIC"}, \code{"BIC"} estimated like default in \code{\link{glm.fit}}}
  \item{family}{the \code{\link{family}} object used}
  \item{linear.predictors}{linear fit on link scale}
  \item{deviance}{minus twice the maximized log-likelihood}
  \item{aic}{a version of Akaike's Information Criterion; minus twice the maximized log-likelihood plus twice the number of parameters. For binomial and Poison families the dispersion is fixed at one and the number of parameters is the number of coefficients. For a gaussian family the dispersion is estimated from the residual deviance, and the number of parameters is the number of coefficients plus one.}
  \item{null.deviance}{the deviance for the null model, comparable with \code{deviance}; the null model includes a non-varying intercept only}
  \item{iter}{number of iterations}
  \item{weights}{working weights of the final iteration}
  \item{df.residual}{the residual degrees of freedom/degrees of freedom error; computed like \code{rank}}
  \item{df.null}{the residual degrees of freedom for the null model}
  \item{converged}{logical; fulfills the \acronym{PIRLS}-algorithm the given convergence conditions?}
  \item{boundary}{logical; is the fitted value on the boundary of the attainable values?}
  \item{offset}{the offset vector used}
  \item{control}{the value of the \code{control} argument used}
  \item{method}{same as input argument \code{method}}
  \item{contrasts}{the contrasts used}
  \item{na.action}{information returned by \code{\link{model.frame}} on the special handling of \code{\link{NA}}s; currently always \code{\link{na.omit}}}
  \item{plot}{if input \code{plot=TRUE}, a list containig two matrixes for plotting}
  \item{tuning}{a list, employed values of \code{lambda} and \code{phi}; if \code{lambda} and/or \code{phi} were cross-validated, these are the optimal values}
  \item{index1}{used index1 argument; see function \code{\link{index}}}
  \item{index2}{used index2 argument; see function \code{\link{index}}}
  \item{index3}{used index3 argument; see function \code{\link{index}}}
  \item{number.selectable.parameters}{number of coefficients that could be selected}
  \item{number.removed.parameters}{number of actual removed coefficients}
  \item{x.reduction}{a matrix; transforms model frame \code{x} into its reduced version; e.g. needed for refitting}
  \item{call}{the matched call}
  \item{formula}{the \code{\link{formula}} supplied}
  \item{terms}{the \code{\link{terms}} object used}
  \item{data}{the data argument}
  \item{x, y}{if requested, the model matrix/the response vector}
  \item{model}{if requested, the model frame}
  \item{xlevels}{a record of the levels of the factors used in fitting}
}
\references{
Fan, J. and R. Li (2001). Variable selection via nonconcave penalized likelihood and its oracle properties. \emph{Journal of the American Statistical Association 96}(456), 1348-1360.\cr
Gertheiss, J. and G. Tutz (2010). Sparse modeling of categorial explanatory variables. \emph{The Annals of Statistics 4}(4), 2150-2180.\cr
Oelker, M.-R., J. Gertheiss and G. Tutz (2012). Regularization and Model Selection with Categorial Predictors and Effect Modifiers in Generalized Linear Models. \emph{Department of Statistics at the University of Munich: Technical Report 122}.\cr
Ulbricht, J. (2010). \emph{Variable Selection in Generalized Linear Models}. Dissertation, Department
of Statistics, University of Munich: Verlag Dr. Hut.\cr
}
\author{Margret-Ruth Oelker (\email{margret.oelker@stat.uni-muenchen.de})}
\note{There is a summary function: \code{summary.gvcm.cat}}
\seealso{
Functions \code{\link{index}}, \code{\link{cat_control}}, \code{\link{plot.gvcm.cat}}, \code{\link{predict.gvcm.cat}}, \code{\link{simulation}}
}
\examples{
## continues example of function simulation 
\donttest{f <- y ~ v(1,u) + v(x1,u) + v(x2,u)
m1 <- gvcm.cat(f, data, binomial(), plot=TRUE)
summary(m1)}
}
\keyword{gvcm.cat}
