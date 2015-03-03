\name{gvcm.cat}
\alias{gvcm.cat}
\alias{pest}
\alias{abc}
\alias{print.gvcm.cat}
\alias{summary.gvcm.cat}
\title{Regularized Categorical Effects/Categorical Effect Modifiers/Continuous/Smooth effects in GLMs}
\description{
The function fits generalized linear models with regularized categorical effects, categorical effect modifiers, continuous effects and smooth effects. 
The model is specified by giving a symbolic description of the linear predictor and a description of the error distribution. 
Estimation employs different regularization and model selection strategies. 
These strategies are either a penalty or a forward selection strategy employing AIC/BIC. 
For non-differentiable penalties, a local quadratic approximation is employed, see Oelker and Tutz (2013). 
}
\usage{
gvcm.cat(formula, data, family = gaussian, method = c("lqa", "AIC", "BIC"), 
tuning = list(lambda=TRUE, specific=FALSE, phi=0.5, grouped.fused=0.5, 
elastic=0.5, vs=0.5, spl=0.5), weights, offset, start, control, 
model = FALSE, x = FALSE, y = FALSE, plot=FALSE, ...)

pest(x, y, indices, family = gaussian, 
tuning = list(lambda=TRUE, specific=FALSE, phi=0.5, grouped.fused=0.5, 
elastic=0.5, vs=0.5, spl=0.5), weights, offset, start = NULL, 
control = cat_control(), plot=FALSE, ...)

abc(x, y, indices, family = gaussian, tuning = c("AIC", "BIC"), 
weights, offset, start, control = cat_control(), plot=FALSE, ...)
}
\arguments{
  \item{formula}{an object of class \code{\link{formula}}: a symbolic description of the model to be fitted. See details}
  \item{data}{a data frame, containing the variables in the model}
  \item{family}{a \code{family} object describing the error distribution and link function to be used in the model; 
       this can be a character string naming a family function, a family function or the result of a call to a family function, 
       see \code{\link{family}} for details; currently only \code{gaussian}, \code{binomial}, \code{poisson}, \code{Gamma} 
       are working}
  \item{method}{fitting method; one out of \code{"lqa"}, \code{"AIC"} or \code{"BIC"}; method \code{"lqa"} induces penalized estimation; 
       it employs a \acronym{PIRLS}-algorithm (see Fan and Li, 2001; Oelker and Tutz, 2013). 
       Methods \code{"AIC"} and \code{"BIC"} employ a forward selection strategy}
  \item{tuning}{a list; tuning parameters for penalized estimation; \code{lambda} is the scalar, overall penalty parameter; 
       if \code{lambda} is a vector of values, these values are cross-validated; 
       if \code{lambda = TRUE}, \code{lambda} is cross-validated on log scale between \code{lambda.lower} and \code{lambda.upper}; see \code{\link{cat_control}}.
       If \code{lambda} is a vector with the same length as elements in the \code{formula} and if \code{specific} equals a vector of proper length, the entries of \code{specific} are interpreted as specific tuning parameters for each entry of the formula. 
       \code{phi}, \code{grouped.fused}, \code{elastic}, \code{vs} and \code{spl} are parameters that weigh the terms of some penalties; 
       must be out of intervall )0,1(; the default \code{0.5} corresponds to equal weights
       }
  \item{weights}{an optional weight vector (for the observations)}
  \item{offset}{an optional offset}
  \item{start}{initial values for the PIRLS algorithm for method \code{lqa}}
  \item{control}{a list of parameters for controlling the fitting process; if emtpy, set to \code{cat_control()}; see \code{\link{cat_control}}}
  \item{model}{for functions \code{gvcm.cat}: a logical value indicating whether the employed model frame shall be returned or not}
  \item{x, y}{for function \code{gvcm.cat}: logical values indicating whether the response vector and the model matrix used in the fitting process shall be returned or not; 
       for functions \code{pest} and \code{abc}: \code{y} must be a response vector, \code{x} a proper coded design matrix}
  \item{plot}{logical; if \code{TRUE}, estimates needed to plot coefficient paths are computed}
  \item{indices}{for \code{pest} and \code{abc} only: the to be used index argument; see function \code{\link{index}}}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
% formula zeug: 
A typical \code{\link{formula}} has the form \code{response ~ 1 + terms}; where \code{response} is the response vector and \code{terms} is a series of terms which specifies a linear predictor. 
There are some special terms for regularized terms: 
  \itemize{
    \item  \code{v(x, u, n="L1", bj=TRUE)} : varying coefficients enter the \code{\link{formula}} as \code{v(x,u)} where \code{u} denotes the categorical effect modifier and \code{x} the modfied covariate. 
           A varying intercept is denoted by \code{v(1,u)}. Varying coefficients with categorical effect modifiers are penalized as described in Oelker et. al. 2012. 
           The argument \code{bj} and the element \code{phi} in argument \code{tuning} allow for the described weights. 
    \item  \code{p(u, n="L1")} : ordinal/nominal covariates \code{u} given as \code{p(u)} are penalized as described in Gertheiss and Tutz (2010). For numeric covariates, \code{p(u)} indicates a pure Lasso penalty.
    \item  \code{grouped(u, ...)}: penalizes a group of covariates with the grouped Lasso penalty of Yuan and Lin (2006); so far, working for categorical covariates only
    \item  \code{sp(x, knots=20, n="L2")}: implents a continuous \code{x} covariate non-parametrically as \eqn{f(x)}; \eqn{f(x)} is represented by centered evaluations of basis functions (cubic B-splines with number of knots = \code{knots}); for \code{n="L2"}, the curvature of \eqn{f(x)} is penalized by a Ridge penalty; see Eilers and Marx (1996)
    \item  \code{SCAD(u)}: penalizes a covariate \code{u} with the SCAD penalty by Fan and Li (2001); for categorical covariates \code{u}, differences of coefficients are penalized by a SCAD penalty, see Gertheiss and Tutz (2010)
    \item  \code{elastic(u)}: penalizes a covariate \code{u} with the elastic net penalty by Zou and Hastie (2005); for categorical covariates \code{u}, differences of coefficients are penalized by the elastic net penalty, see Gertheiss and Tutz (2010)
  }
If the \code{\link{formula}} contains no (varying) intercept, \code{gvcm.cat} assumes a constant intercept. There is no way to avoid an intercept.\cr
For specials \code{p} and \code{v}, there is the special argument \code{n}: 
if \code{n="L1"}, the absolute values in the penalty are replaced by squares of the same terms; 
if \code{n="L2"}, the absolute values in the penalty are replaced by quadratic, Ridge-type terms; 
if \code{n="L0"}, the absolute values in the penalty are replaced by an indicator for non-zero entries of the same terms.\cr
For methods \code{"AIC"} and \code{"BIC"}, the coefficients are not penalized but selected by a forward selection strategy whenever it makes sense; 
for special \code{v(x,u)}, the selection strategy is described in Oelker et. al. 2012; the approach for the other specials corresponds to this idea. \cr
For binomial families the response can also be a success/failure rate or a two-column matrix with the columns giving the numbers of successes and failures. \cr
Function \code{pest} computes \strong{p}enalized \strong{est}imates, that is, it implements method  \code{"lqa"} (\acronym{PIRLS}-algorithm). \cr
Function \code{abc} implements the forward selection strategy employing \strong{A}IC/\strong{B}I\strong{C}. \cr
Categorical effect modifiers and penalized categorical covariates are dummy coded as required by the penalty. If \code{x} in \code{v(x,u)} is binary, it is effect coded (first category refers to -1). Other covariates are coded like given by \code{\link{getOption}}.\cr
There is a summary function: \code{summary.gvcm.cat}
}

\value{
\code{gvcm.cat} returns an object of class \dQuote{\code{gvcm.cat}} which inherits from class \dQuote{\code{\link{glm}}} which inherits from class \dQuote{\code{\link{lm}}}.
An object of class \dQuote{\code{gvcm.cat}} contains:
  \item{coefficients}{named vector of coefficients}
  \item{coefficients.reduced}{reduced vector of coefficients; selected coefficients/differences of coefficients are set to zero}
  \item{coefficients.refitted}{refitted vector of coefficients; i.e. maximum likelihood estimate of that model containing selected covariates only; same length as \code{coefficients.reduced}}
  \item{coefficients.oml}{maximum likelihood estimate of the full model}
  \item{residuals}{deviance residuals}
  \item{fitted.values}{fitted mean values}
  \item{rank}{degrees of freedom model; for \code{method="lqa"} estimated by the trace of the generalized head matrix; for methods \code{"AIC"}, \code{"BIC"} estimated like default in \code{\link{glm.fit}}}
  \item{family}{the \code{\link{family}} object used}
  \item{linear.predictors}{linear fit on link scale}
  \item{deviance}{scaled deviance}
  \item{aic}{a version of Akaike's Information Criterion; minus twice the maximized log-likelihood plus twice the rank. For binomial and Poison families the dispersion is fixed at one. For a gaussian family the dispersion is estimated from the residual deviance, and the number of parameters is the rank plus one.}
  \item{null.deviance}{the deviance for the null model, comparable with \code{deviance}; the null model includes a non-varying intercept only}
  \item{iter}{number of iterations}
  \item{weights}{working weights of the final iteration}
  \item{df.residual}{the residual degrees of freedom/degrees of freedom error; computed like \code{rank}}
  \item{df.null}{the residual degrees of freedom for the null model}
  \item{converged}{logical; fulfills the \acronym{PIRLS}-algorithm the given convergence conditions?}
  \item{boundary}{logical; is the fitted value on the boundary of the attainable values?}
  \item{offset}{the offset vector used}
  \item{control}{the value of the \code{control} argument used}
  \item{contrasts}{the contrasts used}
  \item{na.action}{information returned by \code{\link{model.frame}} on the special handling of \code{\link{NA}}s; currently always \code{\link{na.omit}}}
  \item{plot}{in principle, a list containing two matrixes needed for different types of plots: 
       if input option \code{plot=TRUE}, the first matrix contains estimates needed to plot coefficient paths;
       if \code{lambda} was cross-validated, the second matrix contains the cross-validation scores}
  \item{tuning}{a list, employed tuning parameters; if \code{lambda} was cross-validated, the optimal value is returned}
  \item{indices}{used index argument; see function \code{\link{index}}}
  \item{number.selectable.parameters}{number of coefficients that could be selected}
  \item{number.removed.parameters}{number of actual removed coefficients}
  \item{x.reduction}{a matrix; transforms model frame \code{x} into its reduced version; e.g. needed for refitting}
  \item{beta.reduction}{a matrix; transforms the \code{coefficients} into its reduced version}
  \item{call}{the matched call}
  \item{formula}{the \code{\link{formula}} supplied}
  \item{terms}{the \code{\link{terms}} object used}
  \item{data}{the data argument}
  \item{x, y}{if requested, the model matrix/the response vector}
  \item{model}{if requested, the model frame}
  \item{xlevels}{a record of the levels of the factors used in fitting}
  \item{bootstrap.errors}{experimental}
  \item{method}{same as input argument \code{method}}
In addition, non-empty fits will have components \code{qr}, \code{R} and \code{effects} relating to the final weighted linear fit. 
}

\references{
Eilers, P. H. C. and B. D. Marx (1996). Flexible smoothing with b-splines and penalties. Statist. Sci. 11 (2), 89-121.\cr
Fan, J. and R. Li (2001). Variable selection via nonconcave penalized likelihood and its oracle properties. 
\emph{Journal of the American Statistical Association 96}(456), 1348-1360.\cr
Gertheiss, J. and G. Tutz (2010). Sparse modeling of categorial explanatory variables. 
\emph{The Annals of Statistics 4}(4), 2150-2180.\cr
Oelker, M.-R., J. Gertheiss and G. Tutz (2012). Regularization and model melection with categorial predictors and effect 
modifiers in generalized linear models. \emph{Department of Statistics at the University of Munich: Technical Report 122}.\cr
Oelker, M.-R., J. Gertheiss and G. Tutz (2013). A general family of penalties for combining differing types of penalties in 
generalized structured models. \emph{Department of Statistics at the University of Munich: Technical Report 139}.\cr
Yuan, M. and Y. Lin (2006). Model selection and estimation in regression with grouped variables. R. Stat. Soc. Ser. B Stat. 
Methodol. 68 (1), 49-67.\cr
Zou, H. and T. Hastie (2005). Regularization and variable selection via the Elastic Net. R. Stat. Soc. Ser. B Stat. 
Methodol. 67 (2), 301-320.\cr
}
\note{Pleas note that the functions \code{gvcm.cat}, \code{pest} and the fitting procedure 
for penalized estimation \code{gvcmcatfit} are organized like the functions \code{\link{glm}}/\code{\link{glm.fit}} whenever possible. 
This was done to avoid mistakes and to provide a well-known structure.
}
\author{Margret-Ruth Oelker (\email{margret.oelker@stat.uni-muenchen.de})}
\seealso{
Functions \code{\link{index}}, \code{\link{cat_control}}, \code{\link{plot.gvcm.cat}}, \code{\link{predict.gvcm.cat}}, \code{\link{simulation}}
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
## example for function gvcm.cat()
f <- y ~ v(1,u) + v(x1,u) + v(x2,u)
m1 <- gvcm.cat(f, data, binomial(), plot=TRUE, control=cat_control(lambda.upper=19))
summary(m1)
## example for function predict.gvcm.cat
newdata <- simulation(200, covariates, NULL, true.f, true.coefs , binomial(), seed=789)
prediction <- predict.gvcm.cat(m1, newdata) 
## example for function plot.gvcm.cat 
plot(m1)
plot(m1, type="score")
plot(m1, type="coefs")
}

