\name{plot.gvcm.cat}
\alias{plot.gvcm.cat}
\title{Plot Method for gvcm.cat Objects}
\description{
Function to visualize a \code{\link{gvcm.cat}} object. 
}
\usage{
\method{plot}{gvcm.cat}(x, accuracy = 2, type = "path", individual = FALSE, 
xlim, ylim, main = NULL, indent = 0, color = TRUE, xscale = "lambda", 
label = TRUE, intercept = TRUE, ...)
}
\arguments{
  \item{x}{a \code{\link{gvcm.cat}} object; for \code{type="path"}, a \code{\link{gvcm.cat}} object with value \code{plot} unequal \code{NA} is required}
  \item{accuracy}{integer; number of digits being compared when setting coefficents equal/to zero for plotting}
  \item{type}{one out of \code{"path"}, \code{"score"}, \code{"coefs"}; defines the type of the plot}
  \item{individual}{logical; for \code{type="path"} and \code{type="coefs"} only;
         for \code{type="path"}, it indicates whether the paths of all coefficients shall be plotted 
         into one common figure (default) or in an individual figure per covariate; paths of single 
         covariates can be selected by giving a vector containing the covariates (as characters and as 
         given in the formula, e.g.: \code{individual.paths=c("v(1,u)", "v(x1,u1)")})
         for \code{type="coefs"}, the default is one plot per covariate. \code{individual} allows to select single covariates.
         }
  \item{xlim}{the \code{x} limits \code{(x1, x2)} of the plot}
  \item{ylim}{the \code{y} limits \code{(y1, y2)} of the plot}
  \item{main}{title of the plot}
  \item{indent}{numeric; if larger zero, coefficient names printed on top of each other are adjusted}
  \item{color}{logical; if \code{FALSE}, lines are gray and dotted/dashed}
  \item{xscale}{for \code{type="path"} only; if \code{xscale="lambda"}, the x-axis is scaled as \eqn{1 - \lambda/\lambda_{max}}; if \code{xscale="beta"}, the scale of the x-axis is the scaled L1 norm of the penalized coefficients.}  
  \item{label}{omits addtional information printed in the plot, if \code{FALSE}} 
  \item{intercept}{for \code{type="coefs"} and \code{type="path"} only; if \code{FALSE}, for \code{type="path"}, the path of the intercept is not plotted; if \code{FALSE}, for \code{type="coefs"}, intercept is not added to smooth functions} 
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
Default option \code{type="path"} delivers a graphic with the coefficient paths between 0 (= maximal penalization) and 1 (= no penalization). 
Maximal penalization is defined by the minimal penalty parameter \code{lambda} that sets all penalized coefficients to zero (to constant relating to the intercept and \code{assured.intercept = TRUE}). 
Minimal penalization means no penalization at all, i.e. \code{lambda = 0}. 
Of course the minimal penalty parameter causing maximal penalization depends on how selection and clustering of coefficients is defined (see function \code{\link{gvcm.cat}} and \code{\link{cat_control}}). 
Coefficients belonging to one covariate are plotted in the same color, coefficients that are not modified are plotted as dashed lines.
Paths are drawn by connecting \code{steps} estimates related to different values of \code{lambda}, see \code{\link{cat_control}}.\cr
Option \code{type="score"} plots the cross-validation score (depending on \code{criterion} in \code{\link{cat_control}}) as a function of penalty parameter \code{lambda} and marks the chosen penalty parameter as a dotted line.\cr
Opton \code{type="coefs"} plots the penalized coefficients whenever possible. \cr
So far, there is no plot for methods \code{"AIC"} and \code{"BIC"}.}
\value{A plot.}

\seealso{Function \code{\link{gvcm.cat}}}
\examples{
## see example for function gvcm.cat 
}
\keyword{gvcm.cat}
