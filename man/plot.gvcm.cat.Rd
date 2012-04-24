\name{plot.gvcm.cat}
\alias{plot.gvcm.cat}
\title{Plot Method for gvcm.cat Objects}
\description{
\code{gvcm.cat()} offers two ways of visualizing a \code{\link{gvcm.cat}} object with penalized estimation. 
Default option \code{type="path"} delivers a graphic with the coefficient paths between 0 (= maximal penalization) and 1 (= no penalization). Maximal penalization is defined by the minimal penalty parameter \code{lambda} that sets all penalized coefficients to zero (to constant relating to the intercept and \code{assured.intercept = TRUE}). Minimal penalization means no penalization at all, i.e. \code{lambda = 0}. Of course the minimal penalty parameter causing maximal penalization depends on how selection and clustering of coefficients is defined (see function \code{\link{gvcm.cat}} and \code{\link{cat_control}}). Coefficients belonging to one covariate are plotted in the same color, coefficients that are not modified are plotted as dashed lines.\cr
Option \code{type="score"} plots the cross-validation score (depending on \code{criterion} in \code{\link{cat_control}}) as a function of penalty parameter \code{lambda} and marks the chosen penalty parameter as a dotted line.
}
\usage{
\method{plot}{gvcm.cat}(x, accuracy = 2, type = "path", individual.paths = FALSE, 
xlim, ylim, main = NULL, indent = 0, color = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{\link{gvcm.cat}} object with value \code{plot} unequal \code{NA}}
  \item{accuracy}{integer; number of digits being compared when setting coefficents equal/to zero for plotting}
  \item{type}{one out of \code{"path"}, \code{"score"}; defines the type of plot}
  \item{individual.paths}{logical; indicates whether the paths of all coefficients shall be plotted into one common figure (default) or in an individual figure per covariate; paths of single covariates can be selected by giving a vector containing the covariates' names as characters}
  \item{xlim}{the \code{x} limits \code{(x1, x2)} of the plot}
  \item{ylim}{the \code{y} limits \code{(y1, y2)} of the plot}
  \item{main}{title of the plot}
  \item{indent}{numeric; if larger zero, coefficient names printed on top of each other are adjusted}
  \item{color}{logical; if \code{FALSE}, paths are gray and dotted/dashed}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{Paths are drawn by connecting \code{steps} estimates related to different values of \code{lambda} and constant \code{phi}, see \code{\link{cat_control}}. There is no plot for methods \code{"AIC"} and \code{"BIC"}.}
\value{A plot.}

\seealso{Function \code{\link{gvcm.cat}}}
\examples{
## continues example of function gvcm.cat 
\donttest{plot(m1)}
}
\keyword{gvcm.cat}
