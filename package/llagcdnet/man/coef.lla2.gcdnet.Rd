\name{coef.lla2.gcdnet}
\alias{coef.lla2.gcdnet}
\alias{coef.hsvmpath}
\alias{coef.sqsvmpath}
\alias{coef.logitpath}
\alias{coef.lspath}
\alias{coef.erpath}
\alias{coef.probitpath}
\title{get coefficients or make coefficient predictions from an "lla2.gcdnet" object.}
\description{
Computes the coefficients or returns a list of the indices of the nonzero coefficients at the requested values for \code{lambda} from a fitted \code{\link{lla2.gcdnet}} object.
}
\usage{
\method{coef}{lla2.gcdnet}(object, s = NULL,
type=c("coefficients","nonzero"), ...)
}
\arguments{
		\item{object}{fitted \code{\link{lla.gcdnet}} model object.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the entire sequence used to
		create the model.}
		\item{type}{type \code{"coefficients"} computes the coefficients at the requested
		values for \code{s}. Type
		\code{"nonzero"} returns a list of the indices of the nonzero
		coefficients for each value of \code{s}. Default is \code{"coefficients"}.}
		\item{\dots}{not used. Other arguments to predict. } 
}

\details{
\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{coef} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of coefficients from both left and right \code{lambda} indices. 
}
\value{The object returned depends on type.}


\author{He Zhou and Hui Zou\cr
Maintainer: He Zhou  <zhou1354@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
BugReport: \url{https://github.com/emeryyi/fastcox.git}\cr


Fan, J., Xue, L., & Zou, H. (2014), "Strong oracle optimality of folded concave penalized estimation", \emph{Annals of statistics}, 42(3), 819.\cr
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4295817/}\cr


Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{predict.gcdnet}} method}
\examples{
data(FHT)
fit1 = lla2.gcdnet(x=FHT$x,y=FHT$y, nlambda=20)
coef(fit1,type="coef",s=c(0.1,0.005))
coef(fit1,type="nonzero")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
