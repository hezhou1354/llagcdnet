\name{predict.cv.lla2.gcdnet}
\alias{predict.cv.lla.gcdnet}
\title{make predictions from a "cv.lla2.gcdnet" object.}
\description{
This function makes predictions from a cross-validated gcdnet model,
using the stored \code{"lla2.gcdnet.fit"} object, and the optimal value
chosen for \code{lambda}.
}
\usage{
\method{predict}{cv.lla2.gcdnet}(object, newx, s=c("lambda.1se","lambda.min"),...)
}
\arguments{
		\item{object}{fitted \code{\link{cv.lla2.gcdnet}} object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. Must be a matrix. See documentation for \code{\link{predict.lla2.gcdnet}}.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the value \code{s="lambda.1se"} stored
		on the CV object. Alternatively \code{s="lambda.min"} can be
		used. If \code{s} is numeric, it is taken as the value(s) of
		\code{lambda} to be used.}
		\item{\dots}{not used. Other arguments to predict. } }
\details{This function makes it easier to use the results of
  cross-validation to make a prediction.}
\value{The object returned depends the \dots argument which is passed on
to the \code{\link{predict}} method for \code{\link{lla2.gcdnet}} objects.}

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


\seealso{\code{\link{cv.lla2.gcdnet}}, and \code{\link{coef.cv.lla2.gcdnet}} methods.}
\examples{
data(FHT)
set.seed(2021)
cv = cv.lla2.gcdnet(FHT$x, FHT$y, method="probit", nlambda=20, 
                   pred.loss="misclass", nfolds=5)
pre = predict(cv, newx = FHT$x, 
              s = "lambda.1se", 
              type = "class")
}
\keyword{models}
\keyword{regression}
