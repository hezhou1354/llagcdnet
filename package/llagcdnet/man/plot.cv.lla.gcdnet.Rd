\name{plot.cv.lla.gcdnet}
\alias{plot.cv.lla.gcdnet}
\title{plot the cross-validation curve produced by cv.lla.gcdnet}
\description{
Plots the cross-validation curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv} function from the \code{glmnet} package.}
\usage{
\method{plot}{cv.lla.gcdnet}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{fitted \code{\link{cv.lla.gcdnet}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or
	its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}
\details{A plot is produced.}

\author{He Zhou and Hui Zou\cr
Maintainer: He Zhou <zhou1354@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
BugReport: \url{https://github.com/emeryyi/fastcox.git}\cr


Fan, J., Xue, L., & Zou, H. (2014), "Strong oracle optimality of folded concave penalized estimation", \emph{Annals of statistics}, 42(3), 819.\cr
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4295817/}\cr


Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{cv.lla.gcdnet}}.}
\examples{
# To fit a folded concave (SCAD) penalized logistic 
# regression with nlambda = 20. Use the logistic loss as 
# the cross validation prediction loss. Use five-fold CV 
# to choose the optimal lambda for the regularization 
# parameter.
data(FHT)
set.seed(2021)
cv1 = cv.lla.gcdnet(FHT$x, FHT$y, method ="logit",
                    nlambda=20, pred.loss="loss", nfolds=5)
plot(cv1)


# fit a a folded concave (SCAD) penalized probit
# regression with nlambda = 20. Use the probit loss as 
# the cross validation prediction loss. Use five-fold CV 
# to choose the optimal lambda for the regularization 
# parameter.
data(FHT)
set.seed(2021)
cv2 = cv.lla.gcdnet(FHT$x, FHT$y, method ="probit",
                    nlambda=20, pred.loss="loss", nfolds=5)
plot(cv2)
}
\keyword{models}
\keyword{regression}
