\name{cv.lla.gcdnet}
\alias{cv.lla.gcdnet}
\alias{cv.hsvmpath}
\alias{cv.sqsvmpath}
\alias{cv.logitpath}
\alias{cv.lspath}
\alias{cv.erpath}
\alias{cv.probitpath}
\title{Cross-validation for lla.gcdnet}
\description{Does k-fold cross-validation for lla.gcdnet, produces a plot,
and returns a value for \code{lambda}. This function is modified based on the \code{cv} function from the \code{glmnet} package.}
\usage{
cv.lla.gcdnet(x, y, a, lambda, pred.loss, nfolds, foldid, delta, omega,...)
}
\arguments{
		\item{x}{\code{x} matrix as in \code{\link{lla.gcdnet}}.}
		\item{y}{response variable or class label \code{y} as in \code{\link{lla.gcdnet}}.}
		\item{a}{The constant for the folded concave penalty (SCAD) - default is 3.7.}
		\item{lambda}{optional user-supplied lambda sequence; default is
		\code{NULL}, and \code{\link{lla.gcdnet}} chooses its own sequence.}
		\item{nfolds}{number of folds - default is 5. Although \code{nfolds}
		can be as large as the sample size (leave-one-out CV), it is not
		recommended for large datasets. Smallest value allowable is \code{nfolds=3}.}
		\item{foldid}{an optional vector of values between 1 and \code{nfold}
		identifying what fold each observation is in. If supplied,
		\code{nfold} can be missing.}
		\item{pred.loss}{loss function to use for cross-validation error. Valid options are:
		\itemize{
		\item \code{"loss"} Margin based loss function. When use least square loss \code{"ls"}, it gives mean square error (MSE). When use expectile regression loss \code{"er"}, it gives asymmetric mean square error (AMSE).
		\item \code{"misclass"} only available for classification: it gives misclassification error. }
		Default is \code{"loss"}.}
		\item{delta}{parameter \eqn{\delta}{delta} only used in HHSVM for computing margin based loss function, only available for \code{pred.loss = "loss"}.}
		\item{omega}{parameter \eqn{\omega}{omega} only used in expectile regression. Only available for \code{pred.loss = "loss"}.}
		\item{\dots}{other arguments that can be passed to \code{lla.gcdnet}.}

}
\details{The function runs \code{\link{lla.gcdnet}} \code{nfolds}+1 times; the
first to get the \code{lambda} sequence, and then the remainder to
compute the fit with each of the folds omitted. The average error and standard deviation over the folds are computed. 
}
\value{an object of class \code{\link{cv.lla.gcdnet}} is returned, which is a
list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{a}{the constant for the folded concave penalty (SCAD).}
		\item{cvm}{the mean cross-validated error - a vector of length
		\code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{cvupper}{upper curve = \code{cvm+cvsd}.}
		\item{cvlower}{lower curve = \code{cvm-cvsd}.}
		\item{nzero}{number of non-zero coefficients at each \code{lambda}.}
		\item{name}{a text string indicating type of measure (for plotting
		purposes).}
		\item{lla.gcdnet.fit}{a fitted \code{\link{lla.gcdnet}} object for the full data.}
		\item{lambda.min}{The optimal value of \code{lambda} that gives minimum
		cross validation error \code{cvm}.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is
		within 1 standard error of the minimum.}
}

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


\seealso{
\code{\link{lla.gcdnet}}, 
\code{\link{plot.cv.lla.gcdnet}}, 
\code{\link{predict.cv.lla.gcdnet}}, 
and \code{\link{coef.cv.lla.gcdnet}} methods.}
\examples{
# fit a folded concave (SCAD) penalized HHSVM with 
# nlambda=20. Use the misclassification rate as the 
# cross validation prediction loss. Use five-fold CV 
# to choose the optimal lambda for the regularization 
# parameter.

data(FHT)
set.seed(2021)
cv = cv.lla.gcdnet(FHT$x, FHT$y, method ="hhsvm",
             nlambda=20, pred.loss="misclass",
             nfolds=5, delta=1.5)
plot(cv)


# fit a folded concave (SCAD) penalized least squares 
# with nlambda = 20. Use the least square loss as the 
# cross validation prediction loss. Use five-fold CV to 
# choose the optimal regularization parameter lamdba.

set.seed(2021)
cv1 = cv.lla.gcdnet(FHT$x, FHT$y_reg, method ="ls", 
                    nlambda=20, pred.loss="loss", nfolds=5)
plot(cv1)


# To fit a folded concave (SCAD) penalized logistic 
# regression with nlambda = 20. Use the logistic loss as 
# the cross validation prediction loss. Use five-fold CV 
# to choose the optimal lambda for the regularization 
# parameter.

set.seed(2021)
cv2 = cv.lla.gcdnet(FHT$x, FHT$y, method ="logit",
                    nlambda=20, pred.loss="loss", nfolds=5)
plot(cv2)


# fit a a folded concave (SCAD) penalized probit
# regression with nlambda = 20. Use the probit loss as 
# the cross validation prediction loss. Use five-fold CV 
# to choose the optimal lambda for the regularization 
# parameter.

set.seed(2021)
cv3 = cv.lla.gcdnet(FHT$x, FHT$y, method ="probit", nlambda=20, 
                    pred.loss="loss", nfolds=5)
plot(cv3)
}
\keyword{models}
\keyword{regression}

