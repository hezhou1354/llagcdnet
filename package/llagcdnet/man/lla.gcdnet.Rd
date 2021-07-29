\name{lla.gcdnet}
\alias{lla.gcdnet}
\title{Fits the folded concave regularization (SCAD) paths for large margin classifiers}
\description{Fits a regularization path for the large margin classifiers. The regularization path is computed for the folded concave penalty (SCAD) at a grid of values for the regularization parameter lambda. The three-step local linear approximation (LLA) algorithm initialized by zero is applied.}
\usage{
lla.gcdnet(x, y, lambda = NULL, a = 3.7, delta = 2, omega = 0.5,
           thresholds = NULL, ...)
}
\arguments{
		\item{x}{\code{x} matrix as in \code{\link{gcdnet}}.}

		\item{y}{response variable or class label \code{y} as in \code{\link{gcdnet}}.}
		\item{lambda}{optional user-supplied lambda sequence; default is \code{NULL}, and \code{\link{gcdnet}} chooses its own sequence.}
		\item{a}{The constant for the folded concave penalty (SCAD) - default is 3.7.}
		\item{delta}{parameter \eqn{\delta}{delta} only used in HHSVM for computing margin based loss function. The value must be greater than 0. Default is 2.}
		\item{omega}{parameter \eqn{\omega}{omega} only used in expectile regression. The value must be in (0,1). Default is 0.5.}
		\item{thresholds}{parameter \eqn{thresholds} only used in composite probit regression; default is \code{NULL}, and \code{\link{gcdnet}} chooses its own sequence.}
		\item{\dots}{other arguments that can be passed to \code{\link{gcdnet}}.}


		
}

\details{
Note that the objective function in \code{lla.gcdnet} is
\deqn{Loss(y, X, \beta))/N + P_\lambda(|\beta|)}
where the penalty \eqn{P_\lambda(|\beta|)} is the folded concave penalty whose derivative is
\deqn{P'_\lambda(t) = \lambda*I\{t \le \lambda\} + [\max(a*\lambda - t, 0)/(a - 1)]*I\{t>\lambda\} }
for some \eqn{a > 2}, which is specified by the users. Users can also specify the loss function to use, options include Huberized squared hinge loss, Squared hinge loss, least square loss, logistic regression,expectile regression loss and probit regression loss. 

For computing speed reason, if models are not converging or running slow, consider increasing \code{eps} or decreasing \code{nlambda} before increasing \code{maxit}.


}


\value{
An object with S3 class \code{\link{lla.gcdnet}}.
		\item{call}{the call that produced this object}
		\item{a}{the constant for the folded concave penalty (SCAD).}
		\item{b0}{intercept sequence of length \code{length(lambda)}}
		\item{beta}{a \code{p*length(lambda)} matrix of coefficients, stored as a sparse matrix (\code{dgCMatrix} class, the standard class for sparse numeric matrices in the \code{Matrix} package.). To convert it into normal type matrix use \code{as.matrix()}.}
		\item{lambda}{the actual sequence of \code{lambda} values used}
		\item{df}{the number of nonzero coefficients for each value of
		\code{lambda}.}
		\item{dim}{dimension of coefficient matrix (ices)}
		\item{nthresholds}{the actual number of \code{nthresholds} used for composite probit regression, only available when \code{method="composite.probit"}.}
		\item{thresholds}{the actual sequence of \code{thresholds} values used for composite probit regression, only available when \code{method="composite.probit"}.}
		\item{weights}{the actual weights used for composite probit regression, only available when \code{method="composite.probit"}.}
}

\author{He Zhou and Hui Zou\cr
Maintainer: He Zhou  <zhou1354@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
BugReport: \url{https://github.com/emeryyi/fastcox.git}\cr

Fan, J., Xue, L., & Zou, H. (2014), "Strong oracle optimality of folded concave penalized estimation", \emph{Annals of statistics}, 42(3), 819.\cr
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4295817/}\cr


}

\seealso{\code{\link{plot.lla.gcdnet}}.}

\examples{
data(FHT)
# 1. solution paths for the SCAD penalized least squares.

m1 <- lla.gcdnet(x=FHT$x, y=FHT$y_reg, nlambda=20, method="ls")
plot(m1)


# 2. solution paths for the SCAD penalized HHSVM.

m2 <- lla.gcdnet(x=FHT$x, y=FHT$y, nlambda=20, method="hhsvm")
plot(m2)


# 3. solution paths for the SCAD penalized SVM 
# with the squared hinge loss. 

m3 <- lla.gcdnet(x=FHT$x, y=FHT$y, nlambda=20, method="sqsvm")
plot(m3)


# 4. solution paths for the SCAD penalized logistic regression.

m4 <- lla.gcdnet(x=FHT$x, y=FHT$y, nlambda=20, method="logit")
plot(m4)


# 5. solution paths for the SCAD penalized expectile regression
# with the asymmetric least square parameter omega=0.9.

m5 <- lla.gcdnet(x=FHT$x, y=FHT$y_reg, omega=0.9, nlambda=20, method="er")
plot(m5)


# 6. solution paths for the SCAD penalized probit regression.

m6 <- lla.gcdnet(x=FHT$x, y=FHT$y, nlambda=20, method="probit")
plot(m6)


# 7. solution paths for the SCAD penalized composite-probit regression.

m7 <- lla.gcdnet(x=FHT$x, y=FHT$y_reg, nlambda=20, method="composite.probit")
plot(m7)

}
\keyword{models}
\keyword{regression}
