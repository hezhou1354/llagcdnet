\name{print.lla.gcdnet}
\alias{print.lla.gcdnet}
\title{print an "lla.gcdnet" object}
\description{
Print a summary of the \code{\link{lla.gcdnet}} path at each step along the path.
}
\usage{
\method{print}{lla.gcdnet}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{fitted \code{\link{gcdnet}} object}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
The call that produced the \code{\link{lla.gcdnet}} object is printed, followed by a two-column matrix with columns \code{Df} and \code{Lambda}. The \code{Df} column is the number of nonzero coefficients.
}
\value{
a two-column matrix, the first columns is the number of nonzero coefficients and the second column is \code{Lambda}.}

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

\examples{
data(FHT)
m1 <- lla.gcdnet(x=FHT$x,y=FHT$y,nlambda=20)
print(m1)
}
\keyword{models}
\keyword{regression}

