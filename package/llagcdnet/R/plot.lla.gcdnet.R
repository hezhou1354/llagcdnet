plot.lla.gcdnet <- function(x, xvar = c("norm", "lambda"), 
                        color = FALSE, label = FALSE, ...) {
  beta <- x$beta
  lambda <- x$lambda
  df <- x$df
  xvar <- match.arg(xvar)
  ##beta should be in 'dgCMatrix' format
  which <- nonzero(beta)
  beta <- as.matrix(beta[which, ])
  xvar <- match.arg(xvar)
  switch(xvar, norm = {
    index <- apply(abs(beta), 2, sum)
    iname <- "L1 Norm"
  }, lambda = {
    index <- log(lambda)
    iname <- "Log Lambda"
  })
  xlab <- iname
  ylab <- "Coefficients"
  dotlist <- list(...)
  type <- dotlist$type
  if (is.null(type)) {
    if (color == FALSE) 
      matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, 
              type = "l", pch = 500, col = rainbow(12, 
                                                   start = 0.7, end = 0.95), ...) else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, 
                                                                                               type = "l", pch = 500, ...)
  } else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, 
                 ...)
  atdf <- pretty(index)
  prettydf <- trunc(approx(x = index, y = df, xout = atdf, 
                           rule = 2)$y)
  axis(3, at = atdf, labels = prettydf, cex.axis = 1, tcl = NA)
  if (label) {
    nnz <- length(which)
    xpos <- max(index)
    pos <- 4
    if (xvar == "lambda") {
      xpos <- min(index)
      pos <- 2
    }
    xpos <- rep(xpos, nnz)
    ypos <- beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
} 
