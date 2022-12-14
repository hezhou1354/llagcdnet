lla.gcdnet <- function(x, y, lambda = NULL, a = 3.7, delta = 2, omega = 0.5, 
                       thresholds = NULL, ...){
  
  this.call <- match.call()
  
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  x <- as.matrix(x)
  gcdnet.object <- gcdnet(x, y, lambda = lambda, delta = delta, omega = omega,
                          thresholds = thresholds, ...)
  lambdas <- gcdnet.object$lambda
  thresholds <- gcdnet.object$thresholds
  ws <- matrix(mapply(SCAD_deriv, t(as.matrix(gcdnet.object$beta)), a=a, lambda=lambdas),
              nrow=ncol(x), byrow=TRUE)
  beta <- gcdnet.object$beta
  b0 <- matrix(gcdnet.object$b0, ncol=length(lambdas))
  colnames(b0) <- colnames(gcdnet.object$beta)
  for (s in 1:length(lambdas)){
    lambda <- lambdas[1:s]
    w <- ws[,s]
    second.obj <- gcdnet(x, y, delta = delta, omega = omega,
                         lambda = lambda/lambda[s], pf = w,
                         thresholds = thresholds, ...)
    w <- sapply(second.obj$beta[,s], SCAD_deriv, a=a, lambda=lambda[s])
    third.obj <- gcdnet(x,y, delta = delta, omega = omega,
                        lambda = lambda/lambda[s], pf = w,
                        thresholds = thresholds, ...)
    b0[,s] <- matrix(third.obj$b0, ncol=s)[,s]
    beta[,s] <- third.obj$beta[,s]
    
  }
  df <- apply(abs(beta) > 0, 2, sum)
  out <- list(b0 = b0, beta = beta, 
              df = df,
              dim = gcdnet.object$dim,
              lambda = lambdas,
              nthresholds = length(thresholds),
              thresholds = thresholds,
              weights = gcdnet.object$weights,
              call = this.call)
  out$a <- a
  class(out) <- c("lla.gcdnet", class(gcdnet.object)[2])
  out
}
   