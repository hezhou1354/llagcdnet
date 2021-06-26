lla2.gcdnet <- function(x, y, lambda = NULL, a = 3.7, delta = 2, omega = 0.5, 
                       ini.lambda = NULL, ini.pred.loss = c("misclass", "loss"), ini.nfolds = 5, ini.foldid, 
                       ...){
  
  this.call <- match.call()
  
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  x <- as.matrix(x)

  ini.cv.gcdnet.obj <- cv.gcdnet(x, y, lambda = ini.lambda, 
                             pred.loss = ini.pred.loss, nfolds = ini.nfolds, foldid = ini.foldid,
                             delta = delta, omega = omega, ...)
  ini.coef <- coef(ini.cv.gcdnet.obj, lambda="lambda.1se")
  ini.beta <- ini.coef[-1,]
  ini.b0 <- ini.coef[1,]
  gcdnet.obj <- gcdnet(x, y, lambda = lambda, delta = delta, omega = omega,
                       ...)
  lambdas <- gcdnet.obj$lambda
  ws <- matrix(mapply(SCAD_deriv, t(as.matrix(ini.beta) %*% rep(1, times=length(lambdas))), a=a, lambda=lambdas),
               nrow=ncol(x), byrow=TRUE)
  beta <- gcdnet.obj$beta
  b0 <- gcdnet.obj$b0
  for (s in 1:length(lambdas)){
    lambda.s <- lambdas[1:s]
    w <- ws[,s]
    second.obj <- gcdnet(x, y, delta = delta, omega = omega,
                         lambda = lambda.s/lambda.s[s], pf = w,
                         ...)
    w <- sapply(second.obj$beta[,s], SCAD_deriv, a=a, lambda=lambda.s[s])
    third.obj <- gcdnet(x, y, delta = delta, omega = omega,
                        lambda = lambda.s/lambda.s[s], pf = w,
                        ...)
    b0[s] <- third.obj$b0[s]
    beta[,s] <- third.obj$beta[,s]
    
  }
  df <- apply(abs(beta) > 0, 2, sum)
  out <- list(b0 = b0, beta = beta, 
              df = df,
              dim = gcdnet.obj$dim,
              lambda = lambdas,
              call = this.call)
  out$a <- a
  class(out) <- c("lla2.gcdnet", class(gcdnet.obj)[2])
  out
}
