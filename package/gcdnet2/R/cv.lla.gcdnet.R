cv.lla.gcdnet <- function(x, y, lambda = NULL, pred.loss = c("misclass", "loss"), 
                      nfolds = 5, foldid, delta = 2, omega = 0.5, ...) {
  if (missing(pred.loss)) {
    pred.loss <- "default" 
  } else {
    pred.loss <- match.arg(pred.loss)
  }
  N <- nrow(x)
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  gcdnet.object <- gcdnet(x, y, lambda = lambda, delta = delta, omega = omega,
                          ...)
  lambdas <- gcdnet.object$lambda
  # predict -> coef
  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = N)) 
  } else {
    nfolds <- max(foldid)
  } 
  
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
    
  nz <- rep(0, times=length(lambdas))
  predmat <- matrix(NA, length(y), length(lambdas))
  outlist <- as.list(lambdas)
    
  ###Now fit the nfold models and store them
  for (s in 1:length(lambdas)){
    lambda <- lambdas[1:s]
    first.obj <- gcdnet(x, y, lambda=lambda, delta = delta, omega = omega, ...)
    w <- sapply(first.obj$beta[,s], SCAD_deriv, a=3.7, lambda=lambda[s])
    second.obj <- gcdnet(x, y, delta = delta, omega = omega,
                         lambda = lambda/lambda[s], lambda2 = 0, pf = w,
                         ...)
    w <- sapply(second.obj$beta[,s], SCAD_deriv, a=3.7, lambda=lambda[s])
    third.obj <- gcdnet(x,y, delta = delta, omega = omega,
                        lambda = lambda/lambda[s], lambda2 = 0, pf = w,
                        ...)
    nz[s] <- length(coef(third.obj, type = "nonzero", s=1))
    outlist[[s]] <- third.obj
    for (i in seq(nfolds)){
      which <- foldid == i
      y_sub <- y[!which]
      x_sub <-  x[!which, , drop = FALSE]
      first.iter <- gcdnet(x_sub, y_sub, lambda=lambda, delta = delta, omega = omega, ...)
      w <- sapply(first.iter$beta[,s], SCAD_deriv, a=3.7, lambda=lambda[s])
      second.iter <- gcdnet(x_sub, y_sub, delta = delta, omega = omega,
                            lambda = lambda/lambda[s], lambda2 = 0, pf = w,
                            ...)
      w <- sapply(second.iter$beta[,s], SCAD_deriv, a=3.7, lambda=lambda[s])
      third.iter <- gcdnet(x_sub, y_sub, delta = delta, omega = omega,
                           lambda = lambda/lambda[s], lambda2 = 0, pf = w,
                           ...)
      predmat[which, s] <- predict(third.iter, x[which, , drop = FALSE],
                                   type = "link", s=1)
    }
  }
  ###What to do depends on the pred.loss and the model fit
  fun <- paste("cv.lla", class(gcdnet.object)[[2]], sep = ".")
  cvstuff <- do.call(fun, list(predmat, x, y, foldid, 
                               pred.loss, delta, omega))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambdas, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
              cvlo = cvm - cvsd, nzero = nz, name = cvname, gcdnet.fit = outlist)
  lamin <- getmin(lambdas, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.lla.gcdnet"
  obj
} 
