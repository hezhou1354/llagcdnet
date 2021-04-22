cv.lla.gcdnet <- function(x, y, a = 3.7, lambda = NULL, pred.loss = c("misclass", "loss"), 
                      nfolds = 5, foldid, delta = 2, omega = 0.5, ...) {
  if (missing(pred.loss)) {
    pred.loss <- "default" 
  } else {
    pred.loss <- match.arg(pred.loss)
  }
  N <- nrow(x)
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  # predict -> coef
  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = N)) 
  } else {
    nfolds <- max(foldid)
  } 
  
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  
  lla.gcdnet.object <- lla.gcdnet(x, y, lambda = lambda, delta = delta, 
                                  omega = omega, ...)
  lambda <- lla.gcdnet.object$lambda
  nz <- lla.gcdnet.object$df
  outlist <- as.list(seq(nfolds))
  ###Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    which <- foldid == i
    y_sub <- y[!which]
    outlist[[i]] <- lla.gcdnet(x = x[!which, , drop = FALSE], 
                           y = y_sub, lambda = lambda, delta = delta, omega = omega, ...)
  }
  
  ###What to do depends on the pred.loss and the model fit
  fun <- paste("cv", class(lla.gcdnet.object)[[2]], sep = ".")
  cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, 
                               pred.loss, delta, omega))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambda, a = a, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
              cvlo = cvm - cvsd, nzero = nz, name = cvname, 
              lla.gcdnet.fit = lla.gcdnet.object)
  lamin <- getmin(lambda, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.lla.gcdnet"
  obj
} 
