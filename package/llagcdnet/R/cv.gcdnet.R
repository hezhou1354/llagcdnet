cv.gcdnet <- function(x, y, lambda = NULL, pred.loss = c("misclass", 
    "loss"), nfolds = 5, foldid, delta = 2, omega = 0.5, 
    nthresholds = 3, thresholds = NULL, weights = rep(1, nthrs)/nthrs, 
    ...) {
    if (missing(pred.loss)) 
        pred.loss <- "default" else pred.loss <- match.arg(pred.loss)
    N <- nrow(x)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    gcdnet.object <- gcdnet(x, y, lambda = lambda, delta = delta, omega = omega,
                            nthresholds = nthresholds, thresholds = thresholds, weights = weights,
                            ...)
    lambda <- gcdnet.object$lambda
    nthresholds <- gcdnet.object$nthresholds
    thresholds <- gcdnet.object$thresholds
    weights <- gcdnet.object$weights
    # predict -> coef
    nz <- sapply(coef(gcdnet.object, type = "nonzero"), length)
    if (missing(foldid)) 
        foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist <- as.list(seq(nfolds))
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
        which <- foldid == i
        y_sub <- y[!which]
        outlist[[i]] <- gcdnet(x = x[!which, , drop = FALSE], y = y_sub, 
                               lambda = lambda, delta = delta, omega = omega, 
                               nthresholds = nthresholds, thresholds = thresholds, weights = weights,
                               ...)
    }
    ###What to do depends on the pred.loss and the model fit
    fun <- paste("cv", class(gcdnet.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, 
        pred.loss, delta, omega, nthresholds, thresholds, weights))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + 
        cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, gcdnet.fit = gcdnet.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.gcdnet"
    obj
} 
