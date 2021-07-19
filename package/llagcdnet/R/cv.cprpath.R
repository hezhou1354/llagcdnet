cv.cprpath <- function(outlist, lambda, x, y, foldid, 
    pred.loss, delta, omega,thresholds, weights) {
    typenames <- c(misclass = "Misclassification Error", loss = "Composite Probit Loss")
    if (pred.loss == "default") 
        pred.loss <- "loss"
    if (!match(pred.loss, c("loss"), FALSE)) {
        warning("Only 'loss' available for composite probit regression; 'loss' used")
        pred.loss <- "loss"
    }
    x <- as.matrix(x)
    np <- dim(x)
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    nthrs <- length(thresholds)
    #binary responses generated for composite probit
    cpry <- matrix(0, nrow=nobs, ncol=nthrs)
    for (k in seq(nthrs)){
        cpry[, k] <- ifelse(y > thresholds[k], 1, -1)
    }
    
    prob_min <- 1e-05
    fmax <- qnorm(1-prob_min)
    fmin <- -fmax

    nfolds <- max(foldid)
    predmat <- array(NA, dim = c(length(y), length(lambda), nthrs))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, newx = x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami),] <- preds
        nlams[i] <- nlami
    }
    predmat <- pmin(pmax(predmat, fmin), fmax)
    cvraw <- 2 * cprcls(cpry, predmat, weights)
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 
