predict.cprpath <- function(object, newx, s = NULL, 
    type = c("class", "link"), ...) {
    type <- match.arg(type)
    b0 <- object$b0
    beta <- object$beta
    nbeta <- rbind2(b0, beta)
    nthrs <- nrow(b0)
    rownames(nbeta)[seq(nthrs)] <- paste("(Intercept", seq(nthrs), ")", sep = "")
    if (!is.null(s)) {
        vnames <- dimnames(nbeta)[[1]]
        dimnames(nbeta) <- list(NULL, NULL)
        lambda <- object$lambda
        lamlist <- lambda.interp(lambda, s)
        nbeta=nbeta[,lamlist$left,drop=FALSE]%*%Diagonal(x=lamlist$frac)+nbeta[,lamlist$right,drop=FALSE]%*%Diagonal(x=1-lamlist$frac)
        dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
        b0 <- as.matrix(nbeta[seq(nthrs),])
        beta <- as.matrix(nbeta[-seq(nthrs),])
    }
    nmat <- as.matrix(newx) %*% beta
    np <- dim(nmat)
    nfit <- as.list(seq(nthrs))
    for (i in seq(nthrs)){
        nfit[[i]] <- as.matrix(nmat - matrix(1, nrow=np[1], ncol=1) %*% b0[i,])
    }
    nfit <- abind(nfit, along = 3)
    switch(type, link = nfit, class = ifelse(nfit > 0, 1, -1))
} 
