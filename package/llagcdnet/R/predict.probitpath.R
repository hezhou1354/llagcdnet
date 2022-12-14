predict.probitpath <- function(object, newx, s = NULL, 
    type = c("class", "link"), ...) {
    type <- match.arg(type)
    b0 <- object$b0
    nbeta <- rbind2(b0, object$beta)
    rownames(nbeta)[1] <- "(Intercept)"
    if (!is.null(s)) {
        vnames <- dimnames(nbeta)[[1]]
        dimnames(nbeta) <- list(NULL, NULL)
        lambda <- object$lambda
        lamlist <- lambda.interp(lambda, s)
        nbeta=nbeta[,lamlist$left,drop=FALSE]%*%Diagonal(x=lamlist$frac)+nbeta[,lamlist$right,drop=FALSE]%*%Diagonal(x=1-lamlist$frac)
        dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
    }
    nfit <- as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
    switch(type, link = nfit, class = ifelse(nfit > 0, 1, -1))
} 
