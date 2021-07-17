cprpath <- function(x, y, nlam, flmin, ulam, isd, 
    eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, nobs, nvars, vnames, 
    nthrs, athrs, wt) {
    #################################################################################
    #data setup
    #binary responses generated for composite probit
    cpry <- matrix(0, nrow=nobs, ncol=nthrs)
    for (k in seq(nthrs)){
        cpry[, k] <- ifelse(y > athrs[k], 1, -1)
    }

    #################################################################################
    # call Fortran core
    fit <- .Fortran("cprlassoNET", lam2, nobs, nvars, as.double(x), 
        as.double(cpry), jd, pf, pf2, dfmax, pmax, nlam, flmin, ulam, 
        eps, isd, maxit, nthrs, wt,
        nalam = integer(1), b0 = double(nthrs * nlam), 
        beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
        alam = double(nlam), npass = integer(1), jerr = integer(1))
    #################################################################################
    # output
    outlist <- cpr_getoutput(fit, maxit, pmax, nvars, vnames, nthrs)
    outlist <- c(outlist, list(nthresholds = nthrs, thresholds = athrs, weights = wt, 
                               npasses = fit$npass, jerr = fit$jerr))
    class(outlist) <- c("cprpath")
    outlist
} 
