cv.bandwidth <- function(z, y, cs, nfolds=5, foldid){
  N <- length(z)
  sigma.hat.sq <- mean((z[order(y)[2:N]] - z[order(y)[1:(N-1)]])^2)/2

  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = N)) 
  } else {
    nfolds <- max(foldid)
  } 
  
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  
  predmat <- matrix(NA, length(y), length(cs))
  for (s in 1:length(cs)){
    c <- cs[s]
    hr <- c * (sigma.hat.sq/N)^(1/5)
    hd <- hr^3
    for (i in seq(nfolds)){
      which <- foldid == i
      y_sub <- y[!which]
      z_sub <-  z[!which, , drop = FALSE]
      mon <- monreg(y_sub, z_sub, hd=hd, hr=hr, inverse=0,
                    t = y[which])
      predmat[which, s] <- mon$estimation
    }
  }
  y <- as.double(y)
  cvraw <- (predmat - y)^2
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  cvname <- "Least Square Loss"
  out <- list(c = cs, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
              cvlo = cvm - cvsd, name = cvname)
  lamin <- getmin(cs, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  obj
}
