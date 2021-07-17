hubercls <- function(r, delta) {
  (1 - r - 0.5 * delta) * (r <= 1 - delta) + 0.5 * (1 - r)^2 / delta * (r <= 1) * (r > 1 - delta)
} 

ercls <- function(r, omega) {
  abs(omega - (r < 0)) * r^2
} 

probitcls <- function(r) {
    -log(pnorm(r))
}

cprcls <- function(cpry, predmat, weights){
  dd <- dim(predmat)
  nthrs <- dd[3]
  nlams <- dd[2]
  nobs <- dd[1]
  cprmat <- array(NA, dim = dd)
  for (i in seq(nthrs)){
    cprmat[,,i] <- probitcls((cpry[,i] %*% t(rep(1, times=nlams))) * predmat[,,i])
  }
  apply(cprmat, c(1,2), weighted.mean, w=weights)
}
