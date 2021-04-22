library("llagcdnet")
data(FHT)
set.seed(2011)
x <- FHT$x
y <- FHT$y
nobs <- nrow(x)
nvars <- ncol(x)
cv <- cv.gcdnet(x, y, method ="probit",
               lambda2=0.1, pred.loss="loss", nfolds=5)
plot(cv)
start.time <- Sys.time()
gcdnet.obj <- gcdnet(x, y, nlambda=100, method="probit")
end.time <- Sys.time()
end.time - start.time


gcdnet.beta <- gcdnet.obj$beta[,30]
gcdnet.beta0 <- gcdnet.obj$b0[30]
eta <- y * (gcdnet.beta0 + x%*%gcdnet.beta)
## KKT for beta0
round(sum(y*dnorm(eta)/pnorm(eta)), digits=3)
## KKT for beta
kkt = c()
for (i in 1:nvars){
  kkt[i] <- sum(y*x[,i]*dnorm(eta)/pnorm(eta))/(nobs*gcdnet.obj$lambda[30])
}
round(kkt,digits=2)




