cv.lla.probitpath <- function(predmat, x, y, foldid, 
                              pred.loss, delta, omega){
  typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
  if (pred.loss == "default") 
    pred.loss <- "loss"
  if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
    warning("Only 'misclass' and 'loss' available for probit regression; 'loss' used")
    pred.loss <- "loss"
  }
  prob_min <- 1e-05
  fmax <- qnorm(1-prob_min)
  fmin <- -fmax
  ###Turn y into c(0,1)
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
  predmat <- pmin(pmax(predmat, fmin), fmax)
  cvraw <- switch(pred.loss, loss = 2 * probitcls(y * predmat), 
                  misclass = (y != ifelse(predmat > 0, 1, -1)))
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
  
  
}
