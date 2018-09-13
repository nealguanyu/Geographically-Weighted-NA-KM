HT <- function(time, weight, delta) {
  # Hypothesis test for difference between global and weighted
  #   cumulative hazard.  This has the correct variance formula.
  # Assumes data ordered by time already, from smallest to largest.
  n <- length(time)
  y <- n:1
  yw <- rev(cumsum(rev(weight)))
  gw <- weight / yw
  g <- 1 / y
  gdiff <- gw - g
  stat <- sum(gdiff * delta)
  # Now compute estimated variance of stat.
  vv <- matrix(0, n, n)
  vv[] <- col(vv) <= row(vv) # Lower triangular matrix, all 1's.
  bw <- scale(weight * vv, center = FALSE, scale = yw)
  b <- scale(vv, center = FALSE, scale = y)
  var <- sum((bw - b) ^ 2 %*% as.vector((delta / y)))
  score <- stat / sqrt(var)
  pvalue <- (1 - pnorm(abs(score))) * 2
  c(
    stat = stat,
    var = var,
    score = score,
    pvalue = pvalue
  )
}