Hcov <- function(w1, w2, t1, t2, time, delta) {
  n = length(time)
  tmin = min(t1, t2)
  yw1 = rev(cumsum(rev(w1)))
  yw2 = rev(cumsum(rev(w2)))
  y <- n:1
  w12 <- rev(cumsum(rev(w1 * w2)))
  cov <- cumsum(w12 * delta / (yw1 * yw2 * y))
  cov[findInterval(tmin, time)] # Better: Compute cov matrix for any vector of times.
}
