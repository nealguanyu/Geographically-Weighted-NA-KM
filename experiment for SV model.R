Hest <- function(time, weight, delta, timrec = NULL) {
  # Computes weighted Nelson-Aalen and its variance estimate.
  # Also extracts results for times in timrec.
  # Assume ordered by time already, from smallest to largest.
  n <- length(time)
  if (is.null(weight))
    weight <- rep(1, n) # Compute global estimate.
  yw <- rev(cumsum(rev(weight)))
  v <- rev(cumsum(rev(weight ^ 2)))
  Hest <- cumsum(jump <- weight * delta / yw) # Nelson-Aalen estimate
  Vest <- cumsum(v / yw ^ 2 * jump) # Variance estimate.
  ans <- list(time = time, H = Hest, V = Vest)
  if (length(timrec) > 0) {
    ans <- c(ans, list(extract = extract(time, cbind(
      H = Hest, V = Vest
    ), timrec)))
  }
  ans
}



plotband <- function(Hout, alpha = .95) {
  # Plot estimate and confidence band.
  # Hout is the list produced by Hest().
  z <- qnorm((1 + alpha) / 2)
  H <- Hout$H
  time <- Hout$time
  width <- z * sqrt(Hout$V)
  upper <- H + width
  lower <- H - width
  plot(
    c(time, time),
    c(upper, lower),
    type = "n",
    xlab = "time",
    ylab = "cum. hazard"
  )
  lines(time, H)
  lines(time, upper, col = 2, lty = 2)
  lines(time, lower, col = 2, lty = 2)
  NULL
}
k = 1000 # number of observations
b = 50
l = 200 # location
set.seed(1234)

x.coord <- runif(k, 0, 1000)
y.coord <- runif(k, 0, 1000)
Stime = rep(0, k)
Ctime = rep(0, k)
for (i in 1:k) {
  Stime[i] = rexp(1, 1 / (1000 + x.coord[i]))
  Ctime[i] = rexp(1, 1 / (1000 + x.coord[i]))
}# spatial varying distribution
Time <- pmin(Stime, Ctime)
Delta <- as.numeric(Stime <= Ctime)
distance = rep(0, k)
xy <- cbind(x.coord, y.coord)
distance <- sqrt(colSums((t(xy) - xy[l, ]) ^ 2))

weight <- exp(-distance / b)
D = data.frame(Time, Delta, weight)

D = D[order(Time), ]
H = Hest(D$Time, D$weight, D$Delta)
plotband(H)
lines(sort(Time), sort(Time) * 1 / (1000 + x.coord[l]), col = '4')