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
  lines(time, upper, col = 2)
  lines(time, lower, col = 2)
  NULL
}