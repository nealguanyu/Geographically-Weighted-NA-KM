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