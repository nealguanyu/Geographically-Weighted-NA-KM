Test.multi <- function(time, delta, w, A, tvec) {
  # For each time in tvec and each contrast (given by the columns of
  #   A) among the weighted cumulative hazard estimates (computed
  #   using the columns of w), a test of the null hypothesis that the
  #   true contrast is zero is carried out. A z-score is given for each
  #   contrast and each time in tvec.
  # The variance of the test statistics is computed under the null hypothesis
  # that there is no spatial variation in the survival time distribution.
  # This function sorts by time first.  It does NOT assume the data
  #   is already sorted.
  # time and delta are vectors of length n.
  # w is an n x L matrix of weights.
  # L is the number of locations, n is the number of observations.
  # Each column of w gives the weights for a particular location (perhaps).
  # A is an L x nc matrix, each column giving a contrast among the L
  #   locations.
  # tvec is a vector of times at which the test is carried out.
  library('matrixStats')
  library('MASS')
  tvec <- sort(tvec)
  n <- length(time)
  L <- ncol(w)
  p <- sort.list(time, decreasing = TRUE)
  # Note: It is easier to work with decreasing time than to
  # always have to use rev().
  time <- time[p]
  delta <- delta[p]
  w <- w[p, ]
  rr <- n:1
  yw <- apply(w, 2, cumsum)
  test <- apply(((w * delta) / yw)[rr, ], 2, cumsum) %*% A
  # Note: In the above, time is put back in the proper order.
  ww <- array(rep(w, L), c(n, L, L))
  yww <- apply((ww * aperm(ww, c(1, 3, 2))), c(2, 3), cumsum)
  yw2 <- array(rep(yw, L), c(n, L, L))
  yw2 <- yw2 * aperm(yw2, c(1, 3, 2))
  cv <- apply((yww / yw2 * (delta / (1:n)))[rr, , ], c(2, 3), cumsum)
  it <- findInterval(tvec, rev(time))
  var = array(0, c(ncol(A), ncol(A), length(it)))
  for (i in 1:ncol(A)) {
    for (j in 1:ncol(A)) {
      temp <- rowSums(cv * array(rep(outer(A[, i], A[, j]), each = n), c(n, L, L)))
      temp <- temp[it]
      var[i, j, ] = temp
    }
  }
  stat = rep(0, length(it))
  test <- test[it, ]
  test[2:length(it), ] = colDiffs(test)
  for (i in 2:length(it)) {
    var[, , i] = var[, , i] - var[, , i - 1]
  }
  for (i in 1:length(it)) {
    stat[i] = test[i, ] %*% ginv(var[, , i]) %*% test[i, ]
  }
  return(sum(stat))
}
