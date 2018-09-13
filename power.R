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

nrep <- 1000 # replicates
k = 1000 # sample size
set.seed(1234)
#b=2 # band width
b = 100 # band width
x.coord <- runif(k, 0, 100)
y.coord <- runif(k, 0, 100)
xy <- cbind(x.coord, y.coord)
# p=c(35,65)
p = c(10, 25)
# p=c(35,65)
#p=c(50,50)
distance <- sqrt(colSums((t(xy) - p) ^ 2))
weight <- exp(-distance / b)
pvalue = rep(NA, nrep)
stat = rep(NA, nrep)
# para=1/(1000+x.coord)
#para=1/(500+20*x.coord)
para <- 1 / 1000

qrep <- 20
q90 <- numeric(qrep)
for (i in 1:qrep) {
  Stime = rexp(k, para)
  Ctime = rexp(k, para)
  Time <- pmin(Stime, Ctime)
  q90[i] <- quantile(Time, probs = .90)
}
tau <- median(q90)
print(tau)

for (i in 1:nrep) {
  Stime = rexp(k, para)
  Ctime = pmin(tau, rexp(k, para))
  Time <- pmin(Stime, Ctime)
  Delta <- as.numeric(Stime <= Ctime)
  D <- data.frame(Time, Delta, weight)
  D <- D[order(Time), ]
  result = HT(D$Time, D$weight, D$Delta)
  pvalue[i] = result[4]
  stat[i] = result[3]
}
qqnorm(stat)
qqline(stat, col = 2)
abline(0, 1, col = "red")
sum(pvalue < 0.05)