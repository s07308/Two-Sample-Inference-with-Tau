tau_fixed <- function(x, y, na.rm = TRUE) {
  if(na.rm) {
    x <- na.omit(x)
    y <- na.omit(y)
  }
  
  n1 <- length(x)
  n0 <- length(y)
  
  MW.score <- 0
  for(i in 1:n1) {
    for(j in 1:n0) {
      MW.score <- MW.score + ifelse(x[i] > y[j], 1, -1)
    }
  }
  
  tau.hat <- MW.score / n1 / n0
  
  Fx <- ecdf(x)
  Fy <- ecdf(y)
  sigma10.hat <- 4 * var(Fy(x))
  sigma01.hat <- 4 * var(Fx(y))
  sigma11.hat <- 4 * 0.5 * 0.5
  sigma11.hat.general <- 1 - tau.hat ^ 2
  
  var.tau.fixed <- (n1 + n0 + 1) / 3 / n0 / n1
  var.tau.0 <- ((n1 - 1) * sigma01.hat + (n0 - 1) * sigma10.hat + sigma11.hat) / n1 / n0
  var.tau.general <- ((n1 - 1) * sigma01.hat + (n0 - 1) * sigma10.hat + sigma11.hat.general) / n1 / n0
  
  z.score.fixed <- tau.hat / sqrt(var.tau.fixed)
  z.score.0 <- tau.hat / sqrt(var.tau.0)
  
  p.value.fixed <- 1 - pnorm(abs(z.score.fixed))
  p.value.0 <- 1 - pnorm(abs(z.score.0))
  ci <- c(tau.hat + qnorm(0.025) * sqrt(var.tau.general),
          tau.hat + qnorm(0.975) * sqrt(var.tau.general))
  
  return(list(tau = tau.hat,
              var.tau.fixed = var.tau.fixed,
              var.tau.0 = var.tau.0,
              z.score.fixed = z.score.fixed,
              z.score.0 = z.score.0,
              p.value.fixed = p.value.fixed,
              p.value.0 = p.value.0,
              ci = ci))
}

tau_random <- function(x, y, na.rm = TRUE) {
  if(na.rm) {
    x <- na.omit(x)
    y <- na.omit(y)
  }
  
  n1 <- length(x)
  n0 <- length(y)
  size <- n1 + n0
  
  MW.score <- 0
  for(i in 1:n1) {
    for(j in 1:n0) {
      MW.score <- MW.score + ifelse(x[i] > y[j], 1, -1)
    }
  }
  
  tau.hat <- MW.score / n1 / n0
  
  X1 <- c(rep(1, length(x)), rep(0, length(y)))
  T1 <- c(x, y)
  
  ## Random Variance
  p.hat <- mean(X1 == 1)
  q.hat <- 1 - p.hat
  F1 <- ecdf(x)
  F0 <- ecdf(y)
  eta <- (X1 == 1) * (1 - p.hat) * (2 * F0(T1) - 1) + (X1 == 0) * p.hat * (1 - 2 * F1(T1))
  sigma1 <- var(eta)
  adjust <- p.hat ^ 2 * (1 - p.hat) ^ 2
  
  var.tau.random <- (3 * size * p.hat * q.hat) ^ (-1)
  var.tau <- sigma1 / size / adjust
  
  z.score.random <- tau.hat / sqrt(var.tau.random)
  z.score.0 <- tau.hat / sqrt(var.tau)
  
  p.value.random <- 1 - pnorm(abs(z.score.random))
  p.value.0 <- 1 - pnorm(abs(z.score.0))
  
  ci <- c(tau.hat + qnorm(0.025) * sqrt(var.tau),
          tau.hat + qnorm(0.975) * sqrt(var.tau))
  
  return(list(tau = tau.hat,
              var.tau.random = var.tau.random,
              var.tau = var.tau,
              z.score.random = z.score.random,
              z.score.0 = z.score.0,
              p.value.random = p.value.random,
              p.value.0 = p.value.0,
              ci = ci))
}

tau_ipcw <- function(Y1, delta, X1) {
  size <- length(Y1)
  T1 <- Y1
  df.pair <- matrix(nrow = choose(n = size, k = 2), ncol = 8)
  k <- 1
  for(i in seq(1, size - 1)) {
    for(j in seq(i + 1, size)) {
      df.pair[k, ] <- c(X1[i], Y1[i], delta[i], T1[i], X1[j], Y1[j], delta[j], T1[j])
      k <- k + 1
    }
  }
  
  G.fit <- survfit(Surv(time = Y1, event = 1 - delta) ~ 1)
  G_func <- stepfun(x = G.fit$time, y = c(1, G.fit$surv))
  G.hat.1 <- G_func(df.pair[, 2])
  G.hat.2 <- G_func(df.pair[, 6])
  
  c.hat.1 <- (df.pair[, 2] > df.pair[, 6]) & (df.pair[, 1] > df.pair[, 5]) & (df.pair[, 7] == 1)
  c.hat.1 <- c.hat.1 / (G.hat.2 ^ 2)
  c.hat.1 <- ifelse(is.nan(c.hat.1), 0, c.hat.1)
  
  c.hat.2 <- (df.pair[, 2] < df.pair[, 6]) & (df.pair[, 1] < df.pair[, 5]) & (df.pair[, 3] == 1)
  c.hat.2 <- c.hat.2 / (G.hat.1 ^ 2)
  c.hat.2 <- ifelse(is.nan(c.hat.2), 0, c.hat.2)
  
  d.hat.1 <- (df.pair[, 2] > df.pair[, 6]) & (df.pair[, 1] < df.pair[, 5]) & (df.pair[, 7] == 1)
  d.hat.1 <- d.hat.1 / (G.hat.2 ^ 2)
  d.hat.1 <- ifelse(is.nan(d.hat.1), 0, d.hat.1)
  
  d.hat.2 <- (df.pair[, 2] < df.pair[, 6]) & (df.pair[, 1] > df.pair[, 5]) & (df.pair[, 3] == 1)
  d.hat.2 <- d.hat.2 / (G.hat.1 ^ 2)
  d.hat.2 <- ifelse(is.nan(d.hat.2), 0, d.hat.2)
  
  c.hat <- sum(c.hat.1 + c.hat.2)
  d.hat <- sum(d.hat.1 + d.hat.2)
  
  n0 <- sum(X1 == 0)
  n1 <- sum(X1 == 1)
  
  return((c.hat - d.hat) / n0 / n1)
}