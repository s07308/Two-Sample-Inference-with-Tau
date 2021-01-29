tau_fixed <- function(x, y, na.rm = TRUE) {
  if(na.rm) {
    x <- na.omit(x)
    y <- na.omit(y)
  }
  
  n1 <- length(x)
  n0 <- length(y)
  
  MW.score <- 0
  for(i in 1:n1) {
    MW.score <- MW.score + sum(ifelse(x[i] > y, 1, 0) + ifelse(x[i] < y, -1, 0))
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
              var.tau.general = var.tau.general,
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
    MW.score <- MW.score + sum(ifelse(x[i] > y, 1, 0) + ifelse(x[i] < y, -1, 0))
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

tau_ipcw_res <- function(Y1, delta, X1, t.star) {
  obs.0 <- Y1[X1 == 0]
  obs.1 <- Y1[X1 == 1]
  delta.0 <- delta[X1 == 0]
  delta.1 <- delta[X1 == 1]
  n0 <- sum(1 - X1)
  n1 <- sum(X1)
  
  order.matrix <- ifelse(sapply(obs.1, ">", obs.0), 1, -1)
  delta0.matrix <- matrix(rep(delta.0, n1), ncol = n1)
  delta1.matrix <- matrix(rep(delta.1, n0), ncol = n1, byrow = TRUE)
  orderable <- ifelse(sapply(obs.1, ">", obs.0), delta0.matrix, delta1.matrix)
  
  Y0.matrix <- matrix(rep(obs.0, n1), ncol = n1)
  Y1.matrix <- matrix(rep(obs.1, n0), ncol = n1, byrow = TRUE)
  Y.min <- ifelse(sapply(obs.1, ">", obs.0), Y0.matrix, Y1.matrix)
  
  G.fit <- survfit(Surv(time = Y1, event = 1 - delta) ~ 1)
  G_func <- stepfun(x = G.fit$time, y = c(1, G.fit$surv))
  G <- matrix(G_func(Y.min), ncol = ncol(Y.min))
  
  tau.res <- sum(ifelse(Y.min <= t.star, order.matrix, 0) * orderable * G ^ (-2), na.rm = TRUE) / n0 / n1
  tau.res.rewei <- tau.res * n0 * n1 / sum(ifelse(Y.min <= t.star, abs(order.matrix), 0) * orderable * G ^ (-2), na.rm = TRUE)
  
  return(list(tau.res = tau.res,
              tau.res.rewei = tau.res.rewei))
}

tail_est <- function(Y1, delta, X1, t.star, tail.dist) {
  ## tail estimation
  
  if(tail.dist == "exp") {
    exp.fit <- fit_data(data.frame(time = Y1, censor = delta, X1), dist = "exp", by = "X1")
    rate.0.est <- exp.fit[[1]]$estimate
    rate.1.est <- exp.fit[[2]]$estimate
    
    fc.est <- function(t) return((1 - pexp(t, rate = rate.1.est)) * dexp(t, rate = rate.0.est))
    fd.est <- function(t) return((1 - pexp(t, rate = rate.0.est)) * dexp(t, rate = rate.1.est))
    
    S1 <- 1 - pexp(t.star, rate.1.est)
    S0 <- 1 - pexp(t.star, rate.0.est)
  }
  
  if(tail.dist == "weibull") {
    weibull.fit <- fit_data(data.frame(time = Y1, censor = delta, X1), dist = "weibull", by = "X1")
    param.0.est <- weibull.fit[[1]]$estimate
    param.1.est <- weibull.fit[[2]]$estimate
    
    fc.est <- function(t) {
      f <- dweibull(t, shape = param.0.est["shape"], scale = param.0.est["scale"])
      S <- 1 - pweibull(t, shape = param.1.est["shape"], scale = param.1.est["scale"])
      return(S * f)
    }
    
    fd.est <- function(t) {
      f <- dweibull(t, shape = param.1.est["shape"], scale = param.1.est["scale"])
      S <- 1 - pweibull(t, shape = param.0.est["shape"], scale = param.0.est["scale"])
      
      return(S * f)
    }
    
    S1 <- 1 - pweibull(t.star, shape = param.1.est["shape"], scale = param.1.est["scale"])
    S0 <- 1 - pweibull(t.star, shape = param.0.est["shape"], scale = param.0.est["scale"])
  }
  
  if(tail.dist == "lnorm") {
    lnorm.fit <- fit_data(data.frame(time = Y1, censor = delta, X1), dist = "lnorm", by = "X1")
    param.0.est <- lnorm.fit[[1]]$estimate
    param.1.est <- lnorm.fit[[2]]$estimate
    
    fc.est <- function(t) {
      f <- dlnorm(t, meanlog = param.0.est["meanlog"], sdlog = param.0.est["sdlog"])
      S <- 1 - plnorm(t, meanlog = param.1.est["meanlog"], sdlog = param.1.est["sdlog"])
      return(S * f)
    }
    
    fd.est <- function(t) {
      f <- dlnorm(t, meanlog = param.1.est["meanlog"], sdlog = param.1.est["sdlog"])
      S <- 1 - plnorm(t, meanlog = param.0.est["meanlog"], sdlog = param.0.est["sdlog"])
      
      return(S * f)
    }
    
    S1 <- 1 - plnorm(t.star, meanlog = param.1.est["meanlog"], sdlog = param.1.est["sdlog"])
    S0 <- 1 - plnorm(t.star, meanlog = param.0.est["meanlog"], sdlog = param.0.est["sdlog"])
  }
  
  if(tail.dist == "logis") {
    logis.fit <- fit_data(data.frame(time = Y1, censor = delta, X1), dist = "logis", by = "X1")
    param.0.est <- logis.fit[[1]]$estimate
    param.1.est <- logis.fit[[2]]$estimate
    
    fc.est <- function(t) {
      f <- dlogis(t, location = param.0.est["location"], scale = param.0.est["scale"])
      S <- 1 - plogis(t, location = param.1.est["location"], scale = param.1.est["scale"])
      return(S * f)
    }
    
    fd.est <- function(t) {
      f <- dlogis(t, location = param.1.est["location"], scale = param.1.est["scale"])
      S <- 1 - plogis(t, location = param.0.est["location"], scale = param.0.est["scale"])
      
      return(S * f)
    }
    
    S1 <- 1 - plogis(t.star, location = param.1.est["location"], scale = param.1.est["scale"])
    S0 <- 1 - plogis(t.star, location = param.0.est["location"], scale = param.0.est["scale"])
  }
  
  head.prob <- 1 - S1 * S0
  tail.b.est <- integrate(fc.est, lower = t.star, upper = Inf)$value - integrate(fd.est, lower = t.star, upper = Inf)$value
  
  return(list(tail.b.est = tail.b.est,
              head.prob = head.prob))
}

tau_ipcw2 <- function(Y1, delta, X1, t.star, tail.dist) {
  tau.rewei <- tau_ipcw_res(Y1, delta, X1, t.star)$tau.res.rewei
  tail.result <- tail_est(Y1, delta, X1, t.star, tail.dist)
  
  return(tau.rewei * tail.result$head.prob + tail.result$tail.b.est)
}