tau.bar_func <- function(X, observed.time) {
  n <- length(observed.time)
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  observed.time.1 <- observed.time[X == 1]
  
  ## statistics 
  num.conc <- sum(sapply(observed.time.1, ">", observed.time.0))
  num.disc <- sum(sapply(observed.time.1, "<", observed.time.0))
  
  U <- num.conc - num.disc
  tau.bar <- U / N0 / N1
  
  ## variance estimates
  var.fixed <- tau.bar.var.fixed(X = X, observed.time = observed.time, tau.bar.b = tau.bar)
  var.random <- tau.bar.var.random(X = X, observed.time = observed.time, tau.bar.b = tau.bar)
  
  ## varinace estimates under null
  p1.est <- mean(X)
  p0.est <- 1 - p1.est
  var.null.0 <- 1 / (3 * n * p0.est * p1.est)
  var.null.tau <- tau.bar.var.random(X = X, observed.time = observed.time, tau.bar.b = 0)
  
  ## confidence interval
  ci.fixed <- tau.bar + qnorm(c(0.025, 0.975)) * sqrt(var.fixed)
  ci.random <- tau.bar + qnorm(c(0.025, 0.975)) * sqrt(var.random)
  
  ## standardized statistics under null
  z.val.0 <- tau.bar / sqrt(var.null.0)
  z.val.tau <- tau.bar / sqrt(var.null.tau)
  
  ## p-values
  p.val.0 <- 1 - pchisq(z.val.0 ^ 2, df = 1)
  p.val.tau <- 1 - pchisq(z.val.tau ^ 2, df = 1)
  
  return(list(tau.bar = tau.bar,
              var.fixed = var.fixed,
              var.random = var.random,
              var.null.0 = var.null.0,
              var.null.tau = var.null.tau,
              ci.fixed = ci.fixed,
              ci.random = ci.random,
              z.val.0 = z.val.0,
              z.val.tau = z.val.tau,
              p.val.0 = p.val.0,
              p.val.tau = p.val.tau))
}

tau.hat_func <- function(X, observed.time, delta) {
  n <- length(observed.time)
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  delta.0 <- delta[X == 0]
  observed.time.1 <- observed.time[X == 1]
  delta.1 <- delta[X == 1]
  
  ## estiamte the survival functions of censoring variable
  G0.fit <- survfit(Surv(observed.time.0, 1 - delta.0) ~ 1)
  G1.fit <- survfit(Surv(observed.time.1, 1 - delta.1) ~ 1)
  G0.fit_func <- stepfun(G0.fit$time, c(1, G0.fit$surv))
  G1.fit_func <- stepfun(G1.fit$time, c(1, G1.fit$surv))
  
  ## all combinations
  obs.time.mat.0 <- matrix(observed.time.0, nrow = N0, ncol = N1)
  obs.time.mat.1 <- matrix(observed.time.1, nrow = N0, ncol = N1, byrow = TRUE)
  delta.mat.0 <- matrix(as.logical(delta.0), nrow = N0, ncol = N1)
  delta.mat.1 <- matrix(as.logical(delta.1), nrow = N0, ncol = N1, byrow = TRUE)
  
  ## tilde{Y}
  min.index <- obs.time.mat.0 <= obs.time.mat.1
  obs.min.time.mat <- ifelse(min.index, obs.time.mat.0, obs.time.mat.1)
  
  min.index.1 <- obs.time.mat.1 <= obs.time.mat.0
  obs.min.time.mat.1 <- ifelse(min.index.1, obs.time.mat.1, obs.time.mat.0)
  
  ## orderable indicator
  orderable.indicator <- ifelse(min.index, delta.mat.0, delta.mat.1)
  orderable.indicator <- ifelse(obs.time.mat.0 == obs.time.mat.1 & delta.mat.0 == delta.mat.1,
                                FALSE,
                                orderable.indicator)
  
  ## G0.hat and G1.hat at tilde{Y}
  G0.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G0.fit_func)
  G1.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G1.fit_func)
  
  ## concordance and discordance
  conc <- ifelse(min.index, 1, -1)
  
  ## statistics
  U <- ifelse(orderable.indicator, conc / G0.val / G1.val, 0)
  tau.hat <- sum(U) / N0 / N1
  
  ## variance estimates
  var.1.random <- tau.hat.var1.random(X = X, observed.time = observed.time, delta = delta, tau.hat = tau.hat)
  var.2.random <- tau.hat.var2.random(X = X, observed.time = observed.time, delta = delta)
  var.3.random <- tau.hat.var3.random(X = X, observed.time = observed.time, delta = delta, tau.hat = tau.hat)
  var.random <- (var.3.random - var.2.random - var.1.random) / n
  
  var.2.fixed <- tau.hat.var2.fixed(X = X, observed.time = observed.time, delta = delta)
  var.3.fixed <- tau.hat.var3.fixed(X = X, observed.time = observed.time, delta = delta, tau.hat = tau.hat)
  var.fixed <- (var.3.fixed - var.2.fixed) / n
  
  ## variance estimates under null
  var.3.null.tau <- tau.hat.var3.random(X = X, observed.time = observed.time, delta = delta, tau.hat = 0)
  var.null.tau <- (var.3.null.tau - var.2.random) / n
  
  ## confidence interval
  ci.fixed <- tau.hat + qnorm(c(0.025, 0.975)) * sqrt(var.fixed)
  ci.random <- tau.hat + qnorm(c(0.025, 0.975)) * sqrt(var.random)
  
  ## standardized statistics under null
  z.val.tau <- tau.hat / sqrt(var.null.tau)
  
  ## p-value
  p.val.tau <- 1 - pchisq(z.val.tau ^ 2, df = 1)
  
  return(list(tau.hat = tau.hat,
              U = sum(U),
              var.fixed = var.fixed,
              var.random = var.random,
              var.null.tau = var.null.tau,
              ci.fixed = ci.fixed,
              ci.random = ci.random,
              z.val.tau = z.val.tau,
              p.val.tau = p.val.tau))
}

imputed.tau.hat_func <- function(X, observed.time, delta, t.star) {
  restricted.tau.hat <- restricted.tau.hat_func(X = X, observed.time = observed.time, delta = delta, t.star = t.star)
  
  ## Weibull
  ## estimation distribution (Weibull)
  weibull.fit <- fit_data(data.frame(time = observed.time, censor = delta, X1 = X), dist = "weibull", by = "X1")
  param.0.est <- weibull.fit[[1]]$estimate
  param.1.est <- weibull.fit[[2]]$estimate
  
  ## tail approximate
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
  
  tail.est <- integrate(fc.est, lower = t.star, upper = Inf)$value - integrate(fd.est, lower = t.star, upper = Inf)$value
  imputed.tau.hat.weibull <- restricted.tau.hat$tau.hat + tail.est
  
  ## exponential
  ## estimation distribution (Exponential)
  exp.fit <- fit_data(data.frame(time = observed.time, censor = delta, X1 = X), dist = "exp", by = "X1")
  rate.0.est <- exp.fit[[1]]$estimate
  rate.1.est <- exp.fit[[2]]$estimate
  
  ## tail approximate
  fc.est <- function(t) return((1 - pexp(t, rate = rate.1.est)) * dexp(t, rate = rate.0.est))
  fd.est <- function(t) return((1 - pexp(t, rate = rate.0.est)) * dexp(t, rate = rate.1.est))
  
  tail.est <- integrate(fc.est, lower = t.star, upper = Inf)$value - integrate(fd.est, lower = t.star, upper = Inf)$value
  imputed.tau.hat.exp <- restricted.tau.hat$tau.hat + tail.est
  
  ## log-normal
  ## estimation distribution (log-normal)
  lnorm.fit <- fit_data(data.frame(time = observed.time, censor = delta, X1 = X), dist = "lnorm", by = "X1")
  param.0.est <- lnorm.fit[[1]]$estimate
  param.1.est <- lnorm.fit[[2]]$estimate
  
  ## tail approximate
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
  
  tail.est <- integrate(fc.est, lower = t.star, upper = Inf)$value - integrate(fd.est, lower = t.star, upper = Inf)$value
  imputed.tau.hat.lnorm <- restricted.tau.hat$tau.hat + tail.est
  
  ## logistic
  ## estimation distribution (logistic)
  logis.fit <- fit_data(data.frame(time = observed.time, censor = delta, X1 = X), dist = "logis", by = "X1")
  param.0.est <- logis.fit[[1]]$estimate
  param.1.est <- logis.fit[[2]]$estimate
  
  ## tail approximate
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
  
  tail.est <- integrate(fc.est, lower = t.star, upper = Inf)$value - integrate(fd.est, lower = t.star, upper = Inf)$value
  imputed.tau.hat.logis <- restricted.tau.hat$tau.hat + tail.est
  
  return(list(weibull = imputed.tau.hat.weibull,
              exp = imputed.tau.hat.exp,
              lnorm = imputed.tau.hat.lnorm,
              logis = imputed.tau.hat.logis))
}

tau.bar.var.random <- function(X, observed.time, tau.bar.b) {
  n <- length(observed.time)
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  observed.time.1 <- observed.time[X == 1]
  
  ## limiting variance estimate using estimated F1, F0 and estimated p1, p0
  F1.est <- ecdf(observed.time.1)(observed.time)
  F0.est <- ecdf(observed.time.0)(observed.time)
  p1.est <- mean(X)
  p0.est <- 1 - p1.est
  
  V0 <- (1 - 2 * F1.est) * (X == 0)
  V1 <- (2 * F0.est - 1) * (X == 1)
  
  sigma1.square.est <- mean((p1.est * V0 + p0.est * V1) ^ 2) - (2 * p0.est * p1.est * tau.bar.b) ^ 2
  var.bar.tau.b.est <- sigma1.square.est / (n * p1.est ^ 2 * p0.est ^ 2)
  var.bar.tau.b.est <- var.bar.tau.b.est - tau.bar.b ^ 2 * (p1.est - p0.est) ^ 2 / (n * p1.est * p0.est)
  
  return(var.bar.tau.b.est)
}

tau.bar.var.fixed <- function(X, observed.time, tau.bar.b) {
  n <- length(observed.time)
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  observed.time.1 <- observed.time[X == 1]
  
  ## limiting variance estimate using estimated F1, F0 and estimated p1, p0
  F1.est <- ecdf(observed.time.1)(observed.time.0)
  F0.est <- ecdf(observed.time.0)(observed.time.1)
  p1.est <- mean(X)
  p0.est <- 1 - p1.est
  
  sigma01.square.est <- mean((1 - 2 * F1.est) ^ 2) - tau.bar.b ^ 2
  sigma10.square.est <- mean((2 * F0.est - 1) ^ 2) - tau.bar.b ^ 2
  
  var.bar.tau.b.est <- (sigma01.square.est / p0.est + sigma10.square.est / p1.est) / (N1 + N0)
  
  return(var.bar.tau.b.est)
}

tau.hat.var1.random <- function(X, observed.time, delta, tau.hat) {
  ## var1 estimation
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  var1 <- tau.hat ^ 2 * (p1.hat - p0.hat) ^ 2 / (p0.hat * p1.hat)
  
  return(var1)
}

tau.hat.var2.random <- function(X, observed.time, delta) {
  n <- length(observed.time)
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  delta.0 <- delta[X == 0]
  observed.time.1 <- observed.time[X == 1]
  delta.1 <- delta[X == 1]
  
  ## estiamte the survival functions of censoring variable
  G0.fit <- survfit(Surv(observed.time.0, 1 - delta.0) ~ 1)
  G1.fit <- survfit(Surv(observed.time.1, 1 - delta.1) ~ 1)
  G0.fit_func <- stepfun(G0.fit$time, c(1, G0.fit$surv))
  G1.fit_func <- stepfun(G1.fit$time, c(1, G1.fit$surv))
  
  ## all combinations
  obs.time.mat.0 <- matrix(observed.time.0, nrow = N0, ncol = N1)
  obs.time.mat.1 <- matrix(observed.time.1, nrow = N0, ncol = N1, byrow = TRUE)
  delta.mat.0 <- matrix(as.logical(delta.0), nrow = N0, ncol = N1)
  delta.mat.1 <- matrix(as.logical(delta.1), nrow = N0, ncol = N1, byrow = TRUE)
  
  ## tilde{Y}
  min.index <- obs.time.mat.0 <= obs.time.mat.1
  obs.min.time.mat <- ifelse(min.index, obs.time.mat.0, obs.time.mat.1)
  
  ## orderable indicator
  orderable.indicator <- ifelse(min.index, delta.mat.0, delta.mat.1)
  orderable.indicator <- ifelse(obs.time.mat.0 == obs.time.mat.1 & delta.mat.0 == delta.mat.1,
                                FALSE,
                                orderable.indicator)
  
  ## G0.hat and G1.hat at tilde{Y}
  G0.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G0.fit_func)
  G1.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G1.fit_func)
  
  ## concordance and discordance
  conc <- ifelse(min.index, 1, -1)
  
  ## var2 estimation
  kappa_func <- function(u) {
    trunc.indicator <- (obs.min.time.mat >= u)
    
    kappa.val <- ifelse(orderable.indicator & trunc.indicator,
                        conc / G0.val / G1.val, 0)
    
    return(sum(kappa.val) / choose(n, 2))
  }
  
  kappa.val <- sapply(X = observed.time, FUN = kappa_func)
  R0.val <- apply(sapply(X = observed.time, FUN = "<=", observed.time.0), 2, sum)
  R1.val <- apply(sapply(X = observed.time, FUN = "<=", observed.time.1), 2, sum)
  var2.0 <- n * sum(ifelse((1 - delta) * (1 - X), (kappa.val / R0.val) ^ 2, 0))
  var2.1 <- n * sum(ifelse((1 - delta) * (X), (kappa.val / R1.val) ^ 2, 0))
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  var2 <- (var2.0 + var2.1) / (4 * p0.hat ^ 2 * p1.hat ^ 2)
  
  return(var2)
}

tau.hat.var3.random <- function(X, observed.time, delta, tau.hat) {
  n <- length(observed.time)
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  delta.0 <- delta[X == 0]
  observed.time.1 <- observed.time[X == 1]
  delta.1 <- delta[X == 1]
  
  ## estiamte the survival functions of censoring variable
  G0.fit <- survfit(Surv(observed.time.0, 1 - delta.0) ~ 1)
  G1.fit <- survfit(Surv(observed.time.1, 1 - delta.1) ~ 1)
  G0.fit_func <- stepfun(G0.fit$time, c(1, G0.fit$surv))
  G1.fit_func <- stepfun(G1.fit$time, c(1, G1.fit$surv))
  
  ## var3 estimation
  var3 <- 0
  for(k in 1:length(X)) {
    if(X[k] == 0) {
      min.index <- sapply(observed.time.1, "<=", observed.time[k])
      obs.min.time <- ifelse(min.index, observed.time.1, observed.time[k])
      orderable.indicator <- ifelse(min.index, delta.1, delta[k])
      orderable.indicator <- ifelse(observed.time.1 == observed.time[k] & delta.1 == delta[k],
                                    FALSE,
                                    orderable.indicator)
      G0.val <- sapply(obs.min.time, G0.fit_func)
      G1.val <- sapply(obs.min.time, G1.fit_func)
      conc <- ifelse(min.index, -1, 1)
      U <- ifelse(orderable.indicator, conc / G0.val / G1.val, 0)
      expt <- sum(U) / n
      var3 <- var3 + expt ^ 2
    }
    
    if(X[k] == 1) {
      min.index <- sapply(observed.time.0, "<=", observed.time[k])
      obs.min.time <- ifelse(min.index, observed.time.0, observed.time[k])
      orderable.indicator <- ifelse(min.index, delta.0, delta[k])
      orderable.indicator <- ifelse(observed.time.0 == observed.time[k] & delta.0 == delta[k],
                                    FALSE,
                                    orderable.indicator)
      G0.val <- sapply(obs.min.time, G0.fit_func)
      G1.val <- sapply(obs.min.time, G1.fit_func)
      conc <- ifelse(min.index, 1, -1)
      U <- ifelse(orderable.indicator, conc / G0.val / G1.val, 0)
      expt <- sum(U) / n
      var3 <- var3 + expt ^ 2
    }
  }
  
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  theta1.square <- var3 / n - (2 * p0.hat * p1.hat * tau.hat) ^ 2
  var3 <- theta1.square / (p0.hat * p1.hat) ^ 2
  
  return(var3)
}

tau.hat.var2.fixed <- function(X, observed.time, delta) {
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  delta.0 <- delta[X == 0]
  observed.time.1 <- observed.time[X == 1]
  delta.1 <- delta[X == 1]
  
  ## estiamte the survival functions of censoring variable
  G0.fit <- survfit(Surv(observed.time.0, 1 - delta.0) ~ 1)
  G1.fit <- survfit(Surv(observed.time.1, 1 - delta.1) ~ 1)
  G0.fit_func <- stepfun(G0.fit$time, c(1, G0.fit$surv))
  G1.fit_func <- stepfun(G1.fit$time, c(1, G1.fit$surv))
  
  ## all combinations
  obs.time.mat.0 <- matrix(observed.time.0, nrow = N0, ncol = N1)
  obs.time.mat.1 <- matrix(observed.time.1, nrow = N0, ncol = N1, byrow = TRUE)
  delta.mat.0 <- matrix(as.logical(delta.0), nrow = N0, ncol = N1)
  delta.mat.1 <- matrix(as.logical(delta.1), nrow = N0, ncol = N1, byrow = TRUE)
  
  ## tilde{Y}
  min.index <- obs.time.mat.0 <= obs.time.mat.1
  obs.min.time.mat <- ifelse(min.index, obs.time.mat.0, obs.time.mat.1)
  
  ## orderable indicator
  orderable.indicator <- ifelse(min.index, delta.mat.0, delta.mat.1)
  orderable.indicator <- ifelse(obs.time.mat.0 == obs.time.mat.1 & delta.mat.0 == delta.mat.1,
                                FALSE,
                                orderable.indicator)
  
  ## G0.hat and G1.hat at tilde{Y}
  G0.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G0.fit_func)
  G1.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G1.fit_func)
  
  ## concordance and discordance
  conc <- ifelse(min.index, 1, -1)
  
  ## var1 estimation
  eta_func <- function(u) {
    trunc.indicator <- (obs.min.time.mat >= u)
    
    eta.val <- ifelse(orderable.indicator & trunc.indicator,
                      conc / G0.val / G1.val, 0)
    
    return(sum(eta.val) / N0 / N1)
  }
  
  eta.val <- sapply(X = observed.time, FUN = eta_func)
  R0.val <- apply(sapply(X = observed.time, FUN = "<=", observed.time.0), 2, sum)
  R1.val <- apply(sapply(X = observed.time, FUN = "<=", observed.time.1), 2, sum)
  var2.0 <- N0 * sum(ifelse((1 - delta) * (1 - X), (eta.val / R0.val) ^ 2, 0))
  var2.1 <- N1 * sum(ifelse((1 - delta) * (X), (eta.val / R1.val) ^ 2, 0))
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  var2 <- var2.0 / p0.hat + var2.1 / p1.hat
  
  return(var2)
}

tau.hat.var3.fixed <- function(X, observed.time, delta, tau.hat) {
  n <- length(X)
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  delta.0 <- delta[X == 0]
  observed.time.1 <- observed.time[X == 1]
  delta.1 <- delta[X == 1]
  
  ## estiamte the survival functions of censoring variable
  G0.fit <- survfit(Surv(observed.time.0, 1 - delta.0) ~ 1)
  G1.fit <- survfit(Surv(observed.time.1, 1 - delta.1) ~ 1)
  G0.fit_func <- stepfun(G0.fit$time, c(1, G0.fit$surv))
  G1.fit_func <- stepfun(G1.fit$time, c(1, G1.fit$surv))
  
  ## var3 estimation
  var3.01 <- 0
  for(j in 1:length(observed.time.1)) {
    min.index <- sapply(observed.time.0, "<=", observed.time.1[j])
    obs.min.time <- ifelse(min.index, observed.time.0, observed.time.1[j])
    orderable.indicator <- ifelse(min.index, delta.0, delta.1[j])
    orderable.indicator <- ifelse(observed.time.0 == observed.time.1[j] & delta.0 == delta.1[j],
                                  FALSE,
                                  orderable.indicator)
    G0.val <- sapply(obs.min.time, G0.fit_func)
    G1.val <- sapply(obs.min.time, G1.fit_func)
    conc <- ifelse(min.index, 1, -1)
    U <- ifelse(orderable.indicator, conc / G0.val / G1.val, 0)
    expt <- sum(U) / N0
    var3.01 <- var3.01 + expt ^ 2
  }
  sigma01.square <- var3.01 / N1 - tau.hat ^ 2
  
  var3.10 <- 0
  for(i in 1:length(observed.time.0)) {
    min.index <- sapply(observed.time.1, "<=", observed.time.0[i])
    obs.min.time <- ifelse(min.index, observed.time.1, observed.time.0[i])
    orderable.indicator <- ifelse(min.index, delta.1, delta.0[i])
    orderable.indicator <- ifelse(observed.time.1 == observed.time.0[i] & delta.1 == delta.0[i],
                                  FALSE,
                                  orderable.indicator)
    G0.val <- sapply(obs.min.time, G0.fit_func)
    G1.val <- sapply(obs.min.time, G1.fit_func)
    conc <- ifelse(min.index, -1, 1)
    U <- ifelse(orderable.indicator, conc / G0.val / G1.val, 0)
    expt <- sum(U) / N1
    var3.10 <- var3.10 + expt ^ 2
  }
  sigma10.square <- var3.10 / N0 - tau.hat ^ 2
  
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  var3 <- sigma01.square / p1.hat + sigma10.square / p0.hat
  
  return(var3)
}

restricted.tau.hat_func <- function(X, observed.time, delta, t.star) {
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  delta.0 <- delta[X == 0]
  observed.time.1 <- observed.time[X == 1]
  delta.1 <- delta[X == 1]
  
  ## estiamte the survival functions of censoring variable
  G0.fit <- survfit(Surv(observed.time.0, 1 - delta.0) ~ 1)
  G1.fit <- survfit(Surv(observed.time.1, 1 - delta.1) ~ 1)
  G0.fit_func <- stepfun(G0.fit$time, c(1, G0.fit$surv))
  G1.fit_func <- stepfun(G1.fit$time, c(1, G1.fit$surv))
  
  ## all combinations
  obs.time.mat.0 <- matrix(observed.time.0, nrow = N0, ncol = N1)
  obs.time.mat.1 <- matrix(observed.time.1, nrow = N0, ncol = N1, byrow = TRUE)
  delta.mat.0 <- matrix(as.logical(delta.0), nrow = N0, ncol = N1)
  delta.mat.1 <- matrix(as.logical(delta.1), nrow = N0, ncol = N1, byrow = TRUE)
  
  ## tilde{Y}
  min.index <- obs.time.mat.0 <= obs.time.mat.1
  obs.min.time.mat <- ifelse(min.index, obs.time.mat.0, obs.time.mat.1)
  
  ## orderable indicator
  orderable.indicator <- ifelse(min.index, delta.mat.0, delta.mat.1)
  orderable.indicator <- ifelse(obs.time.mat.0 == obs.time.mat.1 & delta.mat.0 == delta.mat.1,
                                FALSE,
                                orderable.indicator)
  
  ## restricted region indicator
  restricted.indicator <- ifelse(obs.min.time.mat <= t.star, 1, 0)
  
  ## G0.hat and G1.hat at tilde{Y}
  G0.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G0.fit_func)
  G1.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G1.fit_func)
  
  ## concordance and discordance
  conc <- ifelse(min.index, 1, -1)
  
  ## stat
  U <- ifelse(orderable.indicator * restricted.indicator, conc / G0.val / G1.val, 0)
  tau.hat <- sum(U) / N0 / N1
  
  return(list(tau.hat = tau.hat,
              U = sum(U)))
}