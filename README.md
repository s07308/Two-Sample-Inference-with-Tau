# Two-Sample Inference with Tau

The R code in this repository is used to implement the inference procedures based on Kendall’s tau (&tau;<sub>b</sub>) between a binary group indicator and a continuous variable which may be subject to right-censoring. The methods are proposed by Yi-Cheng Tai, Weijing Wang and Martin T. Wells and will be submitted for publication. <br>

### tau.bar_func()
When the observed failure times do not subject to censoring, testing results for H<sub>0</sub>: S<sub>0</sub> = S<sub>1</sub> <br> and H<sub>0</sub>: &tau;<sub>b</sub> = 0 and confidence intervals of &tau;<sub>b</sub> are given. <br>

#### Arguments
`X`: a non-empty numeric vector of group indicators, encoded as 0 or 1 <br>
`observed.time`: a non-empty numeric vector of data <br>

#### Value
A list containing the following components <br>
`tau.bar`: the estimated value of &tau;<sub>b</sub> <br>
`var.fixed`: the variance of the estimator of &tau;<sub>b</sub> when the group indicators are fixed <br>
`var.random`: the variance of the estimator of &tau;<sub>b</sub> when the group indicators are random <br>
`var.null.0`: the variance of the estimator of &tau;<sub>b</sub> under H<sub>0</sub>: S<sub>0</sub> = S<sub>1</sub> <br>
`var.null.tau`: the variance of the estimator of &tau;<sub>b</sub> under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`ci.fixed`: the 95% confidence interval of &tau;<sub>b</sub> when the group indicators are fixed <br>
`ci.random`: the 95% confidence interval of &tau;<sub>b</sub> when the group indicators are random <br>
`z.val.0`: the z-score under H<sub>0</sub>: S<sub>0</sub> = S<sub>1</sub> <br>
`z.val.tau`: the z-score under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`p.value.0`: p-value under H<sub>0</sub>: S<sub>0</sub> = S<sub>1</sub> <br>
`p.value.tau`: p-value under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>

### tau.hat_func()
An IPCW estimator of &tau;<sub>b</sub> is given when the observations are subject to right-censoring. The testing results for H<sub>0</sub>: S<sub>0</sub> = S<sub>1</sub> <br> and H<sub>0</sub>: &tau;<sub>b</sub> = 0 and confidence intervals of &tau;<sub>b</sub> are given.

#### Arguments
`X`: a non-empty numeric vector of group indicators, encoded as 0 or 1 <br>
`observed.time`: a non-empty numeric vector of follow-up time <br>
`delta`: the status indicator. Typically, 0: censored, 1: died <br>

#### Value
A list containing the following components <br>
`tau.hat`: the estimated value of &tau;<sub>b</sub> <br>
`U`: the sum of scores assigned.
`var.fixed`: the variance of the estimator of &tau;<sub>b</sub> when the group indicators are fixed <br>
`var.random`: the variance of the estimator of &tau;<sub>b</sub> when the group indicators are random <br>
`var.null.tau`: the variance of the estimator of &tau;<sub>b</sub> under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`ci.fixed`: the 95% confidence interval of &tau;<sub>b</sub> when the group indicators are fixed <br>
`ci.random`: the 95% confidence interval of &tau;<sub>b</sub> when the group indicators are random <br>
`z.val.tau`: the z-score under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`p.value.tau`: p-value under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>

### imputed.tau.hat_func()
When the upper bound of the support of censoring variable is shorter than the variable of interest, the estimate of &tau;<sub>b</sub> with imputed association pattern in unidentifiable region is provided. Several parametric distributions are used. <br>

#### Arguments
`X`: a non-empty numeric vector of group indicators, encoded as 0 or 1 <br>
`observed.time`: a non-empty numeric vector of follow-up time <br>
`delta`: the status indicator. Typically, 0: censored, 1: died <br>
`t.star`: a pre-specified value sets the identifiable region <br>

#### Value
A list containing the following components <br>
`weibull`: the estimated value of &tau;<sub>b</sub> with imputed weibull tail <br>
`exp`: the estimated value of &tau;<sub>b</sub> with imputed exponential tail <br>
`lnorm`: the estimated value of &tau;<sub>b</sub> with imputed log-normal tail <br>
`logis`: the estimated value of &tau;<sub>b</sub> with imputed logistic tail <br>

## Example
#### Complete Data
The soil water contents (% water by volume) collected from two experimental fields growing bell peppers are under comparison (Gumpertz et al., 1997).

```
tau.bar_func(X = group, observed.time = obs)
$tau.bar
[1] 0.1901042

$var.fixed
[1] 0.0087689

$var.random
[1] 0.0087689

$var.null.0
[1] 0.008796296

$var.null.tau
[1] 0.009722584

$ci.fixed
[1] 0.006568424 0.373639910

$ci.random
[1] 0.006568424 0.373639910

$z.val.0
[1] 2.026944

$z.val.tau
[1] 1.927972

$p.val.0
[1] 0.04266816

$p.val.tau
[1] 0.05385857
```

#### Censored Data
The dataset is obtained from the book "Survival Analysis Techniques for Censored and Truncated Data" (Klein, John P., Moeschberger, Melvin L., 2003). It recoreded the time to infection for patients receiving Kidney Dialysis. <br>

```
tau.hat_func(X = KD$treatment, observed.time = KD$time, delta = KD$delta)
$tau.hat
[1] -0.4437778

$U
[1] -1450.266

$var.fixed
[1] 0.1051801

$var.random
[1] 0.1237697

$var.null.tau
[1] 0.130941

$ci.fixed
[1] -1.0794231  0.1918676

$ci.random
[1] -1.1333112  0.2457556

$z.val.tau
[1] -1.226388

$p.val.tau
[1] 0.2200529
```

```
imputed.tau.hat_func(X = KD$treatment, observed.time = KD$time, delta = KD$delta, t.star = max(KD$time))
$weibull
[1] -0.6192046

$exp
[1] -0.5082488

$lnorm
[1] -0.6497025

$logis
[1] -0.483848
```

## Remark
The dependency packages include `survival`, `parmsurvfit`.
