# Two-Sample Test with Tau

The R code in this repository contains three functions to do two-sample test with &tau; in different scenarios. 

### tau_fixed()
Conduct the two-sample test with &tau; in RCT. <br>

#### Arguments
`x`: a non-empty numeric vector of data <br>
`y`: a non-empty numeric vector of data <br>

#### Value
A list containing the following components <br>
`tau`: the estimated value of &tau;<sub>b</sub> <br>
`var.tau.fixed`: the variance of the estimator of &tau; under H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> <br>
`var.tau.0`: the variance of the estimator of &tau; undet H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`var.tau.general`: the variance of the estimator of &tau; in general <br>
`z.score.fixed`: the z-score under H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> <br>
`z.score.fixed`: the z-score under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`p.value.fixed`: p-value under H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> <br>
`p.value.0`: p-value underH<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`ci`: the 95% confidence interval of &tau;<sub>b</sub> <br>

### tau_random()
Conduct the two-sample test with &tau; in observational study <br>

#### Arguments
`x`: a non-empty numeric vector vector of data <br>
`y`: a non-empty numeric vector vector of data <br>

#### Value
A list containing the following components <br>
`tau`: the estimated value of &tau;<sub>b</sub> <br>
`var.tau.random`: the variance of the estimator of &tau; under H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> <br>
`var.tau`: the variance of the estimator of &tau; in general <br>
`z.score.random`: the z-score under H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> <br>
`z.score.0`: the z-score under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`p.value.random`: p-value under H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> <br>
`p.value.0`: p-value under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`ci`: the 95% confidence interval of &tau;<sub>b</sub> <br>

### tau_ipcw()
Estimate &tau;<sub>b</sub> when the observatiosn are subject to right-censoring. <br>

#### Arguments
`Y1`: follow-up time <br>
`delta`: the status indicator. Typically, 0: censored, 1: died <br>
`X1`: group indicator, coded as 0 and 1 <br>

#### Value
the estimaed value of &tau;<sub>b</sub> <br>

## Example
The soil water contents (% water by volume) collected from two experimental fields growing bell peppers are under comparison (Gumpertz et al., 1997).

```
tau_fixed(X1, X2)
tau_random(X1, X2)
```
