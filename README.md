# Two-Sample Inference with Tau

The R code in this repository is used to implement the inference procedures based on Kendall’s tau (&tau;<sub>b</sub>) between a binary group indicator and a continuous variable which may be subject to right-censoring. The methods are proposed by Yi-Cheng Tai, Weijing Wang and Martin T. Wells and will be submitted for publication. <br>

### tau_fixed()
When the group indicators are known values as in randomized clinical trials, testing results for H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> <br> and H<sub>0</sub>: &tau;<sub>b</sub> = 0 and confidence intervals of &tau;<sub>b</sub> are given. <br>

#### Arguments
`x`: a non-empty numeric vector of data <br>
`y`: a non-empty numeric vector of data <br>

#### Value
A list containing the following components <br>
`tau`: the estimated value of &tau;<sub>b</sub> <br>
`var.tau.fixed`: the variance of the estimator of &tau; under H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> <br>
`var.tau.0`: the variance of the estimator of &tau; under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`var.tau.general`: the variance of the estimator of &tau; in general <br>
`z.score.fixed`: the z-score under H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> <br>
`z.score.fixed`: the z-score under H<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`p.value.fixed`: p-value under H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> <br>
`p.value.0`: p-value underH<sub>0</sub>: &tau;<sub>b</sub> = 0 <br>
`ci`: the 95% confidence interval of &tau;<sub>b</sub> <br>

### tau_random()
When the group indicators are known values as in observational surveys, testing results for H<sub>0</sub>: F<sub>x</sub> = F<sub>y</sub> and H<sub>0</sub>: &tau;<sub>b</sub> = 0 and confidence intervals of &tau;<sub>b</sub> are given. <br> 

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
An IPCW estimator of &tau;<sub>b</sub> is given when the observations are subject to right-censoring. Users can perform bootstrap procedures to obtain variance estimate for further inference problems. For fixed group indicators, bootstrap samples are drawn separately from each group. For random group indicators, bootrstap samples are drawn from the combined dataset. <br>  

#### Arguments
`Y1`: follow-up time <br>
`delta`: the status indicator. Typically, 0: censored, 1: died <br>
`X1`: group indicator, coded as 0 and 1 <br>

#### Value
the estimaed value of &tau;<sub>b</sub> <br>

## Example
#### Complete Data
The soil water contents (% water by volume) collected from two experimental fields growing bell peppers are under comparison (Gumpertz et al., 1997).

```
tau_fixed(X1, X2)
$tau
[1] 0.1901042

$var.tau.fixed
[1] 0.008854167

$var.tau.0
[1] 0.008928687

$var.tau.general
[1] 0.008922413

$z.score.fixed
[1] 2.020309

$z.score.0
[1] 2.01186

$p.value.fixed
[1] 0.02167568

$p.value.0
[1] 0.02211733

$ci
[1] 0.004968861 0.375239472
```

#### Censored Data


