# Two-Sample Test with Tau

The R code in this repository contains three functions to do two-sample test in different scenarios. 

### tau_fixed

#### Arguments
x: a non-empty numeric vector vector of data <br>
y: a non-empty numeric vector vector of data

#### Value
A list containing the following components 
tau: the estimated value of tau_b
var.tau.fixed: the variance of tau under H_0: F_x = F_y
var.tau.0: the variance of tau undet H_0: tau_b = 0
var.tau.general: the variance of tau in general
z.score.fixed: the z-score under H_0: F_x = F_y
z.score.fixed: the z_score under H_0: tau_b = 0
p.value.fixed: p-value under H_0: F_x = F_y
p.value.0: p-value under H_0: tau_b = 0
ci: the 95% confidence interval of tau_b


