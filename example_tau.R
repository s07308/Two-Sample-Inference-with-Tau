#### Example: Water ceontent in Soil
X1 <- c(15.1, 11.2, 10.3, 10.8, 16.6, 8.3, 9.1, 12.3, 9.1,
        14.3, 10.7, 16.1, 10.2, 15.2, 8.9, 9.5, 9.6, 11.3,
        14, 11.3, 15.6, 11.2, 13.8, 9, 8.4, 8.2, 12, 13.9,
        11.6, 16, 9.6, 11.4, 8.4, 8, 14.1, 10.9, 13.2, 13.8,
        14.6, 10.2, 11.5, 13.1, 14.7, 12.5, 10.2, 11.8, 11,
        12.7, 10.3, 10.8, 11, 12.6, 10.8, 9.6, 11.5, 10.6,
        11.7, 10.1, 9.7, 9.7, 11.2, 9.8, 10.3, 11.9, 9.7,
        11.3, 10.4, 12, 11, 10.7, 8.8, 11.1)

X2 <- c(12.1, 10.2, 13.6, 8.1, 13.5, 7.8, 11.8, 7.7, 8.1,
        9.2, 14.1, 8.9, 13.9, 7.5, 12.6, 7.3, 14.9, 12.2,
        7.6, 8.9, 13.9, 8.4, 13.4, 7.1, 12.4, 7.6, 9.9,
        26, 7.3, 7.4, 14.3, 8.4, 13.2, 7.3, 11.3, 7.5,
        9.7, 12.3, 6.9, 7.6, 13.8, 7.5, 13.3, 8, 11.3, 6.8,
        7.4, 11.7, 11.8, 7.7, 12.6, 7.7, 13.2, 13.9, 10.4,
        12.8, 7.6, 10.7, 10.7, 10.9, 12.5, 11.3, 10.7, 13.2,
        8.9, 12.9, 7.7, 9.7, 9.7, 11.4, 11.9, 13.4, 9.2,
        13.4, 8.8, 11.9, 7.1, 8.5, 14, 14.2)

tau_fixed(X1, X2)
tau_random(X1, X2)

#### Example: Time to infection for patients receiving Kidney Dialysis
library(survival)
library(boot)

KD <- read.table("~/Work/TwoSample_Inference_with_Tau/data/kidney_dialysis.txt", header = TRUE)
KD$treatment <- ifelse(KD$treatment == 1, 1, 0)

ipcw.tau.boot <- function(data, indices) {
        d <- data[indices, ]
        tau.hat <- tau_ipcw(d$time, d$delta, d$treatment)
        
        return(tau.hat)
}

## bootstrap under random grouping 
boot.random <- boot(data = KD, statistic = ipcw.tau.boot, R = 2000)
boot.random$t0
quantile(boot.random$t, c(0.025, 0.975))

## bootstrap under fixed grouping 
boot.fixed <- boot(data = KD, statistic = ipcw.tau.boot, strata = KD$treatment, R = 2000)
boot.fixed$t0
quantile(boot.fixed$t, c(0.025, 0.975))
