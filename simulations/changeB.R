library(dplyr)
library(readr)

### Setup
B = 10 # degree of covariate shift between source & target (S_B)
R = 2 # degree of covariate shift between treatment & control (S_R)
c = 1 # complexity of m_a(x) relative to tau(x)
d = 4 # dimension of X
beta = 1 # number of covariates having covariate shift
expit = function(x) 1/(1 + exp(-x))
ps = function(x) expit(sum(x[1:4])/8*R) # propensity score
or0 = function(x) c*(2 * (abs(x[1]) - pi/4) * (abs(x[1]) >= pi/2) + abs(x[1]) * (abs(x[1]) < pi/2)) -
  0.5 * ((sin(x[1]))) # outcome regression f_0
or1 = function(x) c*(2 * (abs(x[1]) - pi/4) * (abs(x[1]) >= pi/2) + abs(x[1]) * (abs(x[1]) < pi/2)) +
  0.5 * ((sin(x[1]))) # outcome regression f_1
sd = 0.5 # sd of Y|A,X
sim = 1 # simulation times
nnew = 10000 # number of new samples from T for MSE estimation

rho = 5 # parameter rho in matern kernel
matern.kernel_kappa2 = function(u, rho) # matern kernel when kappa = 2
{
  u = u + 1e-100
  out = 4 * exp(-2 * sqrt(2) * u/rho)/((pi^0.5) * rho)
  return(out)
}
Kxx = function(x) matern.kernel_kappa2(u = as.matrix(dist(x)), rho = rho) # K(x,x)
Kxy = function(x, y) matern.kernel_kappa2(u = proxy::dist(x, y), rho = rho) # K(x,y)

# read the functions
source("../separate_regression.R")
source("../coke.R")
source("../dr_cate.R")
source("../acw_cate.R")

mse_coke = mse_dr = mse_acw = mse_sr = rep(NA, sim)
header = c("B", "R", "c", "method", "risk")
results_df = data.frame(matrix(ncol = length(header), nrow = 0))
colnames(results_df) = header

for (B in c(1,5,10,15,20,25)) { # vary S_B
  nt = round((70*B^0.5 + 12*R + 5) * 5) # number of samples in source data T
  ns = nt*4 # number of samples in target data S
  
  p = 1/(1 + B^(1/beta)) # p controls degree of covariate shift between S & T

  ### Simulation
  for (time in 1:sim) {
    ### Data Generation
    # source data S
    S = data.frame(matrix(nrow = ns, ncol = d))
    colnames(S) = paste0("x", 1:d)
    # xi (i=2,3,...) ~ Unif(-pi, pi)
    S[, (beta + 1):d] = matrix(runif((d - beta)*ns, min = -pi, max = pi), nrow = ns, ncol = d - beta)
    # x1 ~ a mixture of 2 uniform distributions, more likely to be negative
    for (j in 1:beta) {
      for (i in 1:ns) {
        repeat {
          x_candi = runif(1, min = -pi, max = pi)
          if (x_candi <= 0) {if (runif(1) > p) {break}}
          else {if (runif(1) > 1 - p) {break}}
        }
        S[i, j] = x_candi
      }
    }
    # generate A & Y in S based on true propensity score & outcome regression models
    S$ps = apply(S[,1:d], 1, ps)
    S$a = mapply(function(prob) rbinom(1, size = 1, prob = prob), S$ps) # A: Bernoulli distribution
    S$or[S$a == 0] = apply(S[S$a == 0, 1:d], 1, or0)
    S$or[S$a == 1] = apply(S[S$a == 1, 1:d], 1, or1)
    S$y = mapply(function(mean) rnorm(1, mean = mean, sd = sd), S$or) # Y: normal distribution
    S = subset(S, select = -c(or, ps))
    
    # target data T
    T = data.frame(matrix(nrow = nt, ncol = d))
    colnames(T) = paste0("x", 1:d)
    T[, (beta + 1):d] = matrix(runif((d - beta)*nt, min = -pi, max = pi), nrow = nt, ncol = d - beta)
    # x1 ~ a mixture of 2 uniform distributions, more likely to be positive
    for (j in 1:beta) {
      for (i in 1:nt) {
        repeat {
          x_candi = runif(1, min = -pi, max = pi)
          if (x_candi <= 0) {if (runif(1) > 1 - p) {break}}
          else {if (runif(1) > p) {break}}
        }
        T[i, j] = x_candi
      }
    }
    
    # new data drawn from target population for MSE estimation
    X_new = matrix(nrow = nnew, ncol = d)
    X_new[, (beta + 1):d] = matrix(runif((d - beta)*nnew, min = -pi, max = pi), nrow = nnew, ncol = d - beta)
    for (j in 1:beta) {
      for (i in 1:nnew) {
        repeat {
          x_candi = runif(1, min = -pi, max = pi)
          if (x_candi <= 0) {if (runif(1) > 1 - p) {break}}
          else {if (runif(1) > p) {break}}
        }
        X_new[i, j] = x_candi
      }
    }
    true = apply(X_new, 1, or1) - apply(X_new, 1, or0) # true CATE on X_new
    
    # Estimate CATE on X_new
    est_sr = separate_regression(S, T, X_new)
    est_coke = coke(S, T, X_new)
    est_dr = dr_cate(S, T, X_new)
    est_acw = acw_cate(S, T, X_new)
    
    # Compute MSE
    mse_sr[time] = sum((est_sr - true)^2)/nnew
    mse_coke[time] = sum((est_coke - true)^2)/nnew
    mse_dr[time] = sum((est_dr - true)^2)/nnew
    mse_acw[time] = sum((est_acw - true)^2)/nnew
    
    results_df = rbind(results_df, 
                       c(B, R, c, "SR", mse_sr[time]),
                       c(B, R, c, "DR", mse_dr[time]),
                       c(B, R, c, "ACW", mse_acw[time]),
                       c(B, R, c, "COKE", mse_coke[time]))
}}

colnames(results_df) = header
save(results_df, file = paste("./", "changeB_seed", seedflag, ".Rdata", sep = ""))
