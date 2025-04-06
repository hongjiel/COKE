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

mse_all = mse_perm1 = mse_perm2 = mse_perm3 = rep(NA, sim)
header = c("B", "type", "risk")
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
    
    ### Apply COKE
    lambdas = 2^(0:ceiling(log2(10*(ns/2))))/10/(ns/2) # candidate lambda values
    
    # split S into 3 parts for the cross-fitting methods
    indices = sample(1:ns)
    n1 = ceiling(ns / 3)
    D1 = S[indices[1:n1], ]
    D2 = S[indices[(n1 + 1):(2 * n1)], ]
    D3 = S[indices[(2 * n1 + 1):ns], ]
    S_split = list(D1, D2, D3)
    
    # define the three cyclic permutations
    perms = list(c(1, 2, 3), c(3, 1, 2), c(2, 3, 1))
    est_list = list()
    XT = as.matrix(T)
    
    for (i in seq_along(perms)) {
      perm = perms[[i]]
      
      # nuisance parameter estimation on S_nc (using small lambdas)
      S_nc = S_split[[perm[1]]]
      X_ncor1 = as.matrix(S_nc[S_nc$a == 1,1:d])
      X_ncor0 = as.matrix(S_nc[S_nc$a == 0,1:d])
      
      # create pseudo-outcome on S_tg1
      S_tg1 = S_split[[perm[2]]]
      X_tg1 = as.matrix(S_tg1[,1:d])
      S_tg1$mu0hat = (S_nc$y[S_nc$a == 0] %*% solve(Kxx(X_ncor0) + nrow(X_ncor0) * min(lambdas) * diag(nrow(X_ncor0))) %*% Kxy(X_ncor0, X_tg1))[,]
      S_tg1$mu1hat = (S_nc$y[S_nc$a == 1] %*% solve(Kxx(X_ncor1) + nrow(X_ncor1) * min(lambdas) * diag(nrow(X_ncor1))) %*% Kxy(X_ncor1, X_tg1))[,]
      S_tg1 = S_tg1 %>% 
        mutate(phihat = (a == 1) * (y - mu0hat) + (a == 0) * (mu1hat - y))
      
      # target parameter estimation using:
      ##S_tg1 (candidate models)
      ##S_tg2 (imputation model using a small lambda)
      ##T (pseudo label)
      S_tg2 = S_split[[perm[3]]]
      Kmat_tg1 = Kxx(X_tg1)
      X_tg2or1 = as.matrix(S_tg2[S_tg2$a == 1,1:d])
      X_tg2or0 = as.matrix(S_tg2[S_tg2$a == 0,1:d])
      
      # create pseudo label on XT
      pseudo = (S_tg2$y[S_tg2$a == 1] %*% solve(Kxx(X_tg2or1) + nrow(X_tg2or1) * min(lambdas) * diag(nrow(X_tg2or1))) %*% Kxy(X_tg2or1, XT))[,] -
        (S_tg2$y[S_tg2$a == 0] %*% solve(Kxx(X_tg2or0) + nrow(X_tg2or0) * min(lambdas) * diag(nrow(X_tg2or0))) %*% Kxy(X_tg2or0, XT))[,]
      esttg_lambdas = sapply(lambdas, function(lambda) {
        ((S_tg1$phihat %*% solve(Kmat_tg1 + nrow(Kmat_tg1) * lambda * diag(nrow(Kmat_tg1)))) %*% Kxy(X_tg1, XT))[,]
      })
      ssetg_lambdas = apply(esttg_lambdas, 2, function(est) {
        sum((est - pseudo)^2)
      })
      bestlambda_tg = lambdas[which.min(ssetg_lambdas)]
      
      # final estimation on X_new
      est_list[[i]] = (S_tg1$phihat %*% solve(Kmat_tg1 + nrow(Kmat_tg1) * bestlambda_tg * diag(nrow(Kmat_tg1))) %*% Kxy(X_tg1, X_new))[,]
    }
    
    # return the average estimate over the three permutations
    est_final = Reduce("+", est_list) / length(est_list)
    
    # Compute MSE
    mse_perm1[time] = sum((est_list[[1]] - true)^2)/nnew # no CF (permutation 1)
    mse_perm2[time] = sum((est_list[[2]] - true)^2)/nnew # no CF (permutation 2)
    mse_perm3[time] = sum((est_list[[3]] - true)^2)/nnew # no CF (permutation 3)
    mse_all[time] = sum((est_final - true)^2)/nnew # with CF
    
    results_df = rbind(results_df, 
                       c(B, "perm1", mse_perm1[time]),
                       c(B, "perm2", mse_perm2[time]),
                       c(B, "perm3", mse_perm3[time]),
                       c(B, "all", mse_all[time]))
}}

colnames(results_df) = header
save(results_df, file = paste("./", "changeB_CF_seed", seedflag, ".Rdata", sep = ""))
