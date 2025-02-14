# ACW-CATE (Benchmark)
acw_cate = function(S, T, X_new) {
  d = ncol(X_new) # dimension of X
  ns = nrow(S) # number of samples in S
  lambdas = 2^(0:ceiling(log2(10*(ns/2))))/10/(ns/2) # candidate lambda values
  
  # Split S into 3 parts for the cross-fitting methods
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
  nt = nrow(T)
  
  for (i in seq_along(perms)) {
    perm = perms[[i]]
    S_nc = S_split[[perm[1]]]
    
    # nuisance parameter estimation
    S_nc = S_split[[perm[1]]]
    partS_nc = sample(1:nrow(S_nc), ceiling(nrow(S_nc)/2))
    S_nc1 = S_nc[partS_nc,]
    S_nc2 = S_nc[-partS_nc,]
    X_nc1 = as.matrix(S_nc1[,1:d])
    X_nc1or1 = as.matrix(S_nc1[S_nc1$a == 1,1:d])
    X_nc2or1 = as.matrix(S_nc2[S_nc2$a == 1,1:d])
    X_nc1or0 = as.matrix(S_nc1[S_nc1$a == 0,1:d])
    X_nc2or0 = as.matrix(S_nc2[S_nc2$a == 0,1:d])
    Y_nc1or1 = S_nc1$y[S_nc1$a == 1]
    Y_nc2or1 = S_nc2$y[S_nc2$a == 1]
    Y_nc1or0 = S_nc1$y[S_nc1$a == 0]
    Y_nc2or0 = S_nc2$y[S_nc2$a == 0]
    n_nc1or1 = length(Y_nc1or1)
    n_nc1or0 = length(Y_nc1or0)
    X_ncor1 = as.matrix(S_nc[S_nc$a == 1,1:d])
    X_ncor0 = as.matrix(S_nc[S_nc$a == 0,1:d])
    Kmat_ncor1 = Kxx(X_ncor1)
    Kmat_ncor0 = Kxx(X_ncor0)
    
    # propensity score estimation
    variables = paste0("x", 1:d, collapse = " + ")
    formula = as.formula(paste("a ~", variables))
    model_ps = glm(formula, data = S_nc, family = binomial)
    
    # outcome regression estimation
    # f0
    estmu0_lambdas = sapply(lambdas, function(lambda) {
      ((Y_nc1or0 %*% solve(Kxx(X_nc1or0) + n_nc1or0 * lambda * diag(n_nc1or0))) %*% Kxy(X_nc1or0, X_nc2or0))[,]
    })
    ssemu0_lambdas = apply(estmu0_lambdas, 2, function(est) {
      sum((est - Y_nc2or0)^2)
    })
    bestlambda_mu0 = lambdas[which.min(ssemu0_lambdas)]
    # f1
    estmu1_lambdas = sapply(lambdas, function(lambda) {
      ((Y_nc1or1 %*% solve(Kxx(X_nc1or1) + n_nc1or1 * lambda * diag(n_nc1or1))) %*% Kxy(X_nc1or1, X_nc2or1))[,]
    })
    ssemu1_lambdas = apply(estmu1_lambdas, 2, function(est) {
      sum((est - Y_nc2or1)^2)
    })
    bestlambda_mu1 = lambdas[which.min(ssemu1_lambdas)]
    
    # density ratio estimation
    Stemp = S_nc[, 1:d]
    Ttemp = T
    Stemp$k = 0
    Ttemp$k = 1
    all = rbind(Stemp, Ttemp)
    formula = as.formula(paste("k ~", variables))
    model_dens = glm(formula, data = all, family = binomial)
    
    n2 = nrow(S_split[[perm[2]]])
    n3 = nrow(S_split[[perm[3]]])
    S_tg = rbind(S_split[[perm[2]]], S_split[[perm[3]]])
    weight = predict(model_dens, newdata = S_tg, type = "response") /
      (1 - predict(model_dens, newdata = S_tg, type = "response"))*ns/nt
    weight = weight * nrow(S_tg)/sum(weight)
    weight = pmin(weight, 20) # truncation
    weight = weight * nrow(S_tg)/sum(weight)
    
    # pseudo-outcome phihat
    X_tg = as.matrix(S_tg[,1:d])
    S_tg$what = weight
    S_tg$pihat = predict(model_ps, newdata = S_tg, type = "response")
    S_tg$mu0hat = (S_nc$y[S_nc$a == 0] %*% solve(Kmat_ncor0 + nrow(X_ncor0) * bestlambda_mu0 * diag(nrow(X_ncor0))) %*% Kxy(X_ncor0, X_tg))[,]
    S_tg$mu1hat = (S_nc$y[S_nc$a == 1] %*% solve(Kmat_ncor1 + nrow(X_ncor1) * bestlambda_mu1 * diag(nrow(X_ncor1))) %*% Kxy(X_ncor1, X_tg))[,]
    S_tg = S_tg %>% 
      mutate(phihat = what * (nrow(S_tg) + nt) / nrow(S_tg) * (a * (y - mu1hat) / pihat - (1 - a) * (y - mu0hat) / (1 - pihat))) %>% 
      select(-a, -y, -what, -pihat, -mu0hat, -mu1hat)
    S_tg1 = S_tg[1:n2,]
    S_tg2 = S_tg[(n2 + 1):(n2 + n3),]
    
    T_hat = T
    T_hat$mu0hat = (S_nc$y[S_nc$a == 0] %*% solve(Kmat_ncor0 + nrow(X_ncor0) * bestlambda_mu0 * diag(nrow(X_ncor0))) %*% Kxy(X_ncor0, XT))[,]
    T_hat$mu1hat = (S_nc$y[S_nc$a == 1] %*% solve(Kmat_ncor1 + nrow(X_ncor1) * bestlambda_mu1 * diag(nrow(X_ncor1))) %*% Kxy(X_ncor1, XT))[,]
    T_hat = T_hat %>% 
      mutate(phihat = (nrow(S_tg) + nt) / nt * (mu1hat - mu0hat)) %>% 
      select(-mu0hat, -mu1hat)
    partT = sample(1:nt, ceiling(nt/2))
    T1 = T_hat[partT,]
    T2 = T_hat[-partT,]
    
    # target parameter estimation
      ## S_tg1, T1: training data
      ## S_tg2, T2: test data
    mix1 = rbind(S_tg1, T1)
    mix2 = rbind(S_tg2, T2)
    X_mix1 = as.matrix(mix1[,1:d])
    X_mix2 = as.matrix(mix2[,1:d])
    Kmat_mix1 = Kxx(X_mix1)
    
    estmix_lambdas = sapply(lambdas, function(lambda) {
      ((mix1$phihat %*% solve(Kmat_mix1 + nrow(mix1) * lambda * diag(nrow(mix1)))) %*% Kxy(X_mix1, X_mix2))[,]
    })
    ssemix_lambdas = apply(estmix_lambdas, 2, function(est) {
      sum((est - mix2$phihat)^2)
    })
    bestlambda_mix = lambdas[which.min(ssemix_lambdas)]
    
    # final estimation on X_new
    est_perm = (mix1$phihat %*% solve(Kmat_mix1 + nrow(mix1) * bestlambda_mix * diag(nrow(mix1))) %*% Kxy(X_mix1, X_new))[,]
    est_list[[i]] = est_perm
  }
  
  # return the average estimate over the three permutations
  est_final = Reduce("+", est_list) / length(est_list)
  return(est_final)
}