# DR-CATE (Benchmark)
dr_cate = function(S, T, X_new) {
  d = ncol(X_new) # dimension of X
  ns = nrow(S) # number of samples in S
  lambdas = 2^(0:ceiling(log2(10*(ns/2))))/10/(ns/2) # candidate lambda values
  
  # split S into 2 parts for the cross-fitting methods
  indices = sample(1:ns)
  n1 = ceiling(ns / 2)
  D1 = S[indices[1:n1], ]
  D2 = S[indices[(n1 + 1):ns], ]
  S_split = list(D1, D2)
  
  # define the cyclic permutations
  perms = list(c(1, 2), c(2, 1))
  est_list = list()
  variables = paste0("x", 1:d, collapse = " + ")
  
  for (i in seq_along(perms)) {
    perm = perms[[i]]
    
    # nuisance parameter estimation on S_nc
    S_nc = S_split[[perm[1]]]
    partS_nc = sample(1:nrow(S_nc), ceiling(nrow(S_nc)/2))
    S_nc1 = S_nc[partS_nc,]
    S_nc2 = S_nc[-partS_nc,]
    X_nc1 = as.matrix(S_nc1[,1:d])
    X_nc2 = as.matrix(S_nc2[,1:d])
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
    
    # target parameter estimation using:
      ##S_tg1 (training data)
      ##S_tg2 (test data)
    S_tg = S_split[[perm[2]]]
    partS_tg = sample(1:nrow(S_tg), ceiling(nrow(S_tg)/2))
    S_tg1 = S_tg[partS_tg,]
    S_tg2 = S_tg[-partS_tg,]
    X_tg1 = as.matrix(S_tg1[,1:d])
    X_tg2 = as.matrix(S_tg2[,1:d])
    Kmat_tg1 = Kxx(X_tg1)
    K_tg1tg2 = Kxy(X_tg1, X_tg2)
    
    # pseudo-outcome phihat
    S_tg1$pihat = predict(model_ps, newdata = S_tg1, type = "response")
    S_tg1$mu0hat = (S_nc$y[S_nc$a == 0] %*% solve(Kmat_ncor0 + nrow(X_ncor0) * bestlambda_mu0 * diag(nrow(X_ncor0))) %*% Kxy(X_ncor0, X_tg1))[,]
    S_tg1$mu1hat = (S_nc$y[S_nc$a == 1] %*% solve(Kmat_ncor1 + nrow(X_ncor1) * bestlambda_mu1 * diag(nrow(X_ncor1))) %*% Kxy(X_ncor1, X_tg1))[,]
    S_tg1 = S_tg1 %>% 
      mutate(phihat = mu1hat - mu0hat + a * (y - mu1hat) / pihat - (1 - a) * (y - mu0hat) / (1 - pihat))
    
    S_tg2$pihat = predict(model_ps, newdata = S_tg2, type = "response")
    S_tg2$mu0hat = (S_nc$y[S_nc$a == 0] %*% solve(Kmat_ncor0 + nrow(X_ncor0) * bestlambda_mu0 * diag(nrow(X_ncor0))) %*% Kxy(X_ncor0, X_tg2))[,]
    S_tg2$mu1hat = (S_nc$y[S_nc$a == 1] %*% solve(Kmat_ncor1 + nrow(X_ncor1) * bestlambda_mu1 * diag(nrow(X_ncor1))) %*% Kxy(X_ncor1, X_tg2))[,]
    S_tg2 = S_tg2 %>% 
      mutate(phihat = mu1hat - mu0hat + a * (y - mu1hat) / pihat - (1 - a) * (y - mu0hat) / (1 - pihat))
    
    esttg_lambdas = sapply(lambdas, function(lambda) {
      ((S_tg1$phihat %*% solve(Kmat_tg1 + nrow(Kmat_tg1) * lambda * diag(nrow(Kmat_tg1)))) %*% K_tg1tg2)[,]
    })
    ssetg_lambdas = apply(esttg_lambdas, 2, function(est) {
      sum((est - S_tg2$phihat)^2)
    })
    bestlambda_tg = lambdas[which.min(ssetg_lambdas)]
    
    # final estimation on X_new
    est_perm = (S_tg1$phihat %*% solve(Kmat_tg1 + nrow(Kmat_tg1) * bestlambda_tg * diag(nrow(Kmat_tg1))) %*% Kxy(X_tg1, X_new))[,]
    est_list[[i]] = est_perm
  }
  
  # return the average estimate over the three permutations
  est_final = Reduce("+", est_list) / length(est_list)
  return(est_final)
}