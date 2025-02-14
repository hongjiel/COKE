# COKE (Proposed Method)
coke = function(S, T, X_new) {
  d = ncol(X_new) # dimension of X
  ns = nrow(S) # number of samples in S
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
  return(est_final)
}