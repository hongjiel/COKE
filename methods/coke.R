# COKE (Proposed Method)
coke = function(S, T, X_new) {
  d = ncol(X_new) # dimension of X
  ns = nrow(S) # number of samples in S
  lambdas = 2^(0:ceiling(log2(10*(ns/2))))/10/(ns/2) # candidate lambda values
  
  # split S into 2 parts for the cross-fitting methods
  indices = sample(1:ns)
  n1 = ceiling(ns/2)
  D1 = S[indices[1:n1], ]
  D2 = S[indices[(n1 + 1):ns], ]
  S_split = list(D1, D2)
  
  # define the cyclic permutations
  perms = list(c(1, 2), c(2, 1))
  est_list = list()
  XT = as.matrix(T)
  
  for (i in seq_along(perms)) {
    perm = perms[[i]]
    
    # nuisance parameter estimation on S1 (using small lambdas)
    S1 = S_split[[perm[1]]]
    X1_or1 = as.matrix(S1[S1$a == 1,1:d])
    X1_or0 = as.matrix(S1[S1$a == 0,1:d])
    
    # create pseudo-outcome on S1
    X1 = as.matrix(S1[,1:d])
    S1$mu0hat = (S1$y[S1$a == 0] %*% solve(Kxx(X1_or0) + nrow(X1_or0) * min(lambdas) * diag(nrow(X1_or0))) %*% Kxy(X1_or0, X1))[,]
    S1$mu1hat = (S1$y[S1$a == 1] %*% solve(Kxx(X1_or1) + nrow(X1_or1) * min(lambdas) * diag(nrow(X1_or1))) %*% Kxy(X1_or1, X1))[,]
    S1 = S1 %>% 
      mutate(phihat = (a == 1) * (y - mu0hat) + (a == 0) * (mu1hat - y))
    
    # target parameter estimation using:
      ##S1 (candidate models)
      ##S2 (imputation model using a small lambda)
      ##T (pseudo label)
    S2 = S_split[[perm[2]]]
    Kmat1 = Kxx(X1)
    X2_or1 = as.matrix(S2[S2$a == 1,1:d])
    X2_or0 = as.matrix(S2[S2$a == 0,1:d])
    
    # create pseudo label on XT
    pseudo = (S2$y[S2$a == 1] %*% solve(Kxx(X2_or1) + nrow(X2_or1) * min(lambdas) * diag(nrow(X2_or1))) %*% Kxy(X2_or1, XT))[,] -
      (S2$y[S2$a == 0] %*% solve(Kxx(X2_or0) + nrow(X2_or0) * min(lambdas) * diag(nrow(X2_or0))) %*% Kxy(X2_or0, XT))[,]
    esttg_lambdas = sapply(lambdas, function(lambda) {
      ((S1$phihat %*% solve(Kmat1 + nrow(Kmat1) * lambda * diag(nrow(Kmat1)))) %*% Kxy(X1, XT))[,]
    })
    ssetg_lambdas = apply(esttg_lambdas, 2, function(est) {
      sum((est - pseudo)^2)
    })
    bestlambda_tg = lambdas[which.min(ssetg_lambdas)]
    
    # final estimation on X_new
    est_list[[i]] = (S1$phihat %*% solve(Kmat1 + nrow(Kmat1) * bestlambda_tg * diag(nrow(Kmat1))) %*% Kxy(X1, X_new))[,]
  }
  
  # return the average estimate over the three permutations
  est_final = Reduce("+", est_list) / length(est_list)
  return(est_final)
}