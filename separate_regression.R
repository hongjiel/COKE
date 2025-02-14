# Separate Regression (Benchmark)
separate_regression = function(S, T, X_new) {
  d = ncol(X_new) # dimension of X
  ns = nrow(S) # number of samples in S
  lambdas = 2^(0:ceiling(log2(10*(ns/2))))/10/(ns/2) # candidate lambda values
  
  # Split S into 2 parts
  split_idx = sample(1:ns, ceiling(ns/2))
  S1 = S[split_idx, ]
  S2 = S[-split_idx, ]
  X1or1 = as.matrix(S1[S1$a == 1, 1:d])
  X1or0 = as.matrix(S1[S1$a == 0, 1:d])
  Y1or1 = S1$y[S1$a == 1]
  Y1or0 = S1$y[S1$a == 0]
  n1or1 = length(Y1or1)
  n1or0 = length(Y1or0)
  X2or1 = as.matrix(S2[S2$a == 1, 1:d])
  X2or0 = as.matrix(S2[S2$a == 0, 1:d])
  Y2or1 = S2$y[S2$a == 1]
  Y2or0 = S2$y[S2$a == 0]
  # Compute kernel matrices
  Kmat1or1 = Kxx(X1or1)
  Kmat1or0 = Kxx(X1or0)
  
  # Choose lambda for f0 by cross–validation
  estmu0_lambdas = sapply(lambdas, function(lambda) {
    ((Y1or0 %*% solve(Kmat1or0 + n1or0 * lambda * diag(n1or0))) %*% Kxy(X1or0, X2or0))[,]
  })
  ssemu0_lambdas = apply(estmu0_lambdas, 2, function(est) {
    sum((est - Y2or0)^2)
  })
  bestlambda_mu0 = lambdas[which.min(ssemu0_lambdas)]
  
  # Choose lambda for f1
  estmu1_lambdas = sapply(lambdas, function(lambda) {
    ((Y1or1 %*% solve(Kmat1or1 + n1or1 * lambda * diag(n1or1))) %*% Kxy(X1or1, X2or1))[,]
  })
  ssemu1_lambdas = apply(estmu1_lambdas, 2, function(est) {
    sum((est - Y2or1)^2)
  })
  bestlambda_mu1 = lambdas[which.min(ssemu1_lambdas)]
  
  # Final estimation on X_new
  est_mu1_new = ((Y1or1 %*% solve(Kmat1or1 + n1or1 * bestlambda_mu1 * diag(n1or1))) %*% Kxy(X1or1, X_new))[,]
  est_mu0_new = ((Y1or0 %*% solve(Kmat1or0 + n1or0 * bestlambda_mu0 * diag(n1or0))) %*% Kxy(X1or0, X_new))[,]
  est_sr = est_mu1_new - est_mu0_new
  
  return(est_sr)
}