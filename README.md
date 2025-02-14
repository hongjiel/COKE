This repository contains R functions and simulation scripts used in the paper "Transfer Learning of CATE with Kernel Ridge Regression".

## Methods (R Functions)

- `separate_regression.R`: Separate Regression (`SR`).
- `coke.R`: `COKE`.
- `dr_cate.R`:  DR-Learner for CATE (`DR-CATE`).
- `acw_cate.R`: ACW estimator tailored for CATE estimation (`ACW-CATE`).

## Simulation Scripts

- `changeB.R`: Vary $S_B$ with other parameters fixed under $q=1$.
- `changeR.R`: Vary $S_R$ with other parameters fixed.
- `changeC.R`: Vary $c$ with other parameters fixed.
- `changeN.R`: Vary $n_\mathcal{T} = n/4$ with other parameters fixed.
- `changeB_2dim.R`: Vary $S_B$ with other parameters fixed under $q=2$.
- `changeB_CF.R`: Compare the cross-fitting version of `COKE` with the original Algorithm 3.
