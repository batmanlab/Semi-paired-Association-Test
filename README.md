# Semi-paired-Association-Test

This repository contains Matlab code to reproduce the experiments in our submission "Semi-paired Association Testing".
<p align="center">
  <img width="55%" height="%55" src="fig1.png">
</p>

# Main routines
SAT-rx: semi-paired test in the random X setting. The inputs are kernel matrices for X, Y, and Z (covariates), the outputs are the p-values and test statistics.

``` function [p_val0, p_val, p_valSemi, Sta, StaSemi] = HSIC_Test(Kx, Ky, Kzx, Kzy, pars)```
SAT-fx: semi-paired test in the fixed X setting. The inputs are kernel matrices for X, Y, and Z (covariates), the outputs are the p-values and test statistics.

``` function [p_val0, p_val, p_valSemi, Sta, StaSemi] = LMM_ScoreTest(Kx, Ky, Kzx, Kzy, pars)```

# Simulation Experiments
## Evaluation of type I and type II errors of SAT-rx on simulated data.
- function exp_simul1_typeI_par_func(h2, sigma_n, dim, iter)
- function exp_simul1_typeII_par_func(h2, sigma_n, dim, iter)

## Evaluation of type I and type II errors of SAT-fx on simulated data.
- function exp_simul1_lmm_typeI_par_func(h2, sigma_n, dim, iter)
- function exp_simul1_lmm_typeII_par_func(h2, sigma_n, dim, iter)

## Evaluation of type I and type II errors of SAT-fx on COPD simulated data.
- function exp_copd_simul_typeI_par_func(h2, sig, ydim, iter)
- function exp_copd_simul_typeII_par_func(h2, sig, ydim, iter)

## P-values on Uganda dataset
- function [p_val0All,p_val0,p_val,p_valSemi,Sta,StaSemi,nl,nlSelUp] = exp_exploration_mimic(dataset,pars) 
- function [p_val0,p_val,p_valSemi,Sta,StaSemi,nl,nlSelUp] = exp_exploration_real(dataset,pars)

## P-values on real imaging
- function exp_copd_geneExpTseng_func(pheno_id, gene_id)

## P-values on biomarker data
- function exp_copd_kl_biomarker_func(fev1_cov)

## P-values on gene expression data
- function exp_copd_kl_geneExp_func(fev1_cov)

# Dependencies 
Our code depends on [KCI-test](http://people.tuebingen.mpg.de/kzhang/KCI-test.zip) and [MEGHA-v1](https://scholar.harvard.edu/tge/software/megha).

