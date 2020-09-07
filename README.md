# Semi-paired-Association-Test

This repository contains Matlab code to reproduce the experiments in our submission "Semi-paired Association Testing".

# Main routines

### function [p_val0, p_val, p_valSemi, Sta, StaSemi] = HSIC_Test(Kx, Ky, Kzx, Kzy, pars)
- SAT-rx: Hilbert Schmidt Independence Criterion (HSIC) test using semi-paired data
 Inputs: 
- Kx - kernel matrix on x (NxN, the first np x np block contains paired data)
- Ky - kernel matrix on y (MxM, the first np x np block contains paired data)
- Kzx - kernel matrix on covariate z (for x, empty if not exist)
- Kzy - kernel matrix on covariate z (for y, empty if not exist)
- pars - hyperparameters
Outputs:
- p_val0 - p value of the original HSIC using only paired data
- p_val - p value of our Semi-paired test (SAT), only improve null distribution
- p_valSemi - p value of our SAT, improve both test statistics and null distribution
- Sta - test statistic of the original HSIC using only paired data
- StaSemi - test statistic of our SAT

### function [p_val0, p_val, p_valSemi, Sta, StaSemi] = LMM_ScoreTest(Kx, Ky, Kzx, Kzy, pars)
- SAT-fx: Variance Component Score Test (VCST) using semi-paired data. X-phenotype, Y-genotype
Inputs: 
- Kx - kernel matrix on x (NxN, the first np x np block contains paired data)
- Ky - kernel matrix on y (MxM, the first np x np block contains paired data)
- Kzx - kernel matrix on covariate z (for x, empty if not exist)
- Kzy - kernel matrix on covariate z (for y, empty if not exist)
- pars - hyperparameters
Outputs:
- p_val0 - p value of the original VCST using only paired data
- p_val - p value of our Semi-paired test (SAT), only improve null distribution
- p_valSemi - p value of our SAT, improve both test statistics and null distribution
- Sta - test statistic of the original VCST using only paired data
- StaSemi - test statistic of our SAT

# Experiments scripts
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
Our code depends on [KCI-test](http://people.tuebingen.mpg.de/kzhang/KCI-test.zip), [MEGHA-v1](https://scholar.harvard.edu/tge/software/megha).

