# Semi-paired Association Test (SAT)

This repository contains Matlab code to reproduce the experiments in our submission "Semi-paired Association Testing".
<p align="center">
  <img width="40%" height="40%" src="fig1.png">
</p>

# Main Routines
SAT-rx: semi-paired test in the random X setting. The inputs are kernel matrices for X, Y, and Z (covariates), the outputs are the p-values and test statistics. The detailed information about the inputs and outputs are provided in the code.

```
function [p_val0, p_val, p_valSemi, Sta, StaSemi] = HSIC_Test(Kx, Ky, Kzx, Kzy, pars)
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
```

SAT-fx: semi-paired test in the fixed X setting. The inputs are kernel matrices for X, Y, and Z (covariates), the outputs are the p-values and test statistics. The detailed information about the inputs and outputs are provided in the code.

``` function [p_val0, p_val, p_valSemi, Sta, StaSemi] = LMM_ScoreTest(Kx, Ky, Kzx, Kzy, pars)```

# Simulation Experiments
### Evaluation of type I error and test power of SAT-rx in the simulation setting (1).
```
Simulation (1), SAT-rx
├── run_jobs_simul_typeI_par.sh - script to evaluate type I error of SAT-rx.
├── run_jobs_simul_typeII_par.sh - Script to evaluate test power of SAT-rx.
```
### Evaluation of type I error and test power of SAT-fx in the simulation setting (1).
```
Simulation (1), SAT-rx
├── run_jobs_simul_lmm_typeI_par.sh - script to evaluate type I error of SAT-rx.
├── run_jobs_simul_lmm_typeII_par.sh - Script to evaluate test power of SAT-rx.
```

### Evaluation of type I error and test power of SAT-fx the simulation setting (2). The information about the COPD data can be found in [COPDGene](http://www.copdgene.org/).
```
Simulation (2), SAT-fx
├── run_jobs_copd_simul_typeI_par.sh - script to evaluate type I error of SAT-fx.
├── run_jobs_copd_simul_typeII_par.sh - Script to evaluate test power of SAT-fx.
```

# Real Experiments
### Phenotype discovery on General Population Cohort (GPC), Uganda. The genomic data have been deposited at the European Genome-phenome Archive ([EGA](https://www.ebi.ac.uk/ega/)) under accession number EGAS00001001558. Requests for access to phenotype data may be directed to data@apcdr.org.
```
Uganda data, SAT-rx
├── run_jobs_uganda_mimic.sh - script to calculate p-values of the paired phenotypes.
├── run_jobs_uganda_real.sh - Script to calculate p-values of unpaired phenotypes .
```

# Dependencies 
Our code depends on [KCI-test](http://people.tuebingen.mpg.de/kzhang/KCI-test.zip) and [MEGHA-v1](https://scholar.harvard.edu/tge/software/megha).

