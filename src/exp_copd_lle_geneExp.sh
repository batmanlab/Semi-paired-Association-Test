#!/bin/sh
module load matlab
echo 'locate and run'
echo "pheno_id=$1"
echo "exp_copd_lle_geneExp_func($1)"
matlab -nodisplay -r "exp_copd_lle_geneExp_func($1)"
