#!/bin/sh
module load matlab
echo 'locate and run'
echo "pheno_id=$1"
echo "exp_uganda_mimic_func($1)"
matlab -nodisplay -r "exp_uganda_mimic_func($1)"
