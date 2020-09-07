#!/bin/sh
echo 'locate and run'
echo "pheno_id=$1"
echo "exp_uganda_real_func($1)"
matlab -nodisplay -r "exp_uganda_real_func($1)"
