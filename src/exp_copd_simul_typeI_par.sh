#!/bin/sh
# module load matlab
echo 'locate and run'
echo "h2=$1"
echo "sigma_n=$2"
echo "dim=$3"
echo "iter=$4"
echo "exp_copd_simul_typeI_par_func($1,$2,$3,$4)"
matlab -nodisplay -r "exp_copd_simul_typeI_par_func($1,$2,$3,$4)"
