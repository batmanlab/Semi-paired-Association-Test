#!/bin/sh
module load matlab
echo 'locate and run'
echo "fev1=$1"
echo "exp_copd_lle_biomarker_func($1)"
matlab -nodisplay -r "exp_copd_lle_biomarker_func($1)"
