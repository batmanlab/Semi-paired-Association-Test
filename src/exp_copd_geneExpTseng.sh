#!/bin/sh
module load matlab
echo 'locate and run'
echo "pheno_id=$1"
echo "gene_id=$2"
echo "exp_copd_geneExpTseng_func($1,$2)"
matlab -nodisplay -r "exp_copd_geneExpTseng_func($1,$2)"
