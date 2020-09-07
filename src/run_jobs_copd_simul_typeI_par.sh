#!/bin/sh
jobName=exp_copd_simul_typeI_par
h2=0
for jobId in $(seq 1 1000)
#for jobId in 1000
do
    for sigma_n in 0.1 
    do
        echo "sigma_n=$sigma_n"
        for dim in 20 50 100
        #for dim in 50
        do
            ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done
