#!/bin/sh
jobName=exp_copd_simul_typeII_par
h2=0.5

for jobId in 1000
do
    #for sigma_n in 0.8 0.3
    for sigma_n in 1
    do
        echo "sigma_n=$sigma_n"
        for dim in 60
        do
            ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done

for jobId in 1000
do
    #for sigma_n in 0.4 0.15
    for sigma_n in 1
    do
        echo "sigma_n=$sigma_n"
        for dim in 100
        do
           ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done

for jobId in 1000
do
    #for sigma_n in 2 0.75
    for sigma_n in 1
    do
        echo "sigma_n=$sigma_n"
        for dim in 20
        do
            ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done

for jobId in $(seq 1 999)
do
    #for sigma_n in 0.8 0.3
    for sigma_n in 1
    do
        echo "sigma_n=$sigma_n"
        for dim in 60
        do
            ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done

for jobId in $(seq 1 999)
do
    #for sigma_n in 0.4 0.15
    for sigma_n in 1
    do
        echo "sigma_n=$sigma_n"
        for dim in 100
        do
           ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done

for jobId in $(seq 1 999)
do
    #for sigma_n in 2 0.75
    for sigma_n in 1
    do
        echo "sigma_n=$sigma_n"
        for dim in 20
        do
            ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done