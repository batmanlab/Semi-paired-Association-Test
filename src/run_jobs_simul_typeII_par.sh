#!/bin/sh
jobName=exp_simul1_typeII_par
h2=0.5
for jobId in 1000
do
    for sigma_n in 0.8
    do
        echo "sigma_n=$sigma_n"
        for dim in 50
        do
            ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done

for jobId in 1000
do
    for sigma_n in 0.4
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
    for sigma_n in 0.2
    do
        echo "sigma_n=$sigma_n"
        for dim in 200
        do
            ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done

for jobId in $(seq 1 999)
do
    for sigma_n in 0.8
    do
        echo "sigma_n=$sigma_n"
        for dim in 50
        do
            ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done

for jobId in $(seq 1 999)
do
    for sigma_n in 0.4
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
    for sigma_n in 0.2
    do
        echo "sigma_n=$sigma_n"
        for dim in 200
        do
            ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done
