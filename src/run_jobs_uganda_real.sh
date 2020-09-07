#!/bin/sh
jobName=exp_uganda_real
for jobId in $(seq 25 44)
do
    ${jobName}.sh $jobId 
done
squeue -u mgong1
