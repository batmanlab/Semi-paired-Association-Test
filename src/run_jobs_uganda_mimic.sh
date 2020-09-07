#!/bin/sh
jobName=exp_uganda_mimic
for jobId in $(seq 6 24)
do
    ${jobName}.sh $jobId 
done
