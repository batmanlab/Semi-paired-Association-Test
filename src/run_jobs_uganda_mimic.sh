#!/bin/sh
jobName=exp_uganda_mimic
#logDir=/pghbio/dbmi/batmanlab/mgong1/exp_folder/${jobName}
logDir=/pylon5/ac5616p/mgong1/exp_folder/${jobName}
for jobId in $(seq 6 24)
do
    if [ ! -d "${logDir}" ]; then
         mkdir ${logDir}
    fi
    outFile=${logDir}/output_${jobName}_${jobId}.log
    errorFile=${logDir}/error_${jobName}_${jobId}.log
    echo submit ${jobName}_${jobId}
    echo ${outFile}
    echo ${errorFile}
    sbatch -A bi561ip -p DBMI -N 1 --ntasks-per-node=10 -o ${outFile} -e ${errorFile} -t 48:00:00  ${jobName}.sh $jobId 
    #sbatch -A bi561ip -p DBMI-GPU --gres=gpu:p100:1  -n 10 -o ${outFile} -e ${errorFile} -t 48:00:00  ${jobName}.sh
    #sbatch -A ac5616p -p GPU-shared --gres=gpu:p100:1  -o ${outFile} -e ${errorFile} -t 48:00:00  ${jobName}.sh
done
squeue -u mgong1
