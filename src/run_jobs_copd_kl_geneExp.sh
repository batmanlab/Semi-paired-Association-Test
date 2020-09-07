#!/bin/sh
jobName=exp_copd_kl_geneExp
#logDir=/pghbio/dbmi/batmanlab/mgong1/exp_folder/${jobName}
logDir=/pylon5/ac5616p/mgong1/exp_folder/${jobName}

if [ ! -d "${logDir}" ]; then
     mkdir ${logDir}
fi
outFile=${logDir}/output_${jobName}_${jobId}.log
errorFile=${logDir}/error_${jobName}_${jobId}.log
echo submit ${jobName}_${jobId}
echo ${outFile}
echo ${errorFile}
sbatch -A bi561ip -p DBMI -N 1 --ntasks-per-node=10 -o ${outFile} -e ${errorFile} -t 48:00:00  ${jobName}.sh 0
squeue -u mgong1
