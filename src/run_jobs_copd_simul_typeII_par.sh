#!/bin/sh
jobName=exp_copd_simul_typeII_par
#logDir=/pghbio/dbmi/batmanlab/mgong1/exp_folder/${jobName}
logDir=/pylon5/ac5616p/mgong1/exp_folder/${jobName}
h2=0.5
#for jobId in $(seq 1 1000)
##for jobId in 1000
#do
    ##for sigma_n in 0.8 0.3
    #for sigma_n in 0.8
    #do
        #echo "sigma_n=$sigma_n"
        #for dim in 50
        #do
            #if [ ! -d "${logDir}" ]; then
                #mkdir ${logDir}
            #fi
            #outFile=${logDir}/output_${jobName}_${h2}_${sigma_n}_${dim}_${jobId}.log
            #errorFile=${logDir}/error_${jobName}_${h2}_${sigma_n}_${dim}_${jobId}.log
            #echo submit ${jobName}_${h2}_${sigma_n}_${dim}_${jobId}
            #echo ${outFile}
            #echo ${errorFile}
            #sbatch -A bi561ip -p DBMI -N 1 --ntasks-per-node=10 -o ${outFile} -e ${errorFile} -t 48:00:00  ${jobName}.sh $h2 $sigma_n $dim $jobId
        #done
    #done
#done

for jobId in $(seq 1 999)
#for jobId in 1000
do
    #for sigma_n in 0.4 0.15
    for sigma_n in 0.4
    do
        echo "sigma_n=$sigma_n"
        for dim in 100
        do
            if [ ! -d "${logDir}" ]; then
                mkdir ${logDir}
            fi
            outFile=${logDir}/output_${jobName}_${h2}_${sigma_n}_${dim}_${jobId}.log
            errorFile=${logDir}/error_${jobName}_${h2}_${sigma_n}_${dim}_${jobId}.log
            echo submit ${jobName}_${h2}_${sigma_n}_${dim}_${jobId}
            echo ${outFile}
            echo ${errorFile}
            sbatch -A bi561ip -p DBMI -N 1 --ntasks-per-node=30 -o ${outFile} -e ${errorFile} -t 48:00:00  ${jobName}.sh $h2 $sigma_n $dim $jobId
        done
    done
done

#for jobId in $(seq 1 1000)
##for jobId in 1000
#do
    ##for sigma_n in 2 0.75
    #for sigma_n in 2
    #do
        #echo "sigma_n=$sigma_n"
        #for dim in 20
        #do
            #if [ ! -d "${logDir}" ]; then
                #mkdir ${logDir}
            #fi
            #outFile=${logDir}/output_${jobName}_${h2}_${sigma_n}_${dim}_${jobId}.log
            #errorFile=${logDir}/error_${jobName}_${h2}_${sigma_n}_${dim}_${jobId}.log
            #echo submit ${jobName}_${h2}_${sigma_n}_${dim}_${jobId}
            #echo ${outFile}
            #echo ${errorFile}
            #sbatch -A bi561ip -p DBMI -N 1 --ntasks-per-node=10 -o ${outFile} -e ${errorFile} -t 48:00:00  ${jobName}.sh $h2 $sigma_n $dim $jobId
        #done
    #done
#done
squeue -u mgong1
