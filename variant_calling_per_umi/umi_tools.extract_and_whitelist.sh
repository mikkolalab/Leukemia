#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=24G,h_rt=23:00:00 
#$ -pe shared 1
#$ -r n # job is NOT rerunable
#$ -o SGE_STDOUT.txt 


ja=fastq.job 


PARMS=($(awk "NR==$SGE_TASK_ID" $ja))
FQO=${PARMS[0]}
FQT=$(echo $FQO | sed 's/_R1_/_R2_/g')



log=log/umi_tools.${SGE_TASK_ID}.txt

date > $log


umi_tools whitelist \
      --stdin ${FQO} \
      --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
      --plot-prefix=${FQO}.whitelist \
      --subset-reads=200000000 \
      --log2stderr > ${FQO}.whitelist.tsv  2>> $log

echo "Umitools whitelist returns $?" >> $log 

umi_tools extract \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
    -I ${FQO} \
    -S ${FQO}.extracted.fastq.gz  \
    --read2-in=${FQT} \
    --read2-out=${FQT}.extracted.fastq.gz \
    --whitelist=${FQO}.whitelist.tsv &>> $log 

echo "Umitools extract returns $?" >> $log 

