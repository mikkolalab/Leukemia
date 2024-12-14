#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=48G,h_rt=12:00:00 
#$ -pe shared 1
#$ -r n # job is NOT rerunable
#$ -o SGE_STDOUT.txt 



ja=fastq.job 

PARMS=($(awk "NR==$SGE_TASK_ID" $ja))

FQO=${PARMS[0]}
FQT=${PARMS[1]}

# change for every dataset. `L008` may be different
SM=$(basename $FQT "_L008_R2_001.fastq.gz")
PREFIX=data/shortread_10X_5UTR_singlecell_fetaliver/bam/$SM


. /u/local/Modules/default/init/modules.sh
module load star
module load samtools


index=/u/project/gxxiao/gxxiao2/apps/genomes/hg38/STARIndex
GTF_FILE=/u/project/gxxiao/gxxiao2/apps/genomes/hg38/Gencode/gencode.v37.annotation.gtf #sorted.gtf.gz


log=log/star_remap.$SGE_TASK_ID


date > $log
echo ${SM} ${PREFIX} >> $log



## 1st alignment. (convert bam file to fastqs, then re-align)

STAR --runThreadN 4 --genomeDir $index \
	--readFilesIn ${FQT} \
	--readFilesCommand zcat \
	--outFilterMismatchNmax 3 \
	--outFilterType BySJout \
	--outFilterMultimapNmax 3 \
        --outFilterMismatchNoverReadLmax 0.05 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 8 \
	--alignIntronMin 20 \
	--alignIntronMax 100000 \
	--alignEndsType EndToEnd \
	--outFileNamePrefix ${PREFIX} \
	--outSAMattributes All \
	--outSAMtype BAM Unsorted \
	--outSAMattrRGline ID:${SM} &>> $log 


echo "LOG 1-pass done" >> $log 

## filter splice junctions sj.out

awk '$5 >= 1 && $5 <= 3 && $7 >= 3 {print}' ${PREFIX}SJ.out.tab > ${PREFIX}filtered.SJ.out.tab

## 2nd alignment. (convert bam file to fastqs, then re-align)

STAR --runThreadN 4 --genomeDir $index \
        --readFilesIn ${FQT} \
	--readFilesCommand zcat \
        --outFilterMismatchNmax 3 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 3 \
        --outFilterMismatchNoverReadLmax 0.05 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 8 \
        --alignIntronMin 20 \
        --alignIntronMax 500000 \
	--sjdbFileChrStartEnd ${PREFIX}filtered.SJ.out.tab \
	--sjdbGTFfile $GTF_FILE \
        --outFileNamePrefix ${PREFIX}.2pass \
        --outSAMattributes All \
        --outSAMtype BAM Unsorted \
        --outSAMattrRGline ID:${SM} &>> $log

echo "LOG 2-pass done" >> $log 

samtools sort --threads 16 ${PREFIX}.2passAligned.out.bam -o ${PREFIX}.2passAligned.out.sorted.bam &>> $log
samtools index ${PREFIX}.2passAligned.out.sorted.bam


echo "returns $?" >> $log 




