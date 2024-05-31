#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh
bam=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/OHS/illumina/OHS-bams.txt)
sample="${bam%%.*}"
module load samtools
cd /.mounts/labs/awadallalab/scratch/GCCR_WGS/merged_bams
#samtools index $bam
#if [ $? -ne 0 ]; then fail "samtools index failed"; fi
samtools view -hb $bam chr5:23516673-23537764 > /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/OHS/illumina/PRDM9-bams/$sample-illumina-prdm9-10k-GRCh38.bam
if [ $? -ne 0 ]; then fail "samtools view failed"; fi
samtools index /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/OHS/illumina/PRDM9-bams/$sample-illumina-prdm9-10k-GRCh38.bam
if [ $? -ne 0 ]; then fail "samtools index failed"; fi

