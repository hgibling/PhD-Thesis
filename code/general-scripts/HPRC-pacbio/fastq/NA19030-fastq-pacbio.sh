#!/bin/bash
#$ -cwd
#$ -t 1-100
source ~/scripts/error-handling.sh
bam=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/ont/scratch/HPRC_human_data/NA19030/pacbio/NA19030-list.txt)
module load samtools
cd /.mounts/labs/ont/scratch/HPRC_human_data/NA19030/pacbio

echo $bam

samtools fastq -0 $bam-temp.fq $bam
if [ $? -ne 0 ]; then fail "failed samtools fastq"; fi
