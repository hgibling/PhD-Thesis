#!/bin/bash
#$ -cwd
#$ -t 1-39
source ~/scripts/error-handling.sh

listdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject
savedir=/.mounts/labs/ont/scratch/HPRC_human_data

sample=$(awk "NR==$SGE_TASK_ID" $listdir/illumina-pacbio-prealigned-links-aws.tsv | cut -f1)
module load samtools

cd $savedir/$sample/illumina
samtools fastq -1 $sample-illumina-GRCh38-prdm9-10k.r1.fastq -2 $sample-illumina-GRCh38-prdm9-10k.r2.fastq $sample-illumina-GRCh38-prdm9-10k.bam
if [ $? -ne 0 ]; then fail "failed samtools fastq"; fi

echo "success!"

