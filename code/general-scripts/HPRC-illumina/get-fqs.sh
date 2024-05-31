#!/bin/bash
#$ -cwd
#$ -t 1-52
source ~/scripts/error-handling.sh

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-genotyping
bamdir=$basedir/PRDM9-10k-bams
savedir=$basedir/fqs

sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names.txt)
module load samtools/1.14

samtools fastq -1 $savedir/$sample.illumina-prdm9-10k.r1.fq.gz -2 $savedir/$sample.illumina-prdm9-10k.r2.fq.gz $bamdir/$sample.illumina-prdm9-10k.bam
if [ $? -ne 0 ]; then fail "failed fastq"; fi

echo "success!"
