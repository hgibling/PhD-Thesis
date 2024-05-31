#!/bin/bash
#$ -cwd
#$ -t 1-52
source ~/scripts/error-handling.sh

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-genotyping
bamdir=$basedir/PRDM9-10k-bams
fqdir=$basedir/fqs

sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names.txt)
module load samtools/1.14

samtools view -hb -F 3844 $bamdir/$sample.illumina-prdm9-10k.bam > $bamdir/$sample.illumina-prdm9-10k-filtered.bam
if [ $? -ne 0 ]; then fail "failed view"; fi
samtools index $bamdir/$sample.illumina-prdm9-10k-filtered.bam
if [ $? -ne 0 ]; then fail "failed index"; fi
samtools fastq -1 $fqdir/$sample.illumina-prdm9-10k-filtered.r1.fq -2 $fqdir/$sample.illumina-prdm9-10k-filtered.r2.fq $bamdir/$sample.illumina-prdm9-10k-filtered.bam
if [ $? -ne 0 ]; then fail "failed fastq"; fi

echo "success!"
