#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh
module load samtools

rawdir=/.mounts/labs/awadallalab/scratch/GCCR_WGS/merged_bams
basedir=/.mounts/labs/simpsonlab/users/hgibling/OHS/illumina
savedir=$basedir/PRDM9-bams

sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-sample-names.txt)

#samtools view -hb -F 3844 $rawdir/${sample}_01_LB01-01.markedDuplicates_sorted.bam chr5 > $rawdir/$sample.filtered-chr5-GRCh38.bam
#if [ $? -ne 0 ]; then fail "view failed"; fi
#samtools index $rawdir/$sample.filtered-chr5-GRCh38.bam
#if [ $? -ne 0 ]; then fail "index failed"; fi
cov=$(samtools coverage -Hr chr5 $rawdir/$sample.filtered-chr5-GRCh38.bam | cut -f7)
if [ $? -ne 0 ]; then fail "coverage failed"; fi
num=$(samtools stats $rawdir/$sample.filtered-chr5-GRCh38.bam | grep error | cut -f3)
if [ $? -ne 0 ]; then fail "stats failed"; fi
err=$(printf "%f" $num)
echo -e "$sample\t$cov\t$err" > $savedir/$sample-filtered-chr5-average-coverage.tsv
if [ $? -ne 0 ]; then fail "echo failed"; fi

echo "success!"
