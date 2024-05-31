#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh
module load samtools

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS
bamdir=$basedir/mapped-GRCh38-filtered-maxima

sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-sample-names.txt)

cov=$(samtools coverage -Hr 5:23516673-23537764 $bamdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.bam | cut -f7)
num=$(samtools stats -t $basedir/PRDM9-10k-region-5.tsv $bamdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.bam| grep error | cut -f3)
err=$(printf "%f" $num)

echo -e "$sample\t$cov\t$err" > $bamdir/$sample-filtered-maxima-average-coverage-error.tsv
if [ $? -ne 0 ]; then fail "samtools failed"; fi

echo "success!"
