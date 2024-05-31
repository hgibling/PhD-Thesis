#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh

prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
basedir=/.mounts/labs/simpsonlab/users/hgibling/OHS/illumina
bamdir=$basedir/PRDM9-bams
fqdir=$basedir/PRDM9-bams/filtered-fqs
savedir=$basedir/gbkc-count-calls

sample=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-paired-average-coverage-error-frag.tsv | cut -f1)
cov=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-paired-average-coverage-error-frag.tsv | cut -f2)
err=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-paired-average-coverage-error-frag.tsv | cut -f3)

len=151
k=131

time $gbkcdir/gbkc count -a $prdm9dir/all-alleles-standard.fa -1 $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_paired.r1.fq -2 $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_paired.r2.fq -d -K $k -f $prdm9dir/flank10k.fa -m coverage -e $err -c $cov -l $len -o $savedir/$sample.filtered-paired-coverage-PRDM910k-stats-flank.tsv
if [ $? -ne 0 ]; then fail "gbkc count failed"; fi

echo "success!"
