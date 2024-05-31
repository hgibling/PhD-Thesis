#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh

prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
basedir=/.mounts/labs/simpsonlab/users/hgibling/OHS/illumina
bamdir=$basedir/PRDM9-bams
fqdir=$basedir/PRDM9-bams/filtered-fqs
savedir=$basedir/gbkc-distance-calls

sample=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-paired-average-coverage-error-frag.tsv | cut -f1)
cov=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-paired-average-coverage-error-frag.tsv | cut -f2)
err=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-paired-average-coverage-error-frag.tsv | cut -f3)
frag=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-paired-average-coverage-error-frag.tsv | cut -f4)
std=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-paired-average-coverage-error-frag.tsv | cut -f5)

len=151
k=131

time $gbkcdir/gbkc distance -a $prdm9dir/all-alleles-standard-flank10k.fa -1 $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_paired.r1.fq -2 $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_paired.r2.fq -d -K $k -e $err -c $cov -f $frag -s $std -l $len -m max -o $savedir/$sample.filtered-paired-max-PRDM910k-stats-flank.tsv
if [ $? -ne 0 ]; then fail "gbkc distance failed"; fi

echo "success!"
