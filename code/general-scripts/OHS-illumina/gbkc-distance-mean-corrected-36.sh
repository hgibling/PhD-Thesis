#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh

prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
basedir=/.mounts/labs/simpsonlab/users/hgibling/OHS/illumina
bamdir=$basedir/corrected-realigned-bams
fqdir=$basedir/corrected-reads
savedir=$basedir/gbkc-distance-calls-corrected-36

sample=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-corrected-average-coverage-error-frag.tsv | cut -f1)
cov=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-corrected-average-coverage-error-frag.tsv | cut -f2)
err=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-corrected-average-coverage-error-frag.tsv | cut -f3)
frag=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-corrected-average-coverage-error-frag.tsv | cut -f4)
std=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-corrected-average-coverage-error-frag.tsv | cut -f5)

k=131
len=151

time $gbkcdir/gbkc distance -a $prdm9dir/all-alleles-flank10k.fa -1 $fqdir/$sample.filtered-illumina-prdm9-10k-paired-corrected.r1.fa -2 $fqdir/$sample.filtered-illumina-prdm9-10k-paired-corrected.r2.fa -d -K $k -e $err -c $cov -f $frag -s $std -l $len -m mean -o $savedir/$sample.filtered-corrected-paired-mean-PRDM910k-stats-flank-corrected.tsv
if [ $? -ne 0 ]; then fail "gbkc count failed"; fi

echo "success!"
