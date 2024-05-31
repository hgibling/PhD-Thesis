#!/bin/bash
#$ -cwd
#$ -t 1-52
source ~/scripts/error-handling.sh

prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-genotyping
bamdir=$basedir/PRDM9-10k-bams
fqdir=$basedir/fqs
savedir=$basedir/gbkc-calls

sample=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-PRDM910k-average-coverage-error.tsv | cut -f1)
cov=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-PRDM910k-average-coverage-error.tsv | cut -f2)
err=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-PRDM910k-average-coverage-error.tsv | cut -f3)

if [[ $sample == HG005 ]]; then k=231; len=250; else k=131; len=151; fi

time $gbkcdir/gbkc count -a $prdm9dir/all-alleles-standard.fa -1 $fqdir/$sample.illumina-prdm9-10k-filtered.r1.fq -2 $fqdir/$sample.illumina-prdm9-10k-filtered.r2.fq -d -K $k -f $prdm9dir/flank10k.fa -m coverage -e $err -c $cov -l $len -o $savedir/$sample.filtered-coverage-PRDM910k-stats-flank.tsv
if [ $? -ne 0 ]; then fail "gbkc count failed"; fi

echo "success!"
