#!/bin/bash
#$ -cwd
#$ -t 1-52
source ~/scripts/error-handling.sh

prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-genotyping
bamdir=$basedir/PRDM9-10k-bams
fqdir=$basedir/fqs
savedir=$basedir/gbkc-distance-calls

sample=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-PRDM910k-average-coverage-error-frag.tsv | cut -f1)
cov=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-PRDM910k-average-coverage-error-frag.tsv | cut -f2)
err=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-PRDM910k-average-coverage-error-frag.tsv | cut -f3)
frag=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-PRDM910k-average-coverage-error-frag.tsv | cut -f4)
std=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-PRDM910k-average-coverage-error-frag.tsv | cut -f5)

if [[ $sample == HG005 ]]; then k=231; len=250; else k=131; len=151; fi

time $gbkcdir/gbkc distance -a $prdm9dir/all-alleles-standard-flank10k.fa -1 $fqdir/$sample.illumina-prdm9-10k-filtered-paired.r1.fq -2 $fqdir/$sample.illumina-prdm9-10k-filtered-paired.r2.fq -d -K $k -e $err -c $cov -f $frag -s $std -l $len -m max -o $savedir/$sample.filtered-paired-max-PRDM910k-stats-flank.tsv
if [ $? -ne 0 ]; then fail "gbkc count failed"; fi

echo "success!"
