#!/bin/bash
#$ -cwd
#$ -t 1-52
source ~/scripts/error-handling.sh

prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
paperdir=/.mounts/labs/simpsonlab/users/hgibling/PAPER
savedir=$paperdir/OHS/illumina
errordir=$savedir/error-estimates
fqdir=$savedir/fastqs
calldir=$savedir/rerun-gbkc/count-calls-uncorrected-106

sample=$(awk "NR==$SGE_TASK_ID" $errordir/OHS-stats-error-uncorrected.tsv | cut -f1)
cov=$(awk "NR==$SGE_TASK_ID" $errordir/OHS-stats-error-uncorrected.tsv | cut -f2)
frag=$(awk "NR==$SGE_TASK_ID" $errordir/OHS-stats-error-uncorrected.tsv | cut -f3)
sdev=$(awk "NR==$SGE_TASK_ID" $errordir/OHS-stats-error-uncorrected.tsv | cut -f4)
err=$(awk "NR==$SGE_TASK_ID" $errordir/OHS-stats-error-uncorrected.tsv | cut -f5)

len=151
k=131

time $gbkcdir/gbkc count -a $prdm9dir/all-alleles-standard.fa -1 $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_paired.r1.fq -2 $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_paired.r2.fq -d -K $k -f $prdm9dir/flank10k.fa -m mean -e $err -c $cov -l $len -o $calldir/$sample.filtered-paired-mean-ecount-uncorrected-106.tsv
if [ $? -ne 0 ]; then fail "gbkc count failed"; fi

echo "success!"
