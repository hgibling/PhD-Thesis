#!/bin/bash
#$ -cwd
#$ -t 1-3
source ~/scripts/error-handling.sh

thesisdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36
triodir=$thesisdir/TRIO/redo
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src

sample=$(awk "NR==$SGE_TASK_ID" $triodir/all-samples-stats-2x250.tsv | cut -f1)
cov=$(awk "NR==$SGE_TASK_ID" $triodir/all-samples-stats-2x250.tsv | cut -f2)
err=$(awk "NR==$SGE_TASK_ID" $triodir/all-samples-stats-2x250.tsv | cut -f3)

k=247
len=250

$gbkcdir/gbkc count -a $prdm9dir/all-alleles.fa -f $prdm9dir/flank10k.fa -1 $triodir/$sample-2x250-GIAB-PRDM9-r1.fq -2 $triodir/$sample-2x250-GIAB-PRDM9-r2.fq -d -K $k -m median -e $err -c $cov -l $len -o $triodir/$sample-2x250-median-scores.tsv

if [ $? -ne 0 ]; then fail "gbkc count failed"; fi

echo "success!"
