#!/bin/bash
#$ -cwd
#$ -t 1-3
source ~/scripts/error-handling.sh

thesisdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36
triodir=$thesisdir/TRIO/redo-corrected
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src

sample=$(awk "NR==$SGE_TASK_ID" $triodir/all-samples-stats-corrected-300x.tsv | cut -f1)
cov=$(awk "NR==$SGE_TASK_ID" $triodir/all-samples-stats-corrected-300x.tsv | cut -f2)
err=$(awk "NR==$SGE_TASK_ID" $triodir/all-samples-stats-corrected-300x.tsv | cut -f3)

k=147
len=148

$gbkcdir/gbkc count -a $prdm9dir/all-alleles.fa -f $prdm9dir/flank10k.fa -1 $triodir/$sample-300x-GIAB-PRDM9-r1-aligncorrect.fa -2 $triodir/$sample-300x-GIAB-PRDM9-r2-aligncorrect.fa -d -K $k -m median -e $err -c $cov -l $len -o $triodir/$sample-300x-median-scores-corrected.tsv

if [ $? -ne 0 ]; then fail "gbkc count failed"; fi

echo "success!"
