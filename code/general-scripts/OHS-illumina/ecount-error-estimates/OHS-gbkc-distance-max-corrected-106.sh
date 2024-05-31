#!/bin/bash
#$ -cwd
#$ -t 1-52
source ~/scripts/error-handling.sh

prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
paperdir=/.mounts/labs/simpsonlab/users/hgibling/PAPER
savedir=$paperdir/OHS/illumina
errordir=$savedir/error-estimates
fqdir=$savedir/fastqs-corrected
calldir=$savedir/rerun-gbkc/distance-calls-corrected-106

sample=$(awk "NR==$SGE_TASK_ID" $errordir/OHS-stats-error-corrected.tsv | cut -f1)
cov=$(awk "NR==$SGE_TASK_ID" $errordir/OHS-stats-error-corrected.tsv | cut -f2)
frag=$(awk "NR==$SGE_TASK_ID" $errordir/OHS-stats-error-corrected.tsv | cut -f3)
sdev=$(awk "NR==$SGE_TASK_ID" $errordir/OHS-stats-error-corrected.tsv | cut -f4)
err=$(awk "NR==$SGE_TASK_ID" $errordir/OHS-stats-error-corrected.tsv | cut -f5)

len=151
k=131

time $gbkcdir/gbkc distance -a $prdm9dir/all-alleles-standard.fa -1 $fqdir/$sample.filtered-illumina-prdm9-10k-paired-corrected.r1.fa -2 $fqdir/$sample.filtered-illumina-prdm9-10k-paired-corrected.r2.fa -d -K $k -f $prdm9dir/flank10k.fa -m max -e $err -c $cov -f $frag -s $sdev -l $len -o $calldir/$sample.filtered-paired-max-ecount-corrected-106.tsv
if [ $? -ne 0 ]; then fail "gbkc distance failed"; fi

echo "success!"
