#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh
module load samtools 

prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS
bamdir=$basedir/mapped-GRCh38-filtered-maxima
fqdir=$bamdir/fqs
savedir=$basedir/gbkc-count-calls

sample=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-maxima-nosoft-average-coverage-error-length.tsv | cut -f1)
cov=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-maxima-nosoft-average-coverage-error-length.tsv | cut -f2)
err=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-maxima-nosoft-average-coverage-error-length.tsv | cut -f3)
minlength=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-maxima-nosoft-average-coverage-error-length.tsv | cut -f4)
modelength=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-filtered-maxima-nosoft-average-coverage-error-length.tsv | cut -f5)

samtools fastq $bamdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima-nosoft.bam > $fqdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima-nosoft.fq
if [ $? -ne 0 ]; then fail "samtools astq failed"; fi

time $gbkcdir/gbkc count -a $prdm9dir/all-alleles-standard.fa -1 $fqdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima-nosoft.fq -d -k 21 -K 661 -i 20 -f $prdm9dir/flank10k.fa -m coverage -e $err -c $cov -l $modelength -N 10 -o $savedir/$sample.filtered-maxima-nosoft-coverage-PRDM910k-stats-flank.tsv
if [ $? -ne 0 ]; then fail "gbkc count failed"; fi

echo "success!"
